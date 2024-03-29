#include "lagrange_solver.hpp"
#include "particle.hpp"
#include "thermodynamics.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "omp.h"
#include "euler_droplets.hpp"


// extern double kappa, r;

inline double lagrange_solver::acceleration(double r, double rho_l, double rho_g, double mu, double du)
{
    const double Re = (rho_g*std::abs(du)*2*r)/mu;
    const double Cd = euler_droplets::Ingebo(Re);

    return (3*Cd*rho_g*std::abs(du)*du)/(8*r*rho_l);
}

inline double lagrange_solver::heat_flux(double r, double T_gas, double T_drop, double rho_g, double u_l,
                                         double u_g, double mu, double cp, double k, double h_vap, double phi)
{
    const double Re = (rho_g*std::abs(u_g - u_l)*2*r)/mu;
    const double Pr = cp*mu/k;
    const double BT = cp*(T_gas-T_drop)*phi/h_vap;
    const double Nu = euler_droplets::Nusselt_evaporation(Re,Pr,BT);

    return Nu*2*r*M_PI*k*(T_gas-T_drop);
}

inline double lagrange_solver::radius_change(double r, double Ys, double Yf, double rho_l, double rho_g, double u_l, double u_g, double mu, double D)
{
    const double Re = (rho_g*std::abs(u_g - u_l)*2*r)/mu;
    const double Sc = mu/(rho_g*D);
    const double BM = (Ys-Yf)/(1-Ys);
    const double Sh = euler_droplets::Sherwood_evaporation(Re,Sc,BM);

    return std::min(0.0,(Sh*rho_g*D)/(2*r*rho_l)*(Ys-Yf));
}

double lagrange_solver::integrate_particle(double dt, double V, particle& P, std::vector<double>& W, std::vector<double>& res)
{
    double K1,K2,K3,K4;
    double ap,up,rp,Tp;

    const double rho_l = thermo::species[2].rho_liq;
    const double rho_gas = thermo::density(W);

    const double Tf = thermo::T[P.last_cell_idx];
    const double p = thermo::p[P.last_cell_idx];
    const double u_gas = W[variables::mom_idx]/W[0];

    const double u_drop = P.u;
    const double T_drop = P.T;
    const double T_boil = thermo::boil_temp(2,p);

    const double r = P.r;

    std::vector<double> comp = std::vector<double>(variables::N_comp,0);
    thermo::composition(comp,W);

    const double Ys = comp[2];
    const double Yf = euler_droplets::fuel_mass_fraction(p,P.T);

    const double Df = thermo::difusivity(comp,Tf);
    const double mu = thermo::viscosity(comp,Tf);
    const double cp = thermo::cp_mix_comp(comp,Tf);
    const double k = thermo::thermal_conductivity(comp,Tf);
    const double h_vap = thermo::species[2].h_vap;
    const double C = thermo::species[2].C;

    const double phi = euler_droplets::heat_evap_interp(T_drop,T_boil);

    // Velocity
    if(P.r > 1e-6)
    {
        K1 = dt*acceleration(P.r,rho_l,rho_gas,mu,u_gas - P.u);
        up = P.u + K1/2;
        ap = acceleration(P.r,rho_l,rho_gas,mu,u_gas - up);
        K2 = dt*ap;
        up = P.u + K2/2;
        ap = acceleration(P.r,rho_l,rho_gas,mu,u_gas - up);
        K3 = dt*ap;
        up = P.u + K3;
        ap = acceleration(P.r,rho_l,rho_gas,mu,u_gas - up);
        K4 = dt*ap;
        P.u += K1/6+K2/3+K3/3+K4/6;
    }
    else
    {
        P.u = u_gas;   
    }

    P.x += P.u*dt;

    if(P.x <= 0.005)
    {
        return 0.0;
    }

    // radius change based on evaporation
    K1 = dt*radius_change(r,Ys,Yf,rho_l,rho_gas,u_drop,u_gas,mu,Df);
    rp = P.r + K1/2;    
    K2 = dt*radius_change(rp,Ys,Yf,rho_l,rho_gas,u_drop,u_gas,mu,Df);
    rp = P.r + K2/2;
    K3 = dt*radius_change(rp,Ys,Yf,rho_l,rho_gas,u_drop,u_gas,mu,Df);
    rp = P.r + K3;
    K4 = dt*radius_change(rp,Ys,Yf,rho_l,rho_gas,u_drop,u_gas,mu,Df);
    P.r += K1/6+K2/3+K3/3+K4/6;

    const double Q = heat_flux(r,Tf,T_drop,rho_gas,u_drop,u_gas,mu,cp,k,h_vap,phi);

    // radius boil correction
    P.r += -dt*(Q*phi)/(h_vap*4*M_PI*r*r*rho_l);

    const double Q_corrected = Q*(1-phi);

    // mass
    const double M0 = P.M;
    const double m0 = P.m;
    double md;
    
    if(P.r < 1e-6)
    {
        md = M0/dt/V;

        res[0] += md;
        res[2] += md;
        res[3] += md*P.u;
        res[4] += md*(thermo::enthalpy(P.T,std::vector<double>{0,0,1}) + pow(P.u,2)/2);

        P.reset();

        return md;
    }
    else
    {
        P.m = (4*M_PI*pow(P.r,3)*P.rho)/3;

        P.M = P.N*P.m;

        md = (M0 - P.M)/dt/V;

        res[0] += md;
        res[2] += md;
        res[3] += md*P.u;
        res[4] += md*(thermo::enthalpy(P.T,std::vector<double>{0,0,1}) + pow(P.u,2)/2);
    }

    //Temperature
    const double dm = m0-P.m;

    P.T += (dt*Q_corrected- dm*h_vap)/(P.m*C);
    res[4] += -P.N*Q/V;

    return md;
}

void lagrange_solver::update_particles(double dt, std::vector<particle>& particles, variables& var, mesh const& msh, std::vector<std::vector<double>>& res)
{
    double V;

    omp_set_num_threads(6);
    #pragma omp parallel for shared(dt, particles, var, msh, res) private(V)
    for(int j = 0; j < particles.size();j++)
    {
        if(particles[j].in_use)
        {
            // find where is particle located in mesh
            for(int i = particles[j].last_cell_idx; i < msh.N; i++)
            {
                if(particles[j].x > msh.xf[i-1] && particles[j].x <= msh.xf[i])
                {
                    particles[j].last_cell_idx = i;
                    break;
                }
            }

            if(particles[j].x > msh.xf.back())
            {
                particles[j].reset();
            }

            // update particle speed, mass, temp...
            V = msh.V[particles[j].last_cell_idx];
            var.md[particles[j].last_cell_idx][2] += integrate_particle(dt,V,particles[j],var.W[particles[j].last_cell_idx],res[particles[j].last_cell_idx]);
        }
    }
}

