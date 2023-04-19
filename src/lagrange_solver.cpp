#include "lagrange_solver.hpp"
#include "particle.hpp"
#include "thermodynamics.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "omp.h"


// extern double kappa, r;

inline double lagrange_solver::acceleration(double Cd, double r, double du)
{
    return 3*Cd*std::abs(du)*du/2/r;
}

inline double lagrange_solver::heat_flux(double alfa, double C, double dT)
{
    return alfa*dT/C;
}

inline double lagrange_solver::radius_change(double D, double r, double rho, double dT)
{
    // return std::min(0.0,-D/4/3.14159/((r*r)*rho));   
    return std::min(0.0,-D*log(1+1e-3*dT)/r);
}

inline double lagrange_solver::mass_flux(double r, double dT)
{
    return std::min(0.0,3*r*log(1 + 1e-3*std::max(0.0,dT)));
}

double lagrange_solver::integrate_particle(double dt, double V, particle& P, std::vector<double>& W, std::vector<double>& res)
{
    static double K1,K2,K3,K4;
    static double ap,up,rp;
    static double Tf, uf;

    Tf = thermo::T[P.last_cell_idx];
    uf = W[W.size()-2]/W[0];

    double C = 1e-2;        // momentum transfer constant
    double D = 1e-5;        // mass transfer constant
    double alfa = 10;       // heat transfer constant

    // Velocity

    if(P.r > 2e-6)
    {
        K1 = dt*acceleration(C,P.r,uf - P.u);
        up = P.u + K1/2;
        ap = acceleration(C,P.r,uf - up);
        K2 = dt*ap;
        up = P.u + K2/2;
        ap = acceleration(C,P.r,uf - up);
        K3 = dt*ap;
        up = P.u + K3;
        ap = acceleration(C,P.r,uf - up);
        K4 = dt*ap;
        P.u += K1/6+K2/3+K3/3+K4/6;
    }
    else
    {
        P.u = uf;   
    }
    
    P.x += P.u*dt;

    if(P.x <= 0.005)
    {
        return 0.0;
    }

    P.T += alfa*dt*(Tf - P.T);

    // radius
    K1 = dt*radius_change(D,P.r,P.rho,Tf-200);
    rp = P.r + K1/2;
    K2 = dt*radius_change(D,rp,P.rho,Tf-200);
    rp = P.r + K2/2;
    K3 = dt*radius_change(D,rp,P.rho,Tf-200);
    rp = P.r + K3;
    K4 = dt*radius_change(D,rp,P.rho,Tf-200);
    P.r += K1/6+K2/3+K3/3+K4/6;

    // mass
    double m0 = P.M;
    double md;

    if(P.r < 0)
    {
        md = m0/dt/V;

        res[0] += md;
        res[2] += md;
        res[3] += md*P.u;
        res[4] += md*(thermo::enthalpy(Tf,std::vector<double>{0,0,1}) + pow(P.u,2)/2);

        P.reset();
        return md;
    }

    P.m = 4*M_PI*pow(P.r,3)*P.rho/3;

    P.M = P.N*P.m;

    md = (m0 - P.M)/dt/V;

    res[0] += md;
    res[2] += md;
    res[3] += md*P.u;
    res[4] += md*(thermo::enthalpy(Tf,std::vector<double>{0,0,1}) + pow(P.u,2)/2);

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
            // V = msh.A[particles[j].last_cell_idx]*(msh.xf[particles[j].last_cell_idx] - msh.xf[particles[j].last_cell_idx-1]);
            V = msh.V[particles[j].last_cell_idx];
            var.md[particles[j].last_cell_idx][2] += integrate_particle(dt,V,particles[j],var.W[particles[j].last_cell_idx],res[particles[j].last_cell_idx]);
        }
    }
}

