#include "lagrange_solver.hpp"
#include "particle.hpp"
#include "thermodynamics.hpp"
#include <vector>
#include <iostream>
#include <cmath>

extern double kappa, r;

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
    return std::min(0.0,-D*dT/4/3.14159/((r*r)*rho));   
}

double lagrange_solver::integrate_particle(double dt, particle& P, std::vector<double>& W, std::vector<double>& res)
{
    double K1,K2,K3,K4;
    double ap,up,rp;
    double Tf, uf;

    Tf = thermo::temperature(W,kappa,r);
    uf = W[1]/W[0];

    double C = 1e-2; // momentum transfer constant
    double D = 1e-6; // mass transfer constant
    double alfa = 100; // heat transfer constant

    // double du = W[1]/W[0] - P.u;
    // dt = std::min(dt, std::abs(0.5*du/a));


    // Velocity
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

    // Position
    P.x += P.u*dt;

    // Temperature
    P.T += alfa*dt*(thermo::temperature(W,kappa,r) - P.T);

    // radius
    K1 = dt*radius_change(D,P.r,P.rho,Tf-P.T);
    rp = P.r + K1/2;
    K2 = dt*radius_change(D,rp,P.rho,Tf-P.T);
    rp = P.r + K2/2;
    K3 = dt*radius_change(D,rp,P.rho,Tf-P.T);
    rp = P.r + K3;
    K4 = dt*radius_change(D,rp,P.rho,Tf-P.T);
    P.r += K1/6+K2/3+K3/3+K4/6;

    // mass
    double m0 = P.M;
    P.m = 4/3*3.14159*P.r*P.r*P.r*P.rho;
    P.M = P.N*P.m;

    if(P.r < 0) P.reset();

    // res[2] += (m0 - P.M)*1e10/dt;

    return dt;
}

void lagrange_solver::update_particles(double dt, std::vector<particle>& particles, variables& var, mesh const& msh, std::vector<std::vector<double>>& res)
{
    // double t = 0;

    for(auto& P : particles)
    {
        if(P.in_use)
        {
            // find where is particle located in mesh
            for(int i = P.last_cell_idx; i < msh.N; i++)
            {
                if(P.x < msh.xf[0]) break;
                if(P.x > msh.xf[i-1] && P.x <= msh.xf[i])
                {
                    P.last_cell_idx = i;
                    break;
                }
            }

            if(P.x > msh.xf.back())
            {
                P.reset();
            }

            // update particle speed, mass, temp...
            integrate_particle(dt,P,var.W[P.last_cell_idx],res[P.last_cell_idx]);
        }
    }
}

