#include "lagrange_solver.hpp"
#include "particle.hpp"
#include "thermodynamics.hpp"
#include <vector>
#include <iostream>

extern double kappa, r;

double lagrange_solver::integrate_particle(double dt, particle& P, std::vector<double>& W, std::vector<double>& res)
{
    double C = 1e-8; // konstnta hybnosti
    double alfa = 100; // konstanta teploty

    double a = C/P.m*(W[1]/W[0] - P.u)*(W[1]/W[0] - P.u)*(W[1]/W[0] - P.u)/(std::abs((W[1]/W[0] - P.u)));

    double du = W[1]/W[0] - P.u;

    // dt = std::min(dt, std::abs(0.5*du/a));

    // Semi implicit RK4
    double K1 = dt*a;
    double up = P.u + K1/2;
    double ap = C/P.m*(W[1]/W[0] - up)*(W[1]/W[0] - up)*(W[1]/W[0] - up)/(std::abs((W[1]/W[0] - up)));
    double K2 = dt*ap;
    up = P.u + K2/2;
    ap = C/P.m*(W[1]/W[0] - up)*(W[1]/W[0] - up)*(W[1]/W[0] - up)/(std::abs((W[1]/W[0] - up)));
    double K3 = dt*ap;
    up = P.u + K3;
    ap = C/P.m*(W[1]/W[0] - up)*(W[1]/W[0] - up)*(W[1]/W[0] - up)/(std::abs((W[1]/W[0] - up)));
    double K4 = dt*ap;

    P.u += K1/6+K2/3+K3/3+K4/6;
    P.x += P.u*dt;

    P.T += alfa*dt*(thermo::temperature(W,kappa,r) - P.T);

    return dt;
}

void lagrange_solver::update_particles(double dt, std::vector<particle>& particles, variables& var, mesh const& msh, std::vector<std::vector<double>>& res)
{
    double t = 0;

    for(auto& P : particles)
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

