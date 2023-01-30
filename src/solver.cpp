#include "solver.hpp"
#include "thermodynamics.hpp"
#include <cmath>

extern double kappa;

void solver::compute_wall_flux(double dt, variables& var, mesh const& msh,
                               void(*flux)(variables&,parameters const&))
{
    flux(var, parameters(msh.dx_min,dt,0.9));
}

void solver::compute_exact_flux(variables& var)
{
    for(int i = 0; i < var.N; i++)
    {
        Euler_flux(var.exact_flux[i],var.W[i],kappa);
    }
}

void solver::compute_cell_res(std::vector<std::vector<double>>& res,variables& var, mesh const& msh)
{
    for(int i = 1; i < var.N-1; i++)
    {
        for(int k = 0; k < var.N_var; k++)
        {
            res[i][k] = -1/(msh.A[i]*(msh.xf[i] - msh.xf[i-1]))*(msh.Af[i]*var.flux[i][k] - msh.Af[i-1]*var.flux[i-1][k]);
        }
    }
}

void solver::apply_source_terms(std::vector<std::vector<double>>& res, variables& var, mesh const& msh)
{
    for(int i = 1; i < var.N-1; i++)
    {
        // res[i][0] += var.md[i];
        res[i][var.mom_idx] += ((msh.Af[i]-msh.Af[i-1])/(msh.xf[i]-msh.xf[i-1]))/msh.A[i]*thermo::pressure(var.W[i],kappa);
        res[i][var.eng_idx] += var.q[i]/msh.A[i]; 
    }
}

void solver::Lax_Friedrichs_flux(variables& var, parameters const& par)
{
    compute_exact_flux(var);

    for(int i = 0; i < var.N_walls; i++)
    {
        for(int k = 0; k < var.N_var; k++)
        {
            var.flux[i][k] = 0.5*(var.exact_flux[i+1][k] + var.exact_flux[i][k]) - par.eps*0.5*par.dx/par.dt*(var.W[i+1][k] - var.W[i][k]);
        }
    }
}

inline void solver::Euler_flux(std::vector<double>& flux, std::vector<double> const& W, double kappa)
{
    double p = thermo::pressure(W,kappa);
    int n_var = flux.size()-2;

    flux[0] = W[n_var];

    for(auto idx = 1; idx < n_var; idx++)
    {
        flux[idx] = W[n_var]*W[idx]/W[0];
    }

    flux[n_var] = W[n_var]*W[n_var]/W[0] + p;
    flux[n_var+1] = (W[n_var+1] + p)*W[n_var]/W[0];
}

void solver::Explicit_Euler(variables& var, std::vector<std::vector<double>>& res, double dt)
{
    for(int i = 1; i < var.N-1; i++)
    {
        for(int k = 0; k < var.N_var; k++)
        {
            var.W[i][k] += dt*res[i][k];
        }
    }
}

void solver::Adams_Bashforth(variables& var,std::vector<std::vector<double>>& res_new,std::vector<std::vector<double>>& res_old, double dt)
{
    for(int i = 1; i < var.N-1; i++)
    {
        for(int k = 0; k < var.N_var; k++)
        {
            var.W[i][k] = var.W[i][k] + dt*(1.5*res_new[i][k] - 0.5*res_old[i][k]);
        }
    }
}

double solver::time_step(variables const& var, mesh const& msh, double kappa, double CFL)
{
    double dt = 1e10;

    for(auto const& w : var.W)
    {
        dt = std::min(dt, msh.dx_min/(std::abs(w[var.mom_idx]/w[0] + thermo::speed_of_sound(w,kappa))));
    }

    return CFL*dt;
}

double solver::max_residual(std::vector<std::vector<double>> const& res, variables const& var, int res_idx)
{
    double max_res = 0;

    for(unsigned int i = 1; i < res.size()-1; i++)
    {
        max_res = std::max(max_res, std::abs(res[i][res_idx]));
    }

    return max_res;
}