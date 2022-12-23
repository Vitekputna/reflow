#include "solver.hpp"
#include "thermodynamics.hpp"

void solver::compute_wall_flux(variables& var, mesh const& msh,
                               void(*flux)(variables&))
{
    flux(var);
}

void solver::compute_exact_flux(variables& var)
{
    for(int i = 0; i < var.N; i++)
    {
        Euler_flux(var.exact_flux[i],var.W[i]);
    }
}

void solver::compute_cell_res(std::vector<std::vector<double>>& res,variables& var, mesh const& msh)
{
    for(int i = 1; i < var.N-1; i++)
    {
        for(int k = 0; k < var.N_var; k++)
        {
            res[i][k] = -1/(msh.xf[i] - msh.xf[i-1])*(var.flux[i][k] - var.flux[i-1][k]);
        }
    }
}

void solver::Lax_Friedrichs_flux(variables& var)
{
    compute_exact_flux(var);

    for(int i = 0; i < var.N_walls; i++)
    {
        for(int k = 0; k < var.N_var; k++)
        {
            var.flux[i][k] = 0.5*(var.exact_flux[i+1][k] + var.exact_flux[i][k]) - 0.5*(var.W[i+1][k] - var.W[i][k]);
        }
    }
}

inline void solver::Euler_flux(std::vector<double>& flux, std::vector<double> const& W)
{
    double p = thermo::pressure(W,1.4);

    flux[0] = W[1];
    flux[1] = W[1]*W[1]/W[0] + p;
    flux[2] = (W[2] + p)*W[1]/W[0];
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