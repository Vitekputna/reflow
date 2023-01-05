#pragma once
#include "variables.hpp"
#include "mesh.hpp"
#include "thermodynamics.hpp"

struct parameters
{
    double dx, dt;
    double eps;

    parameters(double _dx, double _dt, double _eps) : dx{_dx}, dt{_dt}, eps{_eps} {}
};

namespace solver
{
    // explicit part
    void compute_wall_flux(double dt, variables& var, mesh const& msh,void(*flux)(variables&,parameters const&));
    void compute_exact_flux(variables& var);
    void compute_cell_res(std::vector<std::vector<double>>& res, variables& var, mesh const& msh);
    void apply_source_terms(std::vector<std::vector<double>>& res, variables& var, mesh const& msh);

    // explicit fluxes
    void Lax_Friedrichs_flux(variables& var, parameters const& par);

    // exact flux function
    inline void Euler_flux(std::vector<double>& flux, std::vector<double> const& W, double kappa);

    // time integration
    void Explicit_Euler(variables& var,std::vector<std::vector<double>>& res, double dt);
    void Adams_Bashforth(variables& var,std::vector<std::vector<double>>& res_new,std::vector<std::vector<double>>& res_old, double dt);

    // runtime functions
    double time_step(variables const& var, mesh const& msh, double kappa, double CFL);
    double max_residual(std::vector<std::vector<double>> const& res, variables const& var, int res_idx);
};