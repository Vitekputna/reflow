#pragma once
#include "variables.hpp"
#include "mesh.hpp"
#include "thermodynamics.hpp"

namespace solver
{
    // explicit part
    void compute_wall_flux(variables& var, mesh const& msh,void(*flux)(variables&));
    void compute_exact_flux(variables& var);
    void compute_cell_res(std::vector<std::vector<double>>& res, variables& var, mesh const& msh);
    void apply_source_terms(variables& var, mesh const& msh);

    // explicit fluxes
    void Lax_Friedrichs_flux(variables& var);

    // exact flux function
    inline void Euler_flux(std::vector<double>& flux, std::vector<double> const& W);

    // time integration
    void Explicit_Euler(variables& var,std::vector<std::vector<double>>& res, double dt);
};