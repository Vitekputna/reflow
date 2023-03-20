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
    void compute_wall_flux(double dt, variables& var, mesh const& msh,void(*flux)(variables&,mesh const&,parameters const&));
    void compute_exact_flux(variables& var);
    void compute_cell_res(std::vector<std::vector<double>>& res, variables& var, mesh const& msh);
    void apply_source_terms(std::vector<std::vector<double>>& res, variables& var, mesh const& msh);
    void chemical_reactions(double dt,std::vector<std::vector<double>>& res, variables& var, mesh const& msh);

    // reconstruction
    void reconstruct(variables& var, mesh const& msh);
    void reconstructed_flux(int i, std::vector<double>& flux, std::vector<double> W, std::vector<double> const& grad, double dx);
    void reconstructed_wave_speed(int i, std::vector<double>& a, std::vector<double> W, std::vector<double> const& grad, double dx);
    inline double minmod(double a, double b);
    inline double van_albada(double a, double b);
    inline double van_leer(double a, double b);

    // explicit fluxes
    void Lax_Friedrichs_flux(variables& var,mesh const& msh, parameters const& par);
    void Kurganov_Tadmore(variables& var, mesh const& msh, parameters const& par);
    void HLL_flux(variables& var, mesh const& msh, parameters const& par);

    // exact flux function
    inline void Euler_flux(int i, std::vector<double>& flux, std::vector<double> const& W);

    // time integration
    void Explicit_Euler(variables& var,std::vector<std::vector<double>>& res, double dt);
    void Adams_Bashforth(variables& var,std::vector<std::vector<double>>& res_new,std::vector<std::vector<double>>& res_old, double dt);

    // runtime functions
    double time_step(variables const& var, mesh const& msh, double CFL);
    double max_residual(std::vector<std::vector<double>> const& res, variables const& var, int res_idx);
};