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
    void compute_wall_flux(double dt, variables& var, mesh const& msh,void(*flux)(variables&,mesh const&,parameters const&,const int, const int,std::vector<int> const&), const int from, const int to, std::vector<int> const& flux_idx);
    void compute_exact_flux(variables& var, const int from, const int to);
    void compute_cell_res(std::vector<std::vector<double>>& res, variables& var, mesh const& msh, const int from, const int to);
    void apply_source_terms(std::vector<std::vector<double>>& res, variables& var, mesh const& msh, const int from, const int to);
    void chemical_reactions(double dt,std::vector<std::vector<double>>& res, variables& var, mesh const& msh, const int from, const int to);

    // reconstruction
    void reconstruct(variables& var, mesh const& msh, const int from, const int to);
    void reconstruct_pressure(variables& var, mesh const& msh, const int from, const int to);
    void wall_pressures(const int wall_idx, std::vector<double>& pressures, variables const& var, mesh const& msh);
    void wall_states(const int wall_idx, std::vector<double>& left_state, std::vector<double>& right_state, variables const& var, mesh const& msh);
    void reconstructed_state(std::vector<double>& state, std::vector<double> const& grad, double dx);
    void reconstructed_flux(int i, std::vector<double>& flux, std::vector<double> W, std::vector<double> const& grad, double dx);
    void reconstructed_wave_speed(int i, std::vector<double>& a, std::vector<double> W, std::vector<double> const& grad, double dx);
    inline double minmod(double a, double b);
    inline double van_albada(double a, double b);
    inline double van_leer(double a, double b);

    // explicit fluxes
    void Lax_Friedrichs_flux(variables& var,mesh const& msh, parameters const& par, const int from, const int to, std::vector<int> const& var_idx);
    void Kurganov_Tadmore(variables& var, mesh const& msh, parameters const& par, const int from, const int to, std::vector<int> const& var_idx);
    void HLL_flux(variables& var, mesh const& msh, parameters const& par, const int from, const int to, std::vector<int> const& var_idx);
    void HLL2_flux(variables& var, mesh const& msh, parameters const& par, const int from, const int to, std::vector<int> const& var_idx);
    void AUSM_flux(variables& var, mesh const& msh, parameters const& par, const int from, const int to, std::vector<int> const& var_idx);
    void AUSM2_flux(variables& var, mesh const& msh, parameters const& par, const int from, const int to, std::vector<int> const& var_idx);
    void upwind(variables& var, mesh const& msh, parameters const& par, const int from, const int to, std::vector<int> const& var_idx);

    double AUSM_wall_mach_number(double M_left, double M_right);
    double AUSM_wall_pressure(double M_left, double M_right, double p_left, double p_right);

    // droplet transport
    void droplet_transport(std::vector<std::vector<double>>& res, variables& var, mesh const& msh, const int from, const int to);

    // exact flux function
    inline void Euler_flux(int i, std::vector<double>& flux, std::vector<double> const& W);
    inline void Euler_flux(const double p, std::vector<double>& flux, std::vector<double> const& W);

    // time integration
    void Explicit_Euler(variables& var,std::vector<std::vector<double>>& res, double dt, const int from, const int to);
    void Adams_Bashforth(variables& var,std::vector<std::vector<double>>& res_new,std::vector<std::vector<double>>& res_old, double dt);

    // runtime functions
    double time_step(variables const& var, mesh const& msh, double CFL);
    double max_residual(std::vector<std::vector<double>> const& res, variables const& var, int res_idx);
    double L2_norm(mesh const& msh, std::vector<std::vector<double>> const& res, int res_idx);
};