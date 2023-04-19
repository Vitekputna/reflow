#pragma once
#include <vector>
#include <random>
#include "mesh.hpp"
#include "variables.hpp"
#include "solver.hpp"
#include "thermodynamics.hpp"
#include "particle.hpp"
#include "specie.hpp"
#include "chem_solver.hpp"

typedef void(*Boundary_func)(variables&,mesh&,std::vector<double>&);

class reflow
{
    // Objects
    public:
    variables var;
    mesh msh;
    particle_manager par_man;
    thermo thermo_manager;
    chem_solver chemistry;

    // Boundary 
    std::vector<Boundary_func> boundary_func_vec;
    std::vector<std::vector<double>> boundary_values_vec;

    // Constatnts
    int n_dt;
    int n_res;
    int n_exp;

    int N, N_var;
    int n_comp = 0;
    int n_drop_frac = 0;
    int n_drop_mom = 0;

    double from, to;

    bool run_w_particles = false;

    double max_res;
    double t_end;
    double CFL;
    
    // Constructors
    reflow();
    reflow(variables& _var, mesh& _msh);
    reflow(int N, int N_var);
    reflow(int N, int N_var, std::vector<double> const& init);
    reflow(int N, int N_var, double from, double to);
    reflow(int N, int N_var, double from, double to, std::vector<double> const& init);
    reflow(int N, int N_var, std::vector<std::vector<double>> const& init);

    void variable_init();

    // Sources
    void apply_heat_source(double Q, double x_from, double x_to);
    void apply_mass_source(double M, double T, double x_from, double x_to, std::vector<double> const& comp);

    // Geometry
    void spline_geometry(std::vector<std::vector<std::vector<double>>> curves, int n);
    void bump_geometry(); 
    void refine_mesh(std::vector<std::vector<double>> ref);

    // Initial conditions
    void initial_conditions(std::vector<double> const& init);
    void initial_conditions(std::vector<std::vector<double>> const& init);
    void initial_conditions(int N_drop, int N_drop_mom, std::vector<double> const& init);
    void initial_conditions(int N_drop, int N_drop_mom, std::vector<std::vector<double>> const& init);

    // Thermodynamics
    void add_specie(double r, double kappa, double Mm, std::vector<double> cp_coeff);

    // Chemistry
    void add_reaction(reaction& R);

    // Lagrangian
    void add_lagrangian_mono_particles(double specie_idx, double mass_flux, double rho, double r, double x,
                                       double u, double T, double T_boil, double vap_heat, double C);

    void add_lagrangian_unif_particles(double specie_idx, double mass_flux, double rho, double r_from, double r_to, 
                                       double x, double u, double T, double T_boil, double vap_heat, double C);

    void apply_lagrangian_particle_inlet(double dt);

    // Boundary
    void add_boundary_function(Boundary_func,std::vector<double> values);
    void apply_boundary_conditions(); 

    // Export
    void export_particles(std::vector<particle>& particles);
    void init_particles(int N_max, int N_particles, int N_per_group);

    // run criteria
    bool maximum_time(double T, double res);
    bool maximum_res(double T, double res);
    bool maximum_time_res(double T, double res);

    void solve(double t_end, double max_residual, double CFL);
};