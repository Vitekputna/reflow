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

class reflow
{
    public:
    variables var;
    mesh msh;
    particle_manager par_man;
    thermo thermo_manager;
    chem_solver chemistry;

    void(*left_boundary)(variables&,mesh&,std::vector<double>&);
    void(*right_boundary)(variables&,mesh&,std::vector<double>&);

    std::vector<double> left_values;
    std::vector<double> right_values;

    int n_dt = 2;
    int n_res = 200;
    int n_exp = 5000;
    int N, N_var;
    double from, to;

    bool run_w_particles = false;

    volatile double max_res = 1;
    
    // Constructors
    reflow();
    reflow(variables& _var, mesh& _msh);
    reflow(int N, int N_var);
    reflow(int N, int N_var, std::vector<double> const& init);
    reflow(int N, int N_var, double from, double to);
    reflow(int N, int N_var, double from, double to, std::vector<double> const& init);
    reflow(int N, int N_var, std::vector<std::vector<double>> const& init);

    // Sources
    void apply_heat_source(double Q, double x_from, double x_to);
    void apply_mass_source(double M, double T, double x_from, double x_to, std::vector<double> const& comp);

    // Geometry
    void spline_geometry(std::vector<std::vector<std::vector<double>>> curves, int n);
    void bump_geometry(); 
    void refine_mesh(std::vector<std::vector<double>> ref);

    // Initial conditions
    void initial_conditions(std::vector<double> const& init);
    void initial_conditions(int N_drop, int N_drop_mom, std::vector<double> const& init);
    void particle_distribution();

    // Thermodynamics
    void add_specie(double r, double kappa, double Mm, std::vector<double> cp_coeff);

    // Chemistry
    void add_reaction(reaction& R);

    // Boundary
    void set_boundary(void(*left)(variables&,mesh&,std::vector<double>&), void(*right)(variables&,mesh&,std::vector<double>&));
    void set_boundary(void(*left)(variables&,mesh&,std::vector<double>&), std::vector<double> _left_values,
                      void(*right)(variables&,mesh&,std::vector<double>&), std::vector<double> _right_values);

    // Export
    void export_particles(std::vector<particle>& particles);
    void init_particles(int N_max, int N_particles, int N_per_group);

    void solve();
};