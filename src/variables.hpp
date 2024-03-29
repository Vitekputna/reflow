#pragma once
#include <vector>
#include <string>

#include "mesh.hpp"
#include "particle.hpp"

struct variables
{
    std::vector<std::vector<double>> W;                 // vector of variables
    std::vector<std::vector<double>> flux;              // vector of numeric fluxes
    std::vector<std::vector<double>> exact_flux;        // vector of exact fluxes
    std::vector<std::vector<double>> grad;              // vector of variable gradients
    std::vector<double> grad_p;

    std::vector<double> q;                              // vector of energy sources
    std::vector<std::vector<double>> md;                // vector of mass sources

    int N_var;          // number of variables
    int N;              // number of cells
    int N_walls;        // number of walls
    int N_ghosts = 2;   // number of ghost cells
    
    static int N_comp;         // number of components

    static int N_drop_eq;                       // total number of droplet eq.
    static int N_drop_frac;                     // number of drop. fractions
    static int N_drop_mom_eq;                   // number of drop. momenta eq.
    static int N_drop_eng_eq;                   // number of drop. energt eq.
    static std::vector<int> drop_mom_idx;       // momentum eq. idx for all droplet fractions
    static std::vector<int> quisc_drop_idx;     // quiscent droplet fractions idx
    static std::vector<int> active_drop_idx;    // active droplet fractions idx (drops that have mom. eq.)

    static int mom_idx;
    static int eng_idx;   // indices for fluid momentum and energy equation

    //default
    variables();
    variables(int _N_var, int _N);

    // wo droplets
    variables(int _N_var, int _N, std::vector<double> const& W_0);
    variables(int _N_var, int _N, std::vector<std::vector<double>> const& W_0);

    // w droplets
    void drop_init(int _N_drop_frac, bool droplet_momenta, bool droplet_energy);

    variables(int _N_var, int _N, int _N_drop_frac, bool droplet_momenta, std::vector<double> const& W_0);
    variables(int _N_var, int _N, int _N_drop_frac, bool droplet_momenta, std::vector<std::vector<double>> const& W_0);

    variables(int _N_var, int _N, int _N_drop_frac, bool droplet_momenta, bool droplet_energy, std::vector<double> const& W_0);
    variables(int _N_var, int _N, int _N_drop_frac, bool droplet_momenta, bool droplet_energy, std::vector<std::vector<double>> const& W_0);


    void apply_heat_source(double Q_tot, double x_from, double x_to, mesh const& msh);
    void apply_mass_source(double M_tot, double T, double x_from, double x_to, mesh const& msh, std::vector<double> const& comp);

    void export_to_file(std::string path, mesh const& msh, std::vector<particle> const& particles);
    void export_timestep(double t, mesh const& msh, std::vector<particle> const& particles);
};