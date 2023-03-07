#pragma once
#include <vector>
#include <string>

#include "mesh.hpp"
#include "particle.hpp"

struct variables
{
    std::vector<std::vector<double>> W;
    std::vector<std::vector<double>> flux;
    std::vector<std::vector<double>> exact_flux;
    std::vector<std::vector<double>> grad;

    std::vector<double> q;
    std::vector<std::vector<double>> md;

    int N_var; // number of variables
    int N; // number of cells
    int N_walls; // number of walls
    int N_ghosts = 2; // number of ghost cells
    int N_comp; // number of components
    int mom_idx, eng_idx;

    variables();
    variables(int _N_var, int _N);
    variables(int _N_var, int _N, std::vector<double> const& W_0);
    variables(int _N_var, int _N, std::vector<std::vector<double>> const& W_0);

    void apply_heat_source(double Q_tot, double x_from, double x_to, mesh const& msh);
    void apply_mass_source(double M_tot, double T, double x_from, double x_to, mesh const& msh, std::vector<double> const& comp);

    void export_to_file(mesh const& msh);
    void export_timestep(double t, mesh const& msh, std::vector<particle> const& particles);
};