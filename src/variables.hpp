#pragma once
#include <vector>
#include <string>

#include "mesh.hpp"

struct variables
{
    std::vector<std::vector<double>> W;
    std::vector<std::vector<double>> flux;
    std::vector<std::vector<double>> exact_flux;

    std::vector<double> q;
    std::vector<double> md;

    int N_var; // number of variables
    int N; // number of cells
    int N_walls; // number of walls
    int N_ghosts = 2; // number of ghost cells

    variables();
    variables(int _N_var, int _N);
    variables(int _N_var, int _N, std::vector<double> const& W_0);
    variables(int _N_var, int _N, std::vector<std::vector<double>> const& W_0);

    void apply_heat_source(double Q_tot, double x_from, double x_to, mesh const& msh);

    void export_to_file(std::string path, mesh const& msh);
};