#pragma once
#include <vector>
#include <string>

struct variables
{
    std::vector<std::vector<double>> W;
    std::vector<std::vector<double>> flux;
    std::vector<std::vector<double>> exact_flux;

    int N_var; // number of variables
    int N; // number of cells
    int N_walls; // number of walls
    int N_ghosts = 2; // number of ghost cells

    variables();
    variables(int _N_var, int _N);
    variables(int _N_var, int _N, std::vector<double> W_0);

    void export_to_file(std::string path);
};