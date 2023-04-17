#pragma once
#include <vector>

struct mesh
{
    std::vector<double> xf;
    std::vector<double> Af;

    std::vector<double> x;
    std::vector<double> A;
    std::vector<double> V;

    double dx_min = 1e15;
    
    int N; // number of cells

    // x bounds
    double x_from = 0;
    double x_to = 1;

    // A bounds
    double A_from = 1;
    double A_to = 1;

    // 
    int N_smooth_cycles = 300;

    mesh();
    mesh(int _N);
    mesh(int _N, double from, double to);
    void refine(std::vector<std::vector<double>> ref);

    void construct_mesh();
    void compute_volumes();

    void fix_cell_centroids();
    void smooth_mesh();

    void bump();
    void cubic(std::vector<std::vector<std::vector<double>>> curves, int n);
    void constant_area(double Area);

    void export_to_file();
};