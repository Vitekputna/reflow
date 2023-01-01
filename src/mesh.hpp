#pragma once
#include <vector>

struct mesh
{
    std::vector<double> x;
    std::vector<double> xf;
    std::vector<double> A;
    std::vector<double> Af;

    double dx_min;

    int N; // number of cells

    // x bounds
    double x_from = 0;
    double x_to = 1;

    // A bounds
    double A_from = 1;
    double A_to = 1;

    mesh();
    mesh(int _N);
    mesh(int _N, double from, double to);

    void bump();
    void cubic(std::vector<std::vector<std::vector<double>>> curves, int n);

    void export_to_file();
};