#pragma once
#include <vector>

struct mesh
{
    std::vector<double> x;
    std::vector<double> xf;
    std::vector<double> A;

    int N; // number of cells

    // x bounds
    double x_from = 0;
    double x_to = 1;

    // A bounds
    double A_from = 1;
    double A_to = 1;

    mesh();
    mesh(int _N);


};