#pragma once
#include <vector>

namespace init
{
    std::vector<std::vector<double>> step(int N, int N_var, std::vector<double> const& L, std::vector<double> const& R);
    std::vector<double> flow(int N_var, double p, double T, double u, std::vector<double> const& comp);
}