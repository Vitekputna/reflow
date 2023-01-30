#include "thermodynamics.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include "specie.hpp"

// Using ideal gas law
double thermodynamics::pressure(std::vector<double> const& W, double kappa)
{
    int mom_idx = W.size()-2;
    return (kappa-1)*(W[mom_idx+1] - 0.5*W[mom_idx]*W[mom_idx]/W[0]);
}

double thermodynamics::speed_of_sound(std::vector<double> const& W, double kappa)
{
    return sqrt(kappa*pressure(W,kappa)/W[0]);
}

double thermodynamics::temperature(std::vector<double> const& W, double kappa, double r)
{
    return thermodynamics::pressure(W,kappa)/r/W[0];
}

std::vector<double> thermodynamics::composition(std::vector<double> const& W)
{
    int n_comp = W.size()-2;
    std::vector<double> comp(n_comp,0.0);

    for(auto idx = 0; idx < n_comp; idx++)
    {
        comp[idx] = W[idx]/W[0];
    }

    return comp;
}

// // Using ideal gas law
// double thermo::pressure(std::vector<double> const& W, double kappa)
// {
//     int mom_idx = W.size()-2;
//     return (kappa-1)*(W[mom_idx+1] - 0.5*W[mom_idx]*W[mom_idx]/W[0]);
// }

// double thermo::speed_of_sound(std::vector<double> const& W, double kappa)
// {
//     return sqrt(kappa*pressure(W,kappa)/W[0]);
// }

// double thermo::temperature(std::vector<double> const& W, double kappa, double r)
// {
//     return thermo::pressure(W,kappa)/r/W[0];
// }

// std::vector<double> thermo::composition(std::vector<double> const& W)
// {
//     int n_comp = W.size()-2;
//     std::vector<double> comp(n_comp,0.0);

//     for(auto idx = 0; idx < n_comp; idx++)
//     {
//         comp[idx] = W[idx]/W[0];
//     }

//     return comp;
// }
