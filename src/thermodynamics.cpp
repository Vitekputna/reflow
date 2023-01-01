#include "thermodynamics.hpp"
#include <vector>
#include <cmath>

double thermo::pressure(std::vector<double> const& W, double kappa)
{
    return (kappa-1)*(W[2] - 0.5*W[1]*W[1]/W[0]);
}

double thermo::speed_of_sound(std::vector<double> const& W, double kappa)
{
    return sqrt(kappa*pressure(W,kappa)/W[0]);
}

double thermo::temperature(std::vector<double> const& W, double kappa, double r)
{
    return thermo::pressure(W,kappa)/r/W[0];
}