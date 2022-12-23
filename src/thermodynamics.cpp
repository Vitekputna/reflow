#include "thermodynamics.hpp"
#include <vector>

double thermo::pressure(std::vector<double> const& W, double gamma)
{
    return (gamma-1)*(W[2] - 0.5*W[1]*W[1]/W[0]);
}