#pragma once
#include <vector>
#include "specie.hpp"

class thermodynamics
{
    std::vector<specie> species;

    double pressure(std::vector<double> const& W, double kappa);
    double speed_of_sound(std::vector<double> const& W, double kappa);
    double temperature(std::vector<double> const& W, double kappa, double r);
    std::vector<double> composition(std::vector<double> const& W);
};

// namespace thermo
// {
//     double pressure(std::vector<double> const& W, double kappa);
//     double speed_of_sound(std::vector<double> const& W, double kappa);
//     double temperature(std::vector<double> const& W, double kappa, double r);
//     std::vector<double> composition(std::vector<double> const& W);
// }