#pragma once
#include <vector>

namespace thermo
{
    double pressure(std::vector<double> const& W, double kappa);
    double speed_of_sound(std::vector<double> const& W, double kappa);
    double temperature(std::vector<double> const& W, double kappa, double r);
}

