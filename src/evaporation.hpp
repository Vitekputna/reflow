#pragma once
#include "thermodynamics.hpp"

namespace evaporation
{
    double fuel_mass_fraction(double p, double T);
    double drop_combustion(const int i, std::vector<double> const& W);
};