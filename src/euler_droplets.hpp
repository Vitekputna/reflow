#pragma once
#include "thermodynamics.hpp"

namespace euler_droplets
{
    void droplet_drag(const int i, std::vector<double> const& W, std::vector<double>& res);

    double Kelbaliyev_Ceylan(const double Re);

    double fuel_mass_fraction(double p, double T);
    double drop_combustion_steady(const int i, std::vector<double> const& W, std::vector<double>& res);
    double drop_combustion_convective(const int i, std::vector<double> const& W, std::vector<double>& res);
};