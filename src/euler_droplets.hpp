#pragma once
#include "thermodynamics.hpp"

namespace euler_droplets
{
    void droplet_drag(const int i, std::vector<double> const& W, std::vector<double>& res);
    void droplet_heat(const int i, std::vector<double> const& W, std::vector<double>& res);

    double Kelbaliyev_Ceylan(const double Re);
    double Ingebo(const double Re);

    double Ranz_Marshall(const double Re, const double Pr);
    double Sherwood_evaporation(const double Re, const double Sc, const double BM);
    double Nusselt_evaporation(const double Re,const double Pr, const double BT);

    double droplet_temperature();

    double heat_evap_interp(double T, double T_boil);

    double fuel_mass_fraction(double p, double T);
    double droplet_evaporation(const int i, std::vector<double>& W, std::vector<double>& res);
    double drop_combustion_steady(const int i, std::vector<double> const& W, std::vector<double>& res);
    double drop_combustion_convective(const int i, std::vector<double> const& W, std::vector<double>& res);
};