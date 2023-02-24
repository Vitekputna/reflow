#pragma once
#include <vector>
#include "specie.hpp"

class thermo
{
public:
    static std::vector<specie> species;

    static int n_comp;
    static double thershold_comp;

    static void load_specie(specie spec);

    static double pressure(std::vector<double> const& W);
    static double speed_of_sound(std::vector<double> const& W);
    static double temperature(std::vector<double> const& W);
    static void composition(std::vector<double>& comp, std::vector<double> const& W);

    static double kappa_mix(std::vector<double> const& W);
    static double kappa_mix_comp(std::vector<double>& comp);
    static double r_mix(std::vector<double> const& W);
    static double r_mix_comp(std::vector<double>& comp);

    static double cp_mix(std::vector<double> const& W);
    static double cp_mix_comp(std::vector<double>& comp, double T);

    static double temp_new(std::vector<double> const& comp, std::vector<double> const& W);

    private:
    static double dF(std::vector<double> const& comp, double T);
};