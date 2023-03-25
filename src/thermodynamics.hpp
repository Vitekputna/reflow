#pragma once
#include <vector>
#include "specie.hpp"

class thermo
{
public:
    static std::vector<specie> species;

    static std::vector<double> T;
    static std::vector<double> p;

    static int n_comp;
    static double thershold_comp;

    thermo();

    static void init(int n);

    static void load_specie(specie spec);

    static void update(std::vector<std::vector<double>> const& W);

    static double density(std::vector<double> const& W);
    static double pressure(std::vector<double> const& W);
    static double speed_of_sound(int i, std::vector<double> const& W);
    static double temperature(std::vector<double> const& W);
    static double enthalpy(double T, std::vector<double> const& comp);
    static double enthalpy(int i, std::vector<double> const& W);
    static double enthalpy_stagnate(int i, std::vector<double> const& W);
    static void composition(std::vector<double>& comp, std::vector<double> const& W);

    static double kappa_mix(std::vector<double> const& W);
    static double kappa_mix_comp(std::vector<double> const& comp);
    static double r_mix(std::vector<double> const& W);
    static double r_mix_comp(std::vector<double> const& comp);

    static double cp_mix(std::vector<double> const& W);
    static double cp_mix_comp(std::vector<double> const& comp, double T);

    static double temp_new(int i, std::vector<double> const& comp, std::vector<double> const& W);

    private:
    static double dF(std::vector<double> const& comp, double r, double T);

    static inline double pressure(int i, std::vector<double> const& W, std::vector<double> const& comp);
};