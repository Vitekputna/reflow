#pragma once
#include <vector>
#include "specie.hpp"

class thermo
{
public:
    static std::vector<specie> species;

    static std::vector<double> T;
    static std::vector<double> p;

    static std::vector<double> Pr;
    static std::vector<double> Sc;

    static int n_comp;
    static double thershold_comp;

    thermo();

    static void init(int n);

    static void load_specie(specie spec);

    static void update(std::vector<std::vector<double>> const& W);
    static void update(std::vector<std::vector<double>> const& W, const int from, const int to);

    static double density(std::vector<double> const& W);
    static double pressure(std::vector<double> const& W);
    static double speed_of_sound(int i, std::vector<double> const& W);
    static double speed_of_sound(const double p, std::vector<double> const& W);
    static double mach_number(int i, std::vector<double> const& W);
    static double mach_number(const double p, std::vector<double> const& W);
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

    static double difusivity(std::vector<double> const& comp, double T);

    static std::vector<double> molar_fraction(std::vector<double> const& mass_fraction);
    static std::vector<double> mass_fraction(std::vector<double> const& molar_fraction);
    static double liquid_fraction(std::vector<double> const& W);

    static double viscosity(std::vector<double> const& comp, double T);

    static double thermal_conductivity(std::vector<double> const& comp, double T);

    private:
    static double temp_new(int i, std::vector<double> const& comp, std::vector<double> const& W);
    static double dF(std::vector<double> const& comp, double r, double T);
    static inline double pressure(int i, std::vector<double> const& W, std::vector<double> const& comp);
};