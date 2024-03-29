#include "evaporation.hpp"
#include "variables.hpp"
#include <iostream>

double evaporation::fuel_mass_fraction(double p, double T)
{
    const double p_ref = thermo::species[2].p_ref;
    const double T_ref = thermo::species[2].T_ref;
    const double h_vap = thermo::species[2].h_vap;
    const double Mm_f = thermo::species[2].Mm;
    const double Mm_p = thermo::species[0].Mm;

    double K = p_ref*exp(h_vap*Mm_f/8314*(1/T_ref - 1/T))/p;

    return std::min(0.99,Mm_f*K/(Mm_p - (Mm_p-Mm_f)*K));
}

double evaporation::drop_combustion_steady(const int i, std::vector<double> const& W, std::vector<double>& res)
{
    const double T_liq = 0.9*thermo::species[2].T_ref;
    const double T_f = thermo::T[i];
    const double T_ref = T_liq + (T_f - T_liq)/3;

    const double h_vap = thermo::species[2].h_vap;
    const double rho_l = thermo::species[2].rho_liq;

    // const double Yf = fuel_mass_fraction(thermo::p[i],T_liq);
    const double Yf = 1-1e-8;
    const double Yf_ref = 2*Yf/3;
    const double Yp_ref = 1-Yf_ref;

    const std::vector<double> comp = {Yp_ref,0,Yf_ref};

    const double BM = Yf/(1-Yf);

    const double k = thermo::thermal_conductivity(comp,T_ref);
    const double cp = thermo::cp_mix_comp(comp,T_ref);

    const double BT = cp*(T_f - T_liq)/h_vap;

    const int N_comp = variables::N_comp;

    int N_idx, Frac_idx;
    double r, dm, D, dm_total = 0;

    for(int j = 0; j < variables::N_drop_frac; j++)
    {
        N_idx = N_comp + 2*j;
        Frac_idx = N_comp + 2*j + 1;

        r = std::pow(3*W[Frac_idx]/(4*W[N_idx]*M_PI*rho_l),0.3333);

        if(W[N_idx] == 0 || r < 0) r = 0;

        dm = W[N_idx]*4*4.1*M_PI*(k/cp)*r*log(1+BM);

        dm = std::max(0.0,dm);

        dm_total += dm;

        res[Frac_idx] -= dm;
        res[2] += dm;
        res[variables::eng_idx] += dm*thermo::enthalpy(thermo::T[i],std::vector<double>{0,0,1});
    }

    return dm_total;
}

double evaporation::drop_combustion_convective(const int i, std::vector<double> const& W, std::vector<double>& res)
{
    const double T_liq = thermo::species[2].T_ref;
    const double T_f = thermo::T[i];
    const double T_ref = T_liq + (T_f - T_liq)/3;
    const double h_vap = thermo::species[2].h_vap;
    const double rho_l = thermo::species[2].rho_liq;

    const std::vector<double> comp = std::vector<double>{0,0,1};

    const double lambda = thermo::thermal_conductivity(comp,T_ref);
    const double cp = thermo::cp_mix_comp(comp,T_ref);
    const double mu = thermo::viscosity(comp,T_ref);
    const double nu = mu/W[0];

    const int N_comp = variables::N_comp;

    int N_idx, Frac_idx;
    double r, dm, dm_total = 0;
    double Nu, Re, Pr;

    for(int j = 0; j < variables::N_drop_frac; j += 2)
    {
        N_idx = N_comp + j;
        Frac_idx = N_comp + j + 1;

        r = std::pow(3*W[Frac_idx]/(4*W[N_idx]*M_PI*rho_l),0.3333);

        if(W[N_idx] == 0) r = 0;

        // Re = 2*r/nu;
        // Pr = nu*cp*W[0]/(lambda);
        // Nu = 2 + 0.6*pow(Re,0.5)*pow(Pr,0.33);  

        // dm = W[N_idx]*2*Nu*M_PI*r*(lambda/cp)*log(1+(cp*(T_f-T_liq))/(h_vap));
        dm = W[N_idx]*4*M_PI*r*(lambda/cp)*log(1+(cp*(T_f-T_liq))/(h_vap));
        dm = std::max(0.0,dm);

        dm_total += dm;

        res[Frac_idx] -= dm;
        res[2] += dm;
        res[variables::eng_idx] += dm*thermo::enthalpy(thermo::T[i],std::vector<double>{0,0,1});
    }

    return dm_total;

}