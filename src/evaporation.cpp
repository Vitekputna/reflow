#include "evaporation.hpp"
#include "variables.hpp"

double evaporation::fuel_mass_fraction(double p, double T)
{
    const double p_ref = thermo::species[2].p_ref;
    const double T_ref = thermo::species[2].T_ref;
    const double h_vap = thermo::species[2].h_vap;
    const double Mm_f = thermo::species[2].Mm;
    const double Mm_p = thermo::species[0].Mm;

    double K = p_ref*exp(h_vap*Mm_f/8314*(1/T_ref - 1/T))/p;

    return std::min(1.0,Mm_f*K/(Mm_p - (Mm_p-Mm_f)*K));
}

double evaporation::drop_combustion(const int i, std::vector<double> const& W)
{
    const double Yf = fuel_mass_fraction(thermo::p[i],340);
    const double Yo = W[1]/W[0];
    const double rho_l = thermo::species[2].rho_liq;

    const double B = (1+0.1515*Yo)/(1-Yf);

    const int N_comp = variables::N_comp;

    int N_idx, Frac_idx;
    double r, dm = 0, D;

    D = 1e-5;

    for(int j = 0; j < variables::N_drop_frac; j += 2)
    {
        if(W[N_idx] == 0) continue;

        N_idx = N_comp + j;
        Frac_idx = N_comp + j + 1;

        r = std::pow(3*W[Frac_idx]/(4*W[N_idx]*3.14159*rho_l),0.3333);

        dm += W[N_idx]*4*M_PI*W[0]*D*r*log(1+B);
    }

    return dm;
}