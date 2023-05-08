#include "euler_droplets.hpp"
#include "variables.hpp"
#include <iostream>

double euler_droplets::fuel_mass_fraction(double p, double T)
{
    const double p_ref = thermo::species[2].p_ref;
    const double T_ref = thermo::species[2].T_ref;
    const double h_vap = thermo::species[2].h_vap;
    const double Mm_f = thermo::species[2].Mm;
    const double Mm_p = thermo::species[0].Mm;

    double K = p_ref*exp(h_vap*Mm_f/8314*(1/T_ref - 1/T))/p;

    return std::min(0.99,Mm_f*K/(Mm_p - (Mm_p-Mm_f)*K));
}

double euler_droplets::drop_combustion_steady(const int i, std::vector<double> const& W, std::vector<double>& res)
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

double euler_droplets::drop_combustion_convective(const int i, std::vector<double> const& W, std::vector<double>& res)
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

// Handbook of atomization and sprays
double euler_droplets::Kelbaliyev_Ceylan(const double Re)
{
    return (24/Re)*pow(1 + 18.5*pow(Re,3.6) + pow(Re/2,11),1/30) + (4/9)*pow(Re,4.5)/(330+pow(Re,4/5));
}


double euler_droplets::Ranz_Marshall(const double Re, const double Pr)
{
    return 2+0.6*pow(Re,0.5)*pow(Pr,1/3);
}

void euler_droplets::droplet_drag(const int i, std::vector<double> const& W, std::vector<double>& res)
{
    int mom_idx, frac_idx, num_idx;
    double A,r,u_drop,u_gas,Cd,Re;

    const double rho_l = thermo::species[2].rho_liq;
    const double rho_gas = thermo::density(W);

    std::vector<double> comp = std::vector<double>(variables::N_comp,0);
    thermo::composition(comp,W);

    const double mu = thermo::viscosity(comp,thermo::T[i]);

    for(int idx = 0; idx < variables::active_drop_idx.size(); idx++)
    {
        num_idx = variables::active_drop_idx[idx]-1;

        if(W[num_idx] == 0) continue;

        frac_idx = variables::active_drop_idx[idx];
        mom_idx = variables::drop_mom_idx[idx];

        r = std::pow(3*W[frac_idx]/(4*W[num_idx]*M_PI*rho_l),0.3333);

        A = M_PI*pow(r,2);
        u_drop = W[mom_idx]/(W[frac_idx]+1e-12);
        u_gas = W[variables::mom_idx]/W[0];

        Re = (rho_gas*abs(u_gas-u_drop)*r)/mu + 1e-12;

        // Cd = Kelbaliyev_Ceylan(Re);
        Cd = 0.45;

        res[mom_idx] += 0.5*Cd*W[num_idx]*thermo::density(W)*A*abs(u_gas-u_drop)*(u_gas-u_drop);
    }
}

void euler_droplets::droplet_heat(const int i, std::vector<double> const& W, std::vector<double>& res)
{
    int mom_idx, frac_idx, num_idx, eng_idx;
    double r,u_drop,T_drop,Nu,Re,Pr;

    const double rho_l = thermo::species[2].rho_liq;
    const double rho_gas = thermo::density(W);

    const double T_gas = thermo::T[i];
    const double u_gas = W[variables::mom_idx]/W[0];

    std::vector<double> comp = std::vector<double>(variables::N_comp,0);
    thermo::composition(comp,W);

    const double mu = thermo::viscosity(comp,thermo::T[i]);
    const double cp = thermo::cp_mix_comp(comp,thermo::T[i]);
    const double k = thermo::thermal_conductivity(comp,thermo::T[i]);

    for(int idx = 0; idx < variables::active_drop_idx.size(); idx++)
    {
        num_idx = variables::active_drop_idx[idx]-1;

        if(W[num_idx] == 0) continue;

        frac_idx = variables::active_drop_idx[idx];
        mom_idx = variables::drop_mom_idx[idx];
        eng_idx = variables::N_comp+2*variables::N_drop_frac+variables::N_drop_eng_eq+idx;

        r = std::pow(3*W[frac_idx]/(4*W[num_idx]*M_PI*rho_l),0.3333);

        T_drop = W[eng_idx]/(W[frac_idx]*thermo::species[2].C + 1e-12);
        u_drop = W[mom_idx]/(W[frac_idx]+1e-12);

        Re = (rho_gas*abs(u_gas-u_drop)*r)/mu + 1e-12;
        Pr = cp*mu/k;

        Nu = Ranz_Marshall(Re,Pr);

        res[eng_idx] += W[num_idx]*Nu*2*r*M_PI*k*(T_gas-T_drop);
    }
}