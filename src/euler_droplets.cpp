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

    return std::min(1.0,Mm_f*K/(Mm_p - (Mm_p-Mm_f)*K));
}

// Handbook of atomization and sprays
double euler_droplets::Kelbaliyev_Ceylan(const double Re)
{
    return (24/Re)*pow(1 + 18.5*pow(Re,3.6) + pow(Re/2,11),1/30) + (4/9)*pow(Re,4.5)/(330+pow(Re,4/5));
}

// Multiphase flows with droplets and particles
double euler_droplets::Ranz_Marshall(const double Re, const double Pr)
{
    return 2+0.6*pow(Re,0.5)*pow(Pr,1/3);
}

double euler_droplets::droplet_evaporation(const int i, std::vector<double>& W, std::vector<double>& res)
{
    int mom_idx, frac_idx, num_idx, eng_idx;
    double r,u_drop,u_gas,T_drop,Ys,Re,Sh,Sc;
    double total_m = 0;
    double md;

    const double rho_l = thermo::species[2].rho_liq;
    const double rho_gas = thermo::density(W);

    std::vector<double> comp = std::vector<double>(variables::N_comp,0);
    thermo::composition(comp,W);

    const double D = thermo::difusivity(comp,thermo::T[i]);

    const double mu = thermo::viscosity(comp,thermo::T[i]);

    const double Y_gas = comp[2];

    for(int idx = 0; idx < variables::active_drop_idx.size(); idx++)
    {
        num_idx = variables::active_drop_idx[idx]-1;
        frac_idx = variables::active_drop_idx[idx];
        mom_idx = variables::drop_mom_idx[idx];

        if(W[num_idx] <= 0 || W[frac_idx] < 1e-12)
        {
            continue;
        }

        
        eng_idx = variables::N_comp+2*variables::N_drop_frac+variables::N_drop_eng_eq+idx;

        r = std::pow(3*W[frac_idx]/(4*W[num_idx]*M_PI*rho_l),0.3333);

        // u_drop = W[mom_idx]/(W[frac_idx]+1e-12);
        u_drop = 50;
        u_gas = W[variables::mom_idx]/W[0];
        T_drop = W[eng_idx]/(W[frac_idx]*thermo::species[2].C + 1e-12);

        Ys = fuel_mass_fraction(thermo::p[i],T_drop);

        Re = (rho_gas*abs(u_gas-u_drop)*r)/mu + 1e-12;
        Sc = mu/(rho_gas*D);

        Sh = 2+0.6*pow(Re,0.5)*pow(Sc,0.333);

        md = W[num_idx]*2*M_PI*r*D*Sh*rho_gas*(Ys-Y_gas);

        total_m += md;

        res[frac_idx] += -md;
        res[mom_idx] += -md*u_drop;
        // res[eng_idx] += -md*thermo::species[2].C*T_drop;
    }
    return total_m;
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

        Cd = Kelbaliyev_Ceylan(Re);
        // Cd = 0.45;

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