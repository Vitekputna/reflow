#include "euler_droplets.hpp"
#include "variables.hpp"
#include <iostream>
#include <cmath>

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

double euler_droplets::heat_evap_interp(double T, double T_boil)
{
    // Cubic
    const double b = 0.5*T_boil;
    const double a = T_boil;

    const double A = a*a*a/3-a*a/2*(a+b)+a*a*b;
    const double B = b*b*b/3-b*b/2*(a+b)+a*b*b;

    const double c = -1/(A-B);
    const double d = -A/(A-B)+1;

    if(T < b)
    {
        return 0.0;
    }
    else if(T >= b && T <= a)
    {
        return -c*((T*T*T)/3-(T*T/2)*(a+b)+a*b*T)+d;
    }
    else
    {
        return 1.0;
    }
}

// Handbook of atomization and sprays
double euler_droplets::Kelbaliyev_Ceylan(const double Re)
{
    return (24/Re)*pow(1 + 18.5*pow(Re,3.6) + pow(Re/2,11),1/30) + (4/9)*pow(Re,4.5)/(330+pow(Re,4/5));
}

// Multiphase flows with droplets and particles str. 73
double euler_droplets::Ingebo(const double Re)
{
    return 27*pow(Re,-0.84);
}

// Multiphase flows with droplets and particles
double euler_droplets::Ranz_Marshall(const double Re, const double Pr)
{
    return 2+0.6*pow(Re,0.5)*pow(Pr,1/3);
}

double euler_droplets::Sherwood_evaporation(const double Re, const double Sc, const double BM)
{
    return (2+0.87*pow(Re,0.5)*pow(Sc,0.33333))/pow(1+BM,0.7);
}

double euler_droplets::Nusselt_evaporation(const double Re,const double Pr, const double BT)
{
    return (2+0.57*pow(Re,0.5)*pow(Pr,0.33333))/pow(1+BT,0.7);
}

double euler_droplets::droplet_evaporation(const int i, std::vector<double>& W, std::vector<double>& res)
{
    int mom_idx, frac_idx, num_idx, eng_idx;
    double r,u_drop,T_drop,Ys,Re,Sh,Sc,Pr,Nu;
    double BT,BM;
    double total_m = 0;
    double md,Q;
    double phi;

    const double rho_l = thermo::species[2].rho_liq;
    const double rho_gas = thermo::density(W);

    const double T_gas = thermo::T[i];
    const double u_gas = W[variables::mom_idx]/W[0];

    std::vector<double> comp = std::vector<double>(variables::N_comp,0);
    thermo::composition(comp,W);

    const double D = thermo::difusivity(comp,thermo::T[i]);
    const double mu = thermo::viscosity(comp,thermo::T[i]);
    const double cp = thermo::cp_mix_comp(comp,thermo::T[i]);
    const double k = thermo::thermal_conductivity(comp,thermo::T[i]);
    const double h_vap = thermo::species[2].h_vap;
    const double Y_gas = comp[2];

    const double T_boil = thermo::boil_temp(2,thermo::p[i]);

    for(int idx = 0; idx < variables::active_drop_idx.size(); idx++)
    {
        num_idx = variables::active_drop_idx[idx]-1;
        frac_idx = variables::active_drop_idx[idx];
        mom_idx = variables::drop_mom_idx[idx];
        eng_idx = variables::N_comp+2*variables::N_drop_frac+variables::N_drop_eng_eq+idx;

        r = std::pow(3*W[frac_idx]/(4*W[num_idx]*M_PI*rho_l),0.3333);

        if(W[num_idx] <= 0 || r < 1e-6)
        {
            continue;
        }

        u_drop = W[mom_idx]/(W[frac_idx]+1e-12);
        T_drop = W[eng_idx]/(W[frac_idx]*thermo::species[2].C + 1e-12);

        Ys = fuel_mass_fraction(thermo::p[i],T_drop);

        Re = (rho_gas*std::abs(u_gas - u_drop)*2*r)/mu;
        Sc = mu/(rho_gas*D);
        Pr = cp*mu/k;

        phi = heat_evap_interp(T_drop,T_boil);

        BT = cp*(T_gas-T_drop)*phi/h_vap;
        BM = (Ys-Y_gas)/(1-Ys);

        Nu = Nusselt_evaporation(Re,Pr,BT);
        Sh = Sherwood_evaporation(Re,Sc,BM);

        // Nu = Ranz_Marshall(Re,Pr);
        // Sh = 2+0.6*pow(Re,0.5)*pow(Sc,0.333);

        md = W[num_idx]*Sh*2*r*M_PI*D*rho_gas*(Ys-Y_gas);
        Q = W[num_idx]*Nu*2*r*M_PI*k*(T_gas-T_drop);
        
        

        md += std::max(0.0,phi*Q/h_vap);
        Q += -std::max(0.0,phi*Q);  
            
        total_m += md;

        res[frac_idx] += -md;
        res[mom_idx] += -md*u_drop;
        res[eng_idx] += Q - md*(h_vap + thermo::species[2].C*T_drop);

        res[0] += md;
        res[2] += md;
        res[variables::mom_idx] += md*u_drop;
        res[variables::eng_idx] += -Q + md*thermo::species[2].h(T_drop) + 0.5*md*u_drop*u_drop;
        // res[variables::eng_idx] += md*thermo::species[2].h(T_drop) + 0.5*md*u_drop*u_drop;
    }

    return total_m;
}

void euler_droplets::droplet_drag(const int i, std::vector<double> const& W, std::vector<double>& res)
{
    int mom_idx, frac_idx, num_idx;
    double A,r,u_drop,Cd,Re,Fd;

    const double rho_l = thermo::species[2].rho_liq;
    const double rho_gas = thermo::density(W);

    const double u_gas = W[variables::mom_idx]/W[0];

    std::vector<double> comp = std::vector<double>(variables::N_comp,0);
    thermo::composition(comp,W);

    const double mu = thermo::viscosity(comp,thermo::T[i]);

    const double gas_vol_frac = 1-thermo::liquid_fraction(W);

    for(int idx = 0; idx < variables::active_drop_idx.size(); idx++)
    {
        num_idx = variables::active_drop_idx[idx]-1;
        frac_idx = variables::active_drop_idx[idx];
        mom_idx = variables::drop_mom_idx[idx];

        r = std::pow(3*W[frac_idx]/(4*W[num_idx]*M_PI*rho_l),0.3333);

        if(W[num_idx] <= 0 || r < 1e-6)
        {
            continue;
        }

        A = M_PI*pow(r,2);
        u_drop = W[mom_idx]/(W[frac_idx]+1e-12);
        
        Re = (rho_gas*std::abs(u_gas-u_drop)*2*r)/mu + 1e-12;

        Cd = Ingebo(Re);
        // Cd = 0.45;
        // Cd = 24/Re;

        Fd = 0.5*Cd*W[num_idx]*rho_gas*A*std::abs(u_gas-u_drop)*(u_gas-u_drop);

        res[mom_idx] += Fd;
        // res[variables::mom_idx] += Fd;
        // res[variables::eng_idx] -= Fd*u_drop;
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

        frac_idx = variables::active_drop_idx[idx];
        mom_idx = variables::drop_mom_idx[idx];
        eng_idx = variables::N_comp+2*variables::N_drop_frac+variables::N_drop_eng_eq+idx;

        r = std::pow(3*W[frac_idx]/(4*W[num_idx]*M_PI*rho_l),0.3333);

        if(W[num_idx] <= 0 || r < 1e-7)
        {
            continue;
        }

        T_drop = W[eng_idx]/(W[frac_idx]*thermo::species[2].C + 1e-12);
        u_drop = W[mom_idx]/(W[frac_idx]+1e-12);

        Re = (rho_gas*std::abs(u_gas-u_drop)*2*r)/mu + 1e-12;
        Pr = cp*mu/k;

        Nu = Ranz_Marshall(Re,Pr);

        res[eng_idx] += W[num_idx]*Nu*2*r*M_PI*k*(T_gas-T_drop);
        // res[variables::eng_idx] += -W[num_idx]*Nu*2*r*M_PI*k*(T_gas-T_drop);
    }
}