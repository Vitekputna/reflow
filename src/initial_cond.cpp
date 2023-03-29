#include <vector>
#include <iostream>
#include "initial_cond.hpp"
#include "thermodynamics.hpp"

std::vector<std::vector<double>> init::step(int N, int N_var, std::vector<double> const& L, std::vector<double> const& R)
{

    N += 2; // add ghost cells

    std::vector<std::vector<double>> W_0 = std::vector<std::vector<double>>(N,std::vector<double>(N_var,0.0));

    for(int i = 0; i < N; i++)
    {
        if(i > N/2)
        {
            W_0[i] = R;
        }
        else
        {
            W_0[i] = L;
        }
    }

    return W_0;
}

std::vector<double> init::flow(int N_var, double p, double T, double u, std::vector<double> const& comp)
{
    double r = thermo::r_mix_comp(comp);

    int n_comp = comp.size();

    int n_other = N_var-n_comp-2;

    if(n_other != 0)
    {
        std::cout << "Detecting other variables than chemical components, momentum and energy.\n";
        std::cout << "Other variables will be set initialy to zero, alernatively use diff. init. func.\n";
    }

    double rho = p/r/T;

    std::vector<double> res = {rho};

    for(auto i = 1; i < n_comp; i++)
    {
        res.push_back(rho*comp[i]);
    }

    for(int i = 0; i < n_other;i++)
    {
        res.push_back(0);
    }

    res.push_back(rho*u);

    double e = (thermo::enthalpy(T,comp) + 0.5*u*u)*rho - p;

    res.push_back(e);

    return res;
}

// flow with dropplet init func...

std::vector<double> init::flow_dropplets(int N_var, double p, double T, double u, std::vector<double> const& comp,
                                         std::vector<double> drp_frac,
                                         std::vector<double> drp_count,
                                         std::vector<double> drp_mom)
{
    double r = thermo::r_mix_comp(comp);

    int n_comp = comp.size();

    int n_other = N_var-n_comp-2;

    if(N_var != n_comp + 2 + drp_frac.size() + drp_count.size() + drp_mom.size())
    {
        std::cout << "Number of parameters is not consistent with supplied parameters!\n";
        std::cout << "Exiting...\n";
        exit(0);
    }

    double rho = p/r/T; // hmotnostní koncentrace plynné fáze (hustota)

    double chi = 0; //hmotnostní koncentrace kondenzované fáze

    for(auto& frac : drp_frac)
    {
        chi += frac;
    }
    std::vector<double> res = {rho + chi};

    for(auto i = 1; i < n_comp; i++)
    {
        res.push_back(rho*comp[i]);
    }

    for(auto& count : drp_count)
    {
        res.push_back(count);
    }

    for(auto& frac : drp_frac)
    {
        res.push_back(frac);
    }

    for(auto& mom : drp_mom)
    {
        res.push_back(mom);
    }

    res.push_back((rho+chi)*u);

    double e = thermo::enthalpy(T,comp)*rho + 0.5*u*u*(rho+chi) - p;

    res.push_back(e);

    return res;
}

std::vector<std::vector<double>> init::nozzle(int N, int N_var,double md, double T0, double p0, double p2, double L_chamber, std::vector<double> const comp, mesh const& msh)
{
    std::vector<std::vector<double>> res = std::vector<std::vector<double>>(N,std::vector<double>(N_var,0.0));

    double r = thermo::r_mix_comp(comp);
    double kappa = thermo::kappa_mix_comp(comp);

    // Chamber values
    int chamber_idx;
    double rho0 = p0/r/T0;
    double u0 = md/msh.A[0]/rho0;

    for(int i = 0; i < N; i++)
    {
        if(msh.x[i] < L_chamber)
        {
            res[i][0] = rho0;    //Density
            
            for(int k = 1; k < thermo::n_comp; k++) res[i][k] = comp[k]*rho0;  //Mass concentration of species

            res[i][thermo::n_comp] = rho0*u0;  // Momentum

            res[i][thermo::n_comp+1] = rho0*thermo::enthalpy(T0,comp) + 0.5*rho0*u0*u0 - p0;
        }
        else
        {
            chamber_idx = i;
            break;
        }
    }

    //Nozzle exit values
    double u2 = sqrt(2*kappa*r*T0/(kappa-1)*(1-pow(p2/p0,(kappa-1)/kappa))); // St Venant
    double T2 = T0*pow(p0/p2,(1-kappa)/kappa);
    double rho2 = p2/r/T2;

    //Nozzle values (linear interpolation)
    double delta = L_chamber-msh.x.back(); //length of nozzle
    double L = msh.x.back();

    //Nozzle interpolated values
    double rho,u,p,T;

    double t;
    for(int i = chamber_idx; i < N; i++)
    {
        t = (L-msh.x[i])/(L-L_chamber);

        rho = t*(rho0-rho2) + rho2;
        u = t*(u0-u2) + u2;
        p = t*(p0-p2) + p2;
        T = t*(T0-T2) + T2;

        res[i][0] = rho;

        for(int k = 1; k < thermo::n_comp; k++) res[i][k] = comp[k]*rho;  //Mass concentration of species

        res[i][thermo::n_comp] = rho0*u0;

        res[i][thermo::n_comp+1] = rho*thermo::enthalpy(T,comp) + 0.5*rho*u*u - p;
    }

    return res;
}