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

// std::vector<std::vector<double>> init::nozzle(int N, int N_var,double md, double T0, double p0, std::vector<double> const comp, mesh const& msh)
// {
//     return std::vector<std::vector<double>>{};
// }