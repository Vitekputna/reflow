#include "thermodynamics.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include "specie.hpp"

std::vector<specie> thermo::species;
int thermo::n_comp = 0;

// Using ideal gas law
double thermo::pressure(std::vector<double> const& W)
{
    int mom_idx = W.size()-2;
    double kappa = thermo::kappa_mix(W);
    return (kappa-1)*(W[mom_idx+1] - 0.5*W[mom_idx]*W[mom_idx]/W[0]);
}

double thermo::speed_of_sound(std::vector<double> const& W)
{
    double kappa = thermo::kappa_mix(W);
    return sqrt(kappa*pressure(W)/W[0]);
}

double thermo::temperature(std::vector<double> const& W)
{
    double r = thermo::r_mix(W);
    return thermo::pressure(W)/r/W[0];
}

void thermo::composition(std::vector<double>& comp, std::vector<double> const& W)
{
    double sum = 0;

    for(auto idx = 1; idx < n_comp; idx++)
    {
        comp[idx] = W[idx]/W[0];
        sum += comp[idx];
    }
    comp[0] = 1-sum;
}

void thermo::load_specie(specie spec)
{
    thermo::species.push_back(spec);
    thermo::n_comp++;
}

double thermo::kappa_mix(std::vector<double> const& W)
{
    double kappa = 0;
    std::vector<double> comp(n_comp); 
    thermo::composition(comp,W);

    for(auto idx = 0; idx < n_comp; idx++)
    {   
        kappa += thermo::species[idx].kappa*comp[idx];
    }

    return kappa;
}

double thermo::kappa_mix_comp(std::vector<double>& comp)
{
    double kappa = 0;

    for(auto idx = 0; idx < n_comp; idx++)
    {
        kappa += thermo::species[idx].kappa*comp[idx];
    }

    return kappa;
}

double thermo::r_mix(std::vector<double> const& W)
{
    double r = 0;
    std::vector<double> comp(n_comp);
    thermo::composition(comp,W);

    for(auto idx = 0; idx < n_comp; idx++)
    {
        r += thermo::species[idx].r*comp[idx];
    }

    // std::cout << r << "\n";

    return r;
}

double thermo::r_mix_comp(std::vector<double>& comp)
{
    double r = 0;
    for(auto idx = 0; idx < n_comp; idx++)
    {
        r += thermo::species[idx].r*comp[idx];
    }
    
    return r;
}