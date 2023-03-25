#include "thermodynamics.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include "specie.hpp"
#include "variables.hpp"

std::vector<specie> thermo::species;

std::vector<double> thermo::T, thermo::p;

int thermo::n_comp = 0;
double thermo::thershold_comp = 1e-4;

thermo::thermo() {}

void thermo::init(int n)
{
    thermo::p = std::vector<double>(n,1e5);
    thermo::T = std::vector<double>(n,300);
}

// Using ideal gas law

double thermo::density(std::vector<double> const& W)
{
    static double chi;
    chi = 0;

    for(int idx = 1; idx < variables::N_drop_frac; idx+=2)
    {
        chi += W[idx + n_comp];
    }

    return W[0]-chi;
}

double thermo::speed_of_sound(int i, std::vector<double> const& W)
{
    double kappa = thermo::kappa_mix(W);
    // return sqrt(kappa*thermo::p[i]/W[0]);
    return sqrt(kappa*thermo::p[i]/density(W));
}

double thermo::enthalpy(int i, std::vector<double> const& W)
{
    static std::vector<double> comp(n_comp);
    thermo::composition(comp,W);

    double T = thermo::T[i];

    double h = 0;
    for(int j = 0; j < thermo::n_comp; j++){h += comp[j]*thermo::species[j].h(T);}
    return h;
}

double thermo::enthalpy(double T, std::vector<double> const& comp)
{
    double h = 0;
    for(int j = 0; j < thermo::n_comp; j++){h += comp[j]*thermo::species[j].h(T);}
    return h;
}

double thermo::enthalpy_stagnate(int i, std::vector<double> const& W)
{
    int n = W.size();

    return thermo::enthalpy(i,W) + pow(W[n-2]/W[0],2)/2;
}

double thermo::temperature(std::vector<double> const& W)
{
    static std::vector<double> comp(n_comp);
    thermo::composition(comp,W);

    return temp_new(0,comp,W);
}

void thermo::composition(std::vector<double>& comp, std::vector<double> const& W)
{
    static double sum;
    static double rho;
    
    rho = density(W);
    sum = 0;
    for(auto idx = 1; idx < n_comp; idx++)
    {
        comp[idx] = W[idx]/rho;
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
    static std::vector<double> comp(n_comp); 
    thermo::composition(comp,W);

    for(auto idx = 0; idx < n_comp; idx++)
    {   
        kappa += thermo::species[idx].kappa*comp[idx];
    }

    return kappa;
}

double thermo::kappa_mix_comp(std::vector<double> const& comp)
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
    static std::vector<double> comp(n_comp);
    thermo::composition(comp,W);

    for(auto idx = 0; idx < n_comp; idx++)
    {
        r += thermo::species[idx].r*comp[idx];
    }

    return r;
}

double thermo::r_mix_comp(std::vector<double> const& comp)
{
    double r = 0;
    for(auto idx = 0; idx < n_comp; idx++)
    {
        r += thermo::species[idx].r*comp[idx];
    }
    
    return r;
}

double thermo::cp_mix(std::vector<double> const& W)
{
    static std::vector<double> comp(n_comp);
    thermo::composition(comp,W);
    double T = thermo::temp_new(0,comp,W);
    double cp = 0;

    for(int i = 0; i < n_comp; i++)
    {
        cp += comp[i]*thermo::species[i].cp(T);
    }

    return cp;
}

double thermo::cp_mix_comp(std::vector<double> const& comp, double T)
{
    double cp = 0;

    for(int i = 0; i < n_comp; i++)
    {
        cp += comp[i]*thermo::species[i].cp(T);
    }

    return cp;
}

double thermo::dF(std::vector<double> const& comp, double r, double T)
{
    double df = 0;

    for(int i = 0; i < n_comp; i++)
    {
        if(comp[i]>thershold_comp)
        {
            df += comp[i]*thermo::species[i].cp(T); 
        }
    }

    return df-r;
}

double thermo::temp_new(int idx, std::vector<double> const& comp, std::vector<double> const& W)
{
    int n = W.size()-1;
    static double rho;
    rho = density(W);
    // double C = W[n]/W[0] - 0.5*W[n-1]*W[n-1]/W[0]/W[0];
    // double C = W[n]/rho - 0.5*W[n-1]*W[n-1]/rho/rho;
    double C = W[n]/rho - 0.5*W[n-1]*W[n-1]/W[0]/rho;
    static double T;
    static double T_last;
    static double F;

    // T = thermo::T[idx];
    T = 300;

    double r = thermo::r_mix_comp(comp);

    double h;

    do
    {
        T_last = T;
        F = 0;

        h = thermo::enthalpy(T,comp);

        for(int i = 0; i < n_comp; i++)
        {
            if(comp[i] > thershold_comp)
            {
                F += comp[i]*h;
            }
        }

        F += -r*T - C;

        T = T - ( F )/dF(comp,r,T);

    } while(std::abs(T-T_last) > 1);

    return T;
}

double thermo::pressure(int i, std::vector<double> const& W, std::vector<double> const& comp)
{
    static int mom_idx = W.size()-2;
    static double rho_g;
    static double rho;
    rho_g = density(W);
    rho = W[0];

    // return (thermo::enthalpy(i,W) + 0.5*W[mom_idx]*W[mom_idx]/W[0]/W[0])*W[0] - W.back();
    // return (thermo::enthalpy(i,W) + 0.5*W[mom_idx]*W[mom_idx]/rho/rho)*rho - W.back();
    return rho_g*thermo::enthalpy(i,W) + 0.5*W[mom_idx]*W[mom_idx]/rho - W.back();
}

void thermo::update(std::vector<std::vector<double>> const& W)
{
    // napocitat tlaky a teploty

    static std::vector<double> comp(n_comp);


    for(int i = 0; i < int(W.size()); i++)
    {
        // std::cout << i << "\n";
        thermo::composition(comp,W[i]);
        thermo::T[i] = thermo::temp_new(i,comp,W[i]);
        thermo::p[i] = thermo::pressure(i,W[i],comp);
    }
}