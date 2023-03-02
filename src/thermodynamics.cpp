#include "thermodynamics.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include "specie.hpp"

std::vector<specie> thermo::species;

std::vector<double> thermo::T, thermo::p;

int thermo::n_comp = 0;
double thermo::thershold_comp = 1e-4;

thermo::thermo() {}

void thermo::init(int n)
{
    thermo::p = std::vector<double>(n,0.0);
    thermo::T = std::vector<double>(n,0.0);
}

// Using ideal gas law
double thermo::pressure(std::vector<double> const& W)
{
    static int mom_idx = W.size()-2;
    double r = thermo::r_mix(W);
    double cp = thermo::cp_mix(W);

    return r/(r-cp)*(0.5*W[mom_idx]*W[mom_idx]/W[0] - W[mom_idx+1]);    
}

double thermo::speed_of_sound(int i, std::vector<double> const& W)
{
    double kappa = thermo::kappa_mix(W);
    return sqrt(kappa*thermo::p[i]/W[0]);
}

double thermo::temperature(std::vector<double> const& W)
{
    static std::vector<double> comp(n_comp);
    thermo::composition(comp,W);

    return temp_new(0,comp,W);
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
    static std::vector<double> comp(n_comp); 
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
    static std::vector<double> comp(n_comp);
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

double thermo::cp_mix_comp(std::vector<double>& comp, double T)
{
    double cp = 0;

    for(int i = 0; i < n_comp; i++)
    {
        cp += comp[i]*thermo::species[i].cp(T);
    }

    return cp;
}

double thermo::dF(std::vector<double> const& comp, double T)
{
    double df = 0;

    for(int i = 0; i < n_comp; i++)
    {
        if(comp[i]>thershold_comp)
        {
            df += comp[i]*(species[i].r - species[i].b - 2*species[i].c*T - 3*species[i].d*T*T - 4*species[i].e*std::pow(T,3) - 5*species[i].f*std::pow(T,4)- 6*species[i].g*std::pow(T,5)); 
        }
    }

    return df;
    // return species[1].r - species[1].b - 2*species[1].c*T - 3*species[1].d*T*T - 4*species[1].e*T*T*T - 5*species[1].f*T*T*T*T - 6*species[1].g*T*T*T*T*T;
}

double thermo::temp_new(int idx, std::vector<double> const& comp, std::vector<double> const& W)
{
    double C = (0.5*W[3]*W[3]/W[0] - W[4])/W[0];
    static double T;
    static double T_last;
    static double F;

    T = thermo::T[idx];

    do
    {
        T_last = T;
        F = 0;

        for(int i = 0; i < n_comp; i++)
        {
            if(comp[i] > thershold_comp)
            {
                F += comp[i]*(thermo::species[i].r*T - thermo::species[i].a - thermo::species[i].b*T - thermo::species[i].c*std::pow(T,2) - thermo::species[i].d*std::pow(T,3) - thermo::species[i].e*std::pow(T,4)
                    -thermo::species[i].f*std::pow(T,5) - thermo::species[i].g*std::pow(T,6) - C);
            }
        }

        T = T - ( F )/dF(comp,T);

    } while(std::abs(T-T_last) > 10);

    return T;
}

double thermo::pressure(double T, std::vector<double> const& W)
{
    static std::vector<double> comp(n_comp);
    thermo::composition(comp,W);

    static int mom_idx = W.size()-2;
    double r = thermo::r_mix(W);
    double cp = thermo::cp_mix_comp(comp,T);

    return r/(r-cp)*(0.5*W[mom_idx]*W[mom_idx]/W[0] - W[mom_idx+1]);    
}

void thermo::update(std::vector<std::vector<double>> const& W)
{
    // napocitat tlaky a teploty

    static std::vector<double> comp(n_comp);

    for(int i = 0; i < W.size(); i++)
    {
        thermo::composition(comp,W[i]);
        thermo::T[i] = thermo::temp_new(i,comp,W[i]);
        thermo::p[i] = thermo::pressure(thermo::T[i],W[i]);
    }
}