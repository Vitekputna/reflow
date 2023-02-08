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

std::vector<double> init::flow(int N_var, double p, double T, double u, std::vector<double> comp)
{
    double r = thermo::r_mix_comp(comp);
    double kappa = thermo::kappa_mix_comp(comp);

    int n_comp = N_var - 2;
    double rho = p/r/T;

    std::vector<double> res = {rho};

    for(auto i = 1; i < n_comp; i++)
    {
        res.push_back(rho*comp[i]);
    }

    res.push_back(rho*u);
    res.push_back(p/(kappa-1) + 0.5*rho*u*u);

    // std::cout << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << " " << res[4] << "\n";

    return res;
}