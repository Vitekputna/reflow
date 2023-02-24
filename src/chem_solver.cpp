#include <iostream>
#include "chem_solver.hpp"
#include "thermodynamics.hpp"

reaction::reaction(){}

reaction::reaction(std::vector<int> R_idx,std::vector<int> P_idx, std::vector<double> R_frac, std::vector<double> P_frac, double dH)
: reactant_idx{R_idx}, product_idx{P_idx}, reactant_fraction{R_frac}, product_fraction{P_frac}, Heat_relased{dH}
{
    n_reaktants = R_idx.size();
    n_products = P_idx.size();
}

chem_solver::chem_solver(){}

void chem_solver::add_reaction(reaction& R)
{
    reactions.push_back(R);
}

void chem_solver::solve(double dt,std::vector<std::vector<double>>& res, variables& var, mesh const& msh)
{
    double K,V;

    std::vector<double> comp(var.N_comp,0.0);

    for(int i = 1; i < var.N-1; i++)
    {
        // V = msh.A[i]*(msh.xf[i] - msh.xf[i-1]);
        // thermo::composition(comp,var.W[i]);

        for(auto const& reaction : reactions)
        {
            K = 1e50;
            
            for(int idx = 0; idx < reaction.n_reaktants; idx++)
            {
                // K = std::min(K, comp[reaction.reactant_idx[idx]]/reaction.reactant_fraction[idx]);
                K = std::min(K, var.W[i][reaction.reactant_idx[idx]]/reaction.reactant_fraction[idx]*reaction.reactant_fraction[0]);
            }

            K = std::max(0.0,K); 

            for(int idx = 0; idx < reaction.n_reaktants; idx++)
            {
                res[i][reaction.reactant_idx[idx]] -= K*reaction.reactant_fraction[idx]/reaction.reactant_fraction[0]/dt;
            }

            for(int idx = 0; idx < reaction.n_products; idx++)
            {
                res[i][reaction.product_idx[idx]] += K*reaction.product_fraction[idx]/reaction.reactant_fraction[0]/dt;
            }

            res[i][var.eng_idx] += K*reaction.Heat_relased/dt;
        }
    }
}