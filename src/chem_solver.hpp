#pragma once
#include <vector>
#include "variables.hpp"
#include "mesh.hpp"

struct reaction
{
    std::vector<int> reactant_idx;
    std::vector<int> product_idx;
    std::vector<double> reactant_fraction;
    std::vector<double> product_fraction;

    double Heat_relased = 0;
    int n_reaktants, n_products;

    reaction();
    reaction(std::vector<int> R_idx,std::vector<int> P_idx, std::vector<double> R_frac, std::vector<double> P_frac, double dH);
};

class chem_solver
{
    private:
    std::vector<reaction> reactions;

    void infinite_rate_chemistry(double dt,std::vector<std::vector<double>>& res, variables& var, mesh const& msh);
    void finite_rate_chemistry(double dt,std::vector<std::vector<double>>& res, variables& var, mesh const& msh);

    public:
    chem_solver();

    void add_reaction(reaction& R);
    void solve(double dt,std::vector<std::vector<double>>& res, variables& var, mesh const& msh);
};