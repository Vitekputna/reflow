#pragma once

#include "mesh.hpp"
#include "variables.hpp"
#include "solver.hpp"
#include "thermodynamics.hpp"

class reflow
{
    variables var;
    mesh msh;
    
    public:
    reflow(variables& _var, mesh& _msh);
    reflow(int N, int N_var, std::vector<double> init);
    void solve();
};