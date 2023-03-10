#include "boundary_cond.hpp"
#include "thermodynamics.hpp"

#include <iostream>

void boundary::_set_value_l(variables& var, int idx, double value)
{
    var.W[0][idx] = value;
}

void boundary::_set_value_r(variables& var, int idx, double value)
{
    var.W.back()[idx] = value;
}

void boundary::_zero_gradient_l(variables& var, int idx, double value)
{
    var.W[0][idx] = var.W[1][idx];
}

void boundary::_zero_gradient_r(variables& var, int idx, double value)
{
    var.W.back()[idx] = var.W.rbegin()[1][idx];
}

void boundary::zero_gradient_l(variables& var, mesh& msh, std::vector<double>& values)
{
    for(int k = 0; k < var.N_var; k++)
    {
        _zero_gradient_l(var,k,0);
    }
}

void boundary::zero_gradient_r(variables& var, mesh& msh, std::vector<double>& values)
{
    for(int k = 0; k < var.N_var; k++)
    {
        _zero_gradient_r(var,k,0);
    }
}

void boundary::set_value_l(variables& var, mesh& msh, std::vector<double>& values)
{
    for(int k = 0; k < var.N_var; k++)
    {
        _set_value_l(var,k,values[k]);
    }
}

// values = (md,p,T,r,kappa,composition_vector)
void boundary::mass_flow_inlet(variables& var, mesh& msh, std::vector<double>& values) 
{
    var.W[0][0] = values[1]/(values[2]*values[3]);
    var.W[0][1] = values[0]/msh.A[0];
    var.W[0][2] = values[1]/(values[4]-1) + 0.5*(var.W[0][1]*var.W[0][1])/var.W[0][0];
}

// values = (p,kappa)
void boundary::subsonic_outlet(variables& var, mesh& msh, std::vector<double>& values)
{
    for(auto idx = 0; idx < var.N_var-1; idx++)
    {
        var.W.back()[idx] = var.W.rbegin()[1][idx];
    }

    double r = thermo::r_mix(var.W.rbegin()[1]);

    double T = values[0]/r/var.W.back()[0];

    static std::vector<double> comp(var.N_comp);
    thermo::composition(comp,var.W.rbegin()[1]);

    double e = (thermo::enthalpy(T,comp) + 0.5*var.W.back()[var.mom_idx]*var.W.back()[var.mom_idx]/var.W.back()[0]/var.W.back()[0]/2)*var.W.back()[0] - values[0]; 

    var.W.back()[var.eng_idx] = e;
}

// values = (md,T,Y0,Y1...)
void boundary::subsonic_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    // p??enos tlaku
    // double p1 = thermo::p[1];
    double p2 = thermo::p[2];

    // double p = 2*p1-p2; //lin extrapolation of 1 and 2
    // double p = 0.5*(p1 + p2); //avg of 1 and 2
    double p = p2; // val of 1
    
    const std::vector<double> comp = {values[2],values[3],values[4]};

    double r = thermo::r_mix_comp(comp);

    var.W[0][0] = p/r/values[1];
    

    for(auto idx = 1; idx <= var.N_comp-1; idx++)
    {
        var.W[0][idx] = var.W[0][0]*comp[idx];
    }

    var.W[0][var.mom_idx] = values[0]/msh.A[0];

    var.W[0][var.eng_idx] = (thermo::enthalpy(values[1],comp) + 0.5*var.W[0][var.mom_idx]*var.W[0][var.mom_idx]/var.W[0][0]/var.W[0][0])*var.W[0][0] - p;
}

void boundary::supersonic_outlet(variables& var, mesh& msh, std::vector<double>& values)
{
    for(int i = 0; i < var.N_var; i++)
    {
        var.W.back()[i] = 0.5*(var.W.rbegin()[2][i] + var.W.rbegin()[1][i]);
    }
}