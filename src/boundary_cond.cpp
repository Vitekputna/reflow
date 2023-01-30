#include "boundary_cond.hpp"
#include "thermodynamics.hpp"

#include <iostream>

extern double kappa;

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

    double e = values[0]/(values[1]-1) + 0.5*var.W.back()[var.mom_idx]*var.W.back()[var.mom_idx]/var.W.back()[0];
    var.W.back()[var.eng_idx] = e;
}

// values = (md,T,r,kappa,Y0,Y1...)
void boundary::subsonic_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    // přenos tlaku
    
    double p = thermo::pressure(var.W[1],kappa);
    var.W[0][0] = p/values[2]/values[1];

    for(auto idx = 1; idx <= var.N_comp-1; idx++)
    {
        var.W[0][idx] = var.W[0][0]*values[idx+4];
    }

    var.W[0][var.mom_idx] = values[0]/msh.A[0];
    var.W[0][var.eng_idx] = p/(values[3]-1) + 0.5*var.W[0][var.mom_idx]*var.W[0][var.mom_idx]/var.W[0][0];
}