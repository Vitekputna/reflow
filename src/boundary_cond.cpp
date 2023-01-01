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

// values = (md,p,T,r,kappa)
void boundary::mass_flow_inlet(variables& var, mesh& msh, std::vector<double>& values) 
{
    var.W[0][0] = values[1]/(values[2]*values[3]);
    var.W[0][1] = values[0]/msh.A[0];
    var.W[0][2] = values[1]/(values[4]-1) + 0.5*(var.W[0][1]*var.W[0][1])/var.W[0][0];
}

// values = (p,kappa)
void boundary::subsonic_outlet(variables& var, mesh& msh, std::vector<double>& values)
{
    var.W.back()[0] = var.W.rbegin()[1][0];
    var.W.back()[1] = var.W.rbegin()[1][1];

    double e = values[0]/(values[1]-1) + 0.5*var.W.back()[1]*var.W.back()[1]/var.W.back()[0];

    var.W.back()[2] = e;
}

// values = (md,T,r,kappa)
void boundary::subsonic_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    // přenos hustoty
    // var.W[0][0] = var.W[1][0];
    // var.W[0][1] = values[0]/msh.A[0];
    // var.W[0][2] = values[1]/(values[2]-1) + 0.5*var.W[0][1]*var.W[0][1]/var.W[0][0];

    // přenos tlaku
    var.W[0][1] = values[0]/msh.A[0];
    double p = thermo::pressure(var.W[1],kappa);
    var.W[0][0] = p/values[2]/values[1];
    var.W[0][2] = p/(values[3]-1) + 0.5*var.W[0][1]*var.W[0][1]/var.W[0][0];

    // první riemannův invariant
    // double K = var.W[1][1]/var.W[1][0] - 2*thermo::speed_of_sound(var.W[1],kappa)/(kappa - 1);
    // // double K2 = std::max(1.0,K + 2*thermo::speed_of_sound(var.W[0],kappa)/(kappa - 1));
    // double K2 = K + 2*thermo::speed_of_sound(var.W[0],kappa)/(kappa - 1);
    // var.W[0][0] = values[0]/(msh.A[0]*(K2));
    // var.W[0][1] = values[0]/msh.A[0];
    // double p0 = var.W[0][0]*values[1]*values[2];
    // var.W[0][2] = p0/(values[3] - 1)+0.5*var.W[0][1]*var.W[0][1]/var.W[0][0];
    //std::cout << K << " " << p0 << "\n";

    // přenos rychlosti
    // var.W[0][1] = values[0]/msh.A[0];
    // var.W[0][0] = var.W[0][1]/var.W[1][0];
    // double p = var.W[0][0]*values[2]*values[1];
    // var.W[0][2] = p/(values[3] - 1) + 0.5*var.W[0][1]*var.W[0][1]/var.W[0][0];
}