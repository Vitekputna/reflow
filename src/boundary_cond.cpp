#include "boundary_cond.hpp"
#include "thermodynamics.hpp"

#include <iostream>
#include <cmath>

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

    double e = (thermo::enthalpy(T,comp) + 0.5*var.W.back()[var.mom_idx]*var.W.back()[var.mom_idx]/var.W.back()[0]/var.W.back()[0])*var.W.back()[0] - values[0]; 

    var.W.back()[var.eng_idx] = e;
}

// values = (md,T,Y0,Y1...)
void boundary::subsonic_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    // přenos tlaku
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

// values = (md,r,rho)
void boundary::quiscent_droplet_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    double dm = 4/3*M_PI*std::pow(values[1],3)*values[2];

    var.W[0][4] = values[0]/(msh.A[0]*(var.W[0][var.mom_idx]/var.W[0][0]));
    var.W[0][3] = var.W[0][4]/dm;
}

// values = (N,md1,r1,md2,r2...,rho) N = number of {md,r} pairs
void boundary::quiscent_droplets_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    int N = values[0];

    double dm;
    double rho = values.back();

    for(int i = 0; i < N; i++)
    {  
        dm = 4/3*M_PI*std::pow(values[2*i+2],3)*rho;

        var.W[0][var.N_comp+i*2+1] = values[2*i+1]/(msh.A[0]*(var.W[0][var.mom_idx]/var.W[0][0]));
        var.W[0][var.N_comp+i*2] = var.W[0][var.N_comp+i*2+1]/dm;

        var.W[0][0] += var.W[0][var.N_comp+i*2+1];
        var.W[0][var.mom_idx] += var.W[0][var.mom_idx]/var.W[0][0]*var.W[0][var.N_comp+i*2+1];
        var.W[0][var.eng_idx] += 0.5*var.W[0][var.N_comp+i*2+1]*(var.W[0][var.mom_idx]/var.W[0][0])*(var.W[0][var.mom_idx]/var.W[0][0]);
    }
}

// values = (md_gas,T,Y0,Y1,Y2,N,md1,r1,md2,r2...,rho) N = number of {md,r} pairs
void boundary::mass_flow_inlet_with_droplets(variables& var, mesh& msh, std::vector<double>& values)
{
    // double p = (thermo::p[1] + thermo::p[2])/2;
    double p = 2*thermo::p[1] - thermo::p[2]; // linear extrapolation

    const std::vector<double> comp = {values[2],values[3],values[4]};
    double r = thermo::r_mix_comp(comp);
    double md_gas = values[0];
    double T_gas = values[1];

    int N = values[5];

    double rho_cond = values.back();

    double rho_gas = p/r/T_gas;

    double droplet_total_mf = 0;

    // Species fractions
    for(auto idx = 1; idx <= var.N_comp-1; idx++)
    {
        var.W[0][idx] = rho_gas*comp[idx];
    }

    // Droplet mass fractions
    double dm,r_drop,md_frac;

    for(int i = 0; i < N; i++)
    {  
        r_drop = values[2*i+2+5];
        md_frac = values[2*i+1+5];

        dm = 4/3*M_PI*std::pow(r_drop,3)*rho_cond;

        droplet_total_mf += md_frac*rho_gas/md_gas;

        var.W[0][var.N_comp+i*2+1] = md_frac*rho_gas/md_gas; 
        var.W[0][var.N_comp+i*2] = var.W[0][var.N_comp+i*2+1]/dm;
    }

    // total mass fraction
    var.W[0][0] = rho_gas+droplet_total_mf;

    // total momentum
    double u = md_gas/msh.A[0]/rho_gas;

    var.W[0][var.mom_idx] = (rho_gas+droplet_total_mf)*u;
    // var.W[0][var.mom_idx] = (values[0]+values[6])/msh.A[0];

    // total energy
    var.W[0][var.eng_idx] = rho_gas*thermo::enthalpy(T_gas,comp) + 0.5*u*u*var.W[0][0] - p;
}


std::vector<double> boundary::normal_distribution(int N_fracs, double mass_flow, double r_mean,double r_var)
{
    double r_max = r_mean + 5*r_var;
    double r_min = r_mean - 5*r_var;

    if(r_min < 0) r_min = 0;

    // Construct linspace from r_min to r_max with N intervals
    int N_points = N_fracs+1;

    double delta_r = (r_max-r_min)/N_fracs;

    double tolerance = 1e-8;

    std::vector<double> r_vector;

    for(double r = r_min; r-r_max < tolerance; r+=delta_r)
    {
        std::cout << r << "\n";
        r_vector.push_back(r);
    }

    // Compute disrete normal distribution
    double f,r;

    for(int i = 0; i < N_fracs; i++)
    {

        for(int j = )

        r = (r_vector[i+1] + r_vector[i])/2;       //mean r inside interval

        f = (1/(r_var*sqrt(2*M_PI)))*exp(-pow(r -r_mean,2)/2/pow(r_var,2));     //


    }

    return std::vector<double>{};
}

// (x,mean,var)
double boundary::normal_distribution(std::vector<double> values)
{
    double x = values[0];
    double mean = values[1];
    double var = values[2];

    return (1/(var*sqrt(2*M_PI)))*exp(-pow(x -mean,2)/2/pow(var,2));
}

// (x,mean,var,rho,N)
double boundary::mass_distribution(double(*distribution)(std::vector<double>),std::vector<double> values)
{
    std::vector<double> params = {values[0],values[1],values[2]};

    double x = values[0];
    double rho = values[3];
    double N = values[4];

    return normal_distribution(params)*(4*M_PI*pow(x,3)/3*rho*N);   
}

double boundary::trapz(double(*func)(std::vector<double>),std::vector<double> values, double x_from, double x_to, int N)
{
    int N_points = N+1;

    double delta_r = (x_to-x_from)/N;

    double tolerance = 1e-8;

    std::vector<double> x;

    for(double r = x_from; r-x_to < tolerance; r+=delta_r)
    {
        x.push_back(r);
    }

    
}