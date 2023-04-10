#include "boundary_cond.hpp"
#include "thermodynamics.hpp"

#include <iostream>
#include <cmath>
#include <vector>

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

    double rho_gas = p/(r*T_gas);

    // Species fractions
    for(auto idx = 1; idx <= var.N_comp-1; idx++)
    {
        var.W[0][idx] = rho_gas*comp[idx];
    }

    // Droplet mass fractions
    double dm,r_drop,md_frac;
    double droplet_total_mf = 0;

    for(int i = 0; i < N; i++)
    {  
        r_drop = values[2*i+2+5];
        md_frac = values[2*i+1+5];

        dm = 4*M_PI*pow(r_drop,3)*rho_cond/3;

        droplet_total_mf += rho_gas*md_frac/md_gas;

        var.W[0][var.N_comp+i*2+1] = rho_gas*md_frac/md_gas;
        var.W[0][var.N_comp+i*2] = var.W[0][var.N_comp+i*2+1]/dm;
    }

    // total mass fraction
    var.W[0][0] = rho_gas+droplet_total_mf;

    // total momentum
    double u = md_gas/msh.A[0]/rho_gas;

    var.W[0][var.mom_idx] = (rho_gas+droplet_total_mf)*u;

    // total energy
    var.W[0][var.eng_idx] = rho_gas*thermo::enthalpy(T_gas,comp) + 0.5*u*u*var.W[0][0] - p;
}

std::vector<double> boundary::flow_with_droplets(double md_gas, double T, std::vector<double> comp, int N_frac, double md_cond, double rho, double mean_r, double var_r)
{
    std::vector<double> return_vec = {md_gas,T};

    for(auto const& Y : comp)
    {
        return_vec.push_back(Y);
    }

    std::vector<double> drop_vec = discretize_distribution(normal_distribution,md_cond,rho,mean_r,var_r,N_frac);

    return_vec.insert(return_vec.end(),drop_vec.begin(),drop_vec.end());

    return return_vec;
}

// (mean,var)
double boundary::normal_distribution(double x, std::vector<double> values)
{
    double mean = values[0];
    double var = values[1];

    return (1/(var*sqrt(2*M_PI)))*exp(-pow(x -mean,2)/2/pow(var,2));
}

// (mean,var,rho,N)
double boundary::mass_distribution(double x, double(*distribution)(double,std::vector<double>),std::vector<double> values)
{
    std::vector<double> params = {values[0],values[1]};

    double rho = values[2];
    double N = values[3];

    return distribution(x,params)*(4*M_PI*pow(x,3)/3*rho*N);   
}

std::vector<double> boundary::trapz(double(*func)(double,std::vector<double>),std::vector<double> values, double x_from, double x_to, int N)
{
    int N_points = N+1;

    double delta_r = (x_to-x_from)/N;

    double tolerance = 1e-12;

    std::vector<double> x;

    for(double r = x_from; r-x_to < tolerance; r+=delta_r)
    {
        x.push_back(r);
    }

    double x_middle, y_middle;
    double den = 0;
    double num = 0;

    double integral = 0;

    for(int i = 0; i < N; i++)
    {
        x_middle = 0.5*(x[i+1] + x[i]);
        
        y_middle = mass_distribution(x_middle,func,values);

        integral += y_middle*(x[i+1]-x[i]);

        den += y_middle;
        num += y_middle*x_middle;
    }

    double center = num/den;

    return std::vector<double>{integral,center};
}

std::vector<double> boundary::discretize_distribution(double(*distribution)(double,std::vector<double>), double md, double rho, double mean, double var, int N_intervals)
{
    std::cout << "##########################################\n";
    std::cout << "Discretizing droplet distribution\n";

    std::cout << "Mean radius:\t" << mean << "\n";
    std::cout << "Variance:\t" << var << "\n";
    std::cout << "Mass flux:\t" << md << "\n";
    std::cout << "Density:\t" << rho << "\n";
    std::cout << "N intervals:\t" << N_intervals << "\n";

    // compute parameters of distribution
    double nu3 = 3*mean*pow(var,2) + pow(mean,3);           // Compute third uncentered moment of normal distribution

    double N = 3*md/(4*M_PI*rho*nu3);                       // Compute number of droplets based on mean radius

    std::cout << "Number of droplets:\t" << N << "\n";

    double x_from = mean - 5*var;                                   // minimum
    double x_to = mean + 5*var;                                     // maximum

    std::cout << "r min:\t" << x_from << "\n";
    std::cout << "r max:\t" << x_to << "\n";

    if(x_from < 0) x_from = 0;

    // construct linspace
    int N_points = N_intervals+1;

    double delta_x = (x_to-x_from)/N_intervals;

    double tolerance = 1e-8;

    std::vector<double> x_vector;

    for(double x = x_from; x-x_to < tolerance; x+=delta_x)
    {
        x_vector.push_back(x);
    }

    // Discretize distribution
    double total_integral = 0;

    std::vector<double> params = {mean,var,rho,N};
    std::vector<double> return_vec = {(double)N_intervals};

    for(int i = 0; i < N_intervals; i++)
    {
        auto res = trapz(distribution,params,x_vector[i],x_vector[i+1],51);

        return_vec.push_back(res[0]);
        return_vec.push_back(res[1]);

        std::cout << i << " r:  " << res[1] << "\tmd:  " << res[0] << "\n";

        total_integral += res[0];
    }

    return_vec.push_back(rho);

    std::cout << "Discretized mass flux:\t" << total_integral << "\n";
    std::cout << "##########################################\n";

    return return_vec;
}