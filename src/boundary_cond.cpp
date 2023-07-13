#include "boundary_cond.hpp"
#include "thermodynamics.hpp"

#include <iostream>
#include <cmath>
#include <vector>

// values = (p)
void boundary::subsonic_outlet(variables& var, mesh& msh, std::vector<double>& values)
{
    for(auto idx = 0; idx < var.N_var-1; idx++)
    {
        var.W.back()[idx] = var.W.rbegin()[1][idx];
    }

    const double alfa = 1;

    const double p_flex = 0.5*(thermo::p.rbegin()[1] + thermo::p.rbegin()[2]);
    const double p_fixed = values[0];

    const double p = (1-alfa)*p_flex + alfa*p_fixed;

    const double r = thermo::r_mix(var.W.rbegin()[1]);

    const double rho = thermo::density(var.W.rbegin()[1]);

    const double T = p/r/rho;

    const double u_gas = var.W.back()[var.mom_idx]/var.W.back()[0];

    static std::vector<double> comp(var.N_comp);
    thermo::composition(comp,var.W.rbegin()[1]);

    const double gas_vol_frac = 1-thermo::liquid_fraction(var.W.rbegin()[1]);

    double e = gas_vol_frac*((thermo::enthalpy(T,comp) + 0.5*u_gas*u_gas)*rho - p);

    var.W.back()[var.eng_idx] = e;
}

// values = (md,T,Y0,Y1...)
void boundary::mass_flow_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    // přenos tlaku
    double p1 = thermo::p[1];
    double p2 = thermo::p[2];

    // double p = 2*p1-p2; //lin extrapolation of 1 and 2
    // double p = 0.5*(p1 + p2); //avg of 1 and 2
    // double p = p2; // val of 1
    const double p = p1; // val of 1
    
    const std::vector<double> comp = {values[2],values[3],values[4]};

    const double r = thermo::r_mix_comp(comp);

    var.W[0][0] = p/r/values[1];

    for(auto idx = 1; idx <= var.N_comp-1; idx++)
    {
        var.W[0][idx] = var.W[0][0]*comp[idx];
    }

    var.W[0][var.mom_idx] = values[0]/msh.A[0];

    const double u = var.W[0][var.mom_idx]/var.W[0][0];

    var.W[0][var.eng_idx] = (thermo::enthalpy(values[1],comp) + 0.5*u*u)*var.W[0][0] - p;
}

// values = (p_0,T_0,Y0,Y1...)
void boundary::pressure_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    const double p0 = values[0];
    const double T0 = values[1];

    const std::vector<double> comp = {values[2],values[3],values[4]};
    const double r = thermo::r_mix_comp(comp);
    const double k = thermo::kappa_mix_comp(comp);

    const double M = thermo::mach_number(1,var.W[1]);

    const double A = (1+ pow(M,2)*(k-1)/2);

    const double p = pow(A,-k/(k-1))*p0;
    const double T = pow(A,-1)*T0;

    const double rho = p/(r*T);

    const double c = sqrt(k*r*T);
    const double u = c*M;

    // Total density
    var.W[0][0] = rho;

    // Specie density
    for(auto idx = 1; idx <= var.N_comp-1; idx++)
    {
        var.W[0][idx] = rho*comp[idx];
    }

    // Momentum
    var.W[0][var.mom_idx] = rho*u;

    // Energy
    const double h = thermo::enthalpy(T,comp);
    var.W[0][var.eng_idx] = (h + 0.5*u*u)*rho - p;
}

// (N,md1,r1,md2,r2...,rho) N = number of {md,r} pairs
void boundary::quiscent_droplet_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    const int N = values[0];
    const double md_gas = var.W[0][variables::mom_idx]*msh.A[0];
    const double rho_gas = var.W[0][0];
    const double rho_cond = values.back();
    const double u_gas = var.W[0][variables::mom_idx]/rho_gas;

    double m,r_drop,md_frac;
    double droplet_total_mf = 0;
    
    for(int i = 0; i < N; i++)
    {  
        r_drop = values[2*i+2];
        md_frac = values[2*i+1];

        m = 4*M_PI*pow(r_drop,3)*rho_cond/3;

        droplet_total_mf += rho_gas*md_frac/md_gas;

        var.W[0][var.N_comp+i*2+1] = rho_gas*md_frac/md_gas;
        var.W[0][var.N_comp+i*2] = var.W[0][var.N_comp+i*2+1]/m;
    }

    var.W[0][0] += droplet_total_mf;
    var.W[0][variables::mom_idx] += droplet_total_mf*u_gas;
    var.W[0][variables::eng_idx] += 0.5*droplet_total_mf*u_gas*u_gas;

}

// (N,md1,r1,md2,r2...,rho,u) N = number of {md,r} pairs
void boundary::active_drop_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    const int N = values[0];
    const double md_gas = var.W[0][variables::mom_idx]*msh.A[0];
    const double rho_gas = var.W[0][0];
    const double rho_cond = values.rbegin()[1];
    const double u_cond = values.back();
    const double u_gas = var.W[0][variables::mom_idx]/rho_gas;

    double m,r_drop,md_frac;
    double droplet_total_mf = 0;
    
    for(int i = 0; i < N; i++)
    {  
        r_drop = values[2*i+2];
        md_frac = values[2*i+1];

        m = 4*M_PI*pow(r_drop,3)*rho_cond/3;

        droplet_total_mf += rho_gas*md_frac/md_gas;

        var.W[0][var.N_comp+i*2+1] = md_frac/(msh.A[0]*u_cond);
        var.W[0][var.N_comp+i*2] = var.W[0][var.N_comp+i*2+1]/m;
        var.W[0][variables::drop_mom_idx[i]] = var.W[0][var.N_comp+i*2+1]*u_cond;
    }
}

// (N,md1,r1,md2,r2...,rho,u,T) N = number of {md,r} pairs
void boundary::active_thermal_drop_inlet(variables& var, mesh& msh, std::vector<double>& values)
{
    const int N = values[0];
    const double rho_cond = values.rbegin()[2];
    const double u_cond = values.rbegin()[1];
    const double T_cond = values.back();

    // Dispersed fraction
    double m,r_drop,md_frac;
    double droplet_total_mf = 0;
    
    for(int i = 0; i < N; i++)
    {  
        r_drop = values[2*i+2];

        md_frac = values[2*i+1];

        m = (4*M_PI*pow(r_drop,3)*rho_cond)/3;

        //mass concentraion
        var.W[0][var.N_comp+i*2+1] = md_frac/(msh.A[0]*u_cond);

        //number concentration
        var.W[0][var.N_comp+i*2] = var.W[0][var.N_comp+i*2+1]/m;

        //momentum
        var.W[0][variables::drop_mom_idx[i]] = var.W[0][var.N_comp+i*2+1]*u_cond;

        //energy
        var.W[0][var.N_comp + 3*variables::N_drop_frac + i] = var.W[0][var.N_comp+i*2+1]*thermo::species[2].C*T_cond;
    }

    // Continuous phase correction
    const double gas_vol_frac = 1- thermo::liquid_fraction(var.W[0]);

    const double gas_density = var.W[0][0];

    for(int i = 0; i < variables::N_comp; i++) 
    {
        var.W[0][i] *= gas_vol_frac; 
    }

    const double gas_momentum = var.W[0][variables::mom_idx];
    const double u_gas = gas_momentum/gas_density;
    const double u_gas_corrected = u_gas/gas_vol_frac;

    var.W[0][variables::mom_idx] = gas_vol_frac*gas_density*u_gas_corrected;

    const double gas_energy = var.W[0][variables::eng_idx];

    var.W[0][variables::eng_idx] = (gas_energy - 0.5*gas_density*u_gas*u_gas + 0.5*gas_density*u_gas_corrected*u_gas_corrected)*gas_vol_frac;
}

void boundary::supersonic_outlet(variables& var, mesh& msh, std::vector<double>& values)
{
    for(int i = 0; i < var.N_var; i++)
    {
        var.W.back()[i] = 0.5*(var.W.rbegin()[2][i] + var.W.rbegin()[1][i]);
    }
}

std::vector<double> boundary::quiscent_droplets(double(*distribution)(double,std::vector<double>),int N_frac, double md_cond, double rho, double mean_r, double var_r)
{
    std::vector<double> drop_vec = discretize_distribution(normal_distribution,md_cond,rho,mean_r,var_r,N_frac);
    return drop_vec;
}

std::vector<double> boundary::active_droplets(double(*distribution)(double,std::vector<double>),int N_frac, double md_cond, double rho, double u_drop, double mean_r, double var_r)
{
    std::vector<double> drop_vec = discretize_distribution(normal_distribution,md_cond,rho,mean_r,var_r,N_frac);
    drop_vec.push_back(u_drop);
    return drop_vec;
}

std::vector<double> boundary::active_thermal_droplets(double(*distribution)(double,std::vector<double>),int N_frac, double md_cond, double rho, double T, double u_drop, double mean_r, double var_r)
{
    std::vector<double> drop_vec = discretize_distribution(normal_distribution,md_cond,rho,mean_r,var_r,N_frac);
    drop_vec.push_back(u_drop);
    drop_vec.push_back(T);
    return drop_vec;
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