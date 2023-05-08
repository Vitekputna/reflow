#pragma once
#include "variables.hpp"
#include "mesh.hpp"

namespace boundary
{
    // abstractions
    void mass_flow_inlet(variables& var, mesh& msh, std::vector<double>& values);
    void quiscent_droplet_inlet(variables& var, mesh& msh, std::vector<double>& values);
    void active_drop_inlet(variables& var, mesh& msh, std::vector<double>& values);
    void active_thermal_drop_inlet(variables& var, mesh& msh, std::vector<double>& values);

    void subsonic_outlet(variables& var, mesh& msh, std::vector<double>& values);
    void supersonic_outlet(variables& var, mesh& msh, std::vector<double>& values);

    // initializations
    std::vector<double> quiscent_droplets(double(*distribution)(double,std::vector<double>),int N_frac, double md_cond, double rho, double mean_r, double var_r);
    std::vector<double> active_droplets(double(*distribution)(double,std::vector<double>),int N_frac, double md_cond, double rho, double u_drop, double mean_r, double var_r);
    std::vector<double> active_thermal_droplets(double(*distribution)(double,std::vector<double>),int N_frac, double md_cond, double rho, double T, double u_drop, double mean_r, double var_r);

    //helper functions
    double normal_distribution(double x, std::vector<double> values);
    double mass_distribution(double x, double(*distribution)(double,std::vector<double>),std::vector<double> values);
    std::vector<double> trapz(double(*func)(double,std::vector<double>),std::vector<double> values,double x_from, double x_to, int N);
    std::vector<double> discretize_distribution(double(*distribution)(double,std::vector<double>), double md, double rho, double mean, double var, int N_intervals);
}