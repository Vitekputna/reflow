#pragma once
#include "variables.hpp"
#include "mesh.hpp"

namespace boundary
{
    // fundamental
    void _zero_gradient_r(variables& var, int idx, double value);
    void _set_gradient_r(variables& var, int idx, double value);
    void _set_value_r(variables& var, int idx, double value);

    void _zero_gradient_l(variables& var, int idx, double value);
    void _set_gradient_l(variables& var, int idx, double value);
    void _set_value_l(variables& var, int idx, double value);

    // abstractions
    void zero_gradient_l(variables& var, mesh& msh, std::vector<double>& values);
    void zero_gradient_r(variables& var, mesh& msh, std::vector<double>& values);
    void set_value_l(variables& var, mesh& msh, std::vector<double>& values);
    void mass_flow_inlet(variables& var, mesh& msh, std::vector<double>& values);
    void subsonic_inlet(variables& var, mesh& msh, std::vector<double>& values);
    void subsonic_outlet(variables& var, mesh& msh, std::vector<double>& values);
    void supersonic_outlet(variables& var, mesh& msh, std::vector<double>& values);

    // drop + flow
    void mass_flow_inlet_with_droplets(variables& var, mesh& msh, std::vector<double>& values);
    std::vector<double> flow_with_droplets(double md_gas, double T, std::vector<double> comp, int N_frac, double md_cond, double rho, double mean_r, double var_r);

    //helper functions
    double normal_distribution(double x, std::vector<double> values);
    double mass_distribution(double x, double(*distribution)(double,std::vector<double>),std::vector<double> values);
    std::vector<double> trapz(double(*func)(double,std::vector<double>),std::vector<double> values,double x_from, double x_to, int N);
    std::vector<double> discretize_distribution(double(*distribution)(double,std::vector<double>), double md, double rho, double mean, double var, int N_intervals);
}