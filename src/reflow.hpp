#pragma once
#include <vector>
#include <random>
#include "mesh.hpp"
#include "variables.hpp"
#include "solver.hpp"
#include "thermodynamics.hpp"
#include "particle.hpp"

class reflow
{
    public:
    variables var;
    mesh msh;
    std::vector<particle> particles;

    void(*left_boundary)(variables&,mesh&,std::vector<double>&);
    void(*right_boundary)(variables&,mesh&,std::vector<double>&);

    std::vector<double> left_values;
    std::vector<double> right_values;

    int n_dt = 2;
    int n_res = 5000;

    volatile double max_res = 50;
    
    
    reflow(variables& _var, mesh& _msh);
    reflow(int N, int N_var, std::vector<double> const& init);
    reflow(int N, double from, double to, int N_var, std::vector<double> const& init);
    reflow(int N, int N_var, std::vector<std::vector<double>> const& init);

    void apply_heat_source(double Q, double x_from, double x_to);

    void spline_geometry(std::vector<std::vector<std::vector<double>>> curves, int n);
    void bump_geometry();

    void set_boundary(void(*left)(variables&,mesh&,std::vector<double>&), void(*right)(variables&,mesh&,std::vector<double>&));
    void set_boundary(void(*left)(variables&,mesh&,std::vector<double>&), std::vector<double> _left_values,
                      void(*right)(variables&,mesh&,std::vector<double>&), std::vector<double> _right_values);


    void spawn_particles(int n_particles);
    void export_particles(std::vector<particle>& particles);

    void solve();
};