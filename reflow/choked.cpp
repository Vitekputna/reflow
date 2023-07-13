#include <iostream>
#include <vector>
#include <memory>

#include "../src/reflow.hpp"
#include "../src/boundary_cond.hpp"
#include "../src/initial_cond.hpp"
#include "../src/params.hpp"

const auto init_comp = std::vector<double>{1};
double T0 = 300;
double p0 = 1e6;
double p2 = 101325;
double md = 200;

int main(int argc, char** argv)
{
    // Sim object
    reflow Simulation;

    // Mesh generation
    std::vector<std::vector<std::vector<double>>> curves;
    std::vector<std::vector<double>> curve = {{0,1},{1.5,1},{1,0},{1,0}};
    curves.push_back(curve);
    curve = {{1.5,1},{2,0.5},{1,0},{1,0}};
    curves.push_back(curve);
    curve = {{2,0.5},{2.5,1},{1,0},{1,0}};
    curves.push_back(curve);
    curve = {{2.5,1},{4,1},{1,0},{1,0}};
    curves.push_back(curve);

    Simulation.refine_mesh(std::vector<std::vector<double>>{{0,4,500}});
    Simulation.spline_geometry(curves,100);

    Simulation.msh.export_to_file();

    // Species
    Simulation.add_specie(287,1.3,29,air_cp,air_k,air_mu);          //Air

    // Initial conditions
    double rho = p0/(thermo::r_mix_comp(init_comp)*T0);
    double u = 0;
    double momentum = rho*u;
    double energy = rho*thermo::enthalpy(T0,init_comp) -p0 + rho*u*u/2;

    // Initial conditions
    Simulation.initial_conditions(std::vector<double>{rho,momentum,energy});

    // Boundary functions and values
    Simulation.add_boundary_function(boundary::pressure_inlet,std::vector<double>{p0,300,1});
    Simulation.add_boundary_function(boundary::subsonic_outlet,std::vector<double>{0.97*p0});

    // Solving
    Simulation.solve(1,1000,0.02);

    // Exporting all lagrangian particles
    Simulation.var.export_to_file(Simulation.msh,Simulation.par_man.particles);

    return 0;
}
