#include <iostream>
#include <vector>
#include <memory>

#include "../src/reflow.hpp"
#include "../src/boundary_cond.hpp"
#include "../src/initial_cond.hpp"

std::vector<double> prod_cp = {9.95096576e+02, 4.87796972e-01, -1.33106793e-04, 1.32096919e-08, 3.86499580e-13, -9.79511577e-17};
std::vector<double> fuel_cp = {2.95260582e+02, 4.95204053e+00, -2.60871236e-03, 7.03837556e-07, -9.39772178e-11, 4.90497393e-15};
std::vector<double> oxi_cp =  {6.27400878e+02, 1.09090162e+00, -6.21904849e-04, 1.77914259e-07, -2.46557076e-11, 1.31958533e-15};

const auto init_comp = std::vector<double>{0,1,0};
double T0 = 300;
double p2 = 101325;
double md = 100;

int main(int argc, char** argv)
{
    // Sim object
    reflow Simulation;

    // Mesh generation
    Simulation.refine_mesh(std::vector<std::vector<double>>{{0,0.5,1000}});
    Simulation.msh.constant_area(0.002);

    // Species
    Simulation.add_specie(311.39,1.30,26.7,prod_cp);        //Products
    Simulation.add_specie(188,1.31,44,oxi_cp);              //Oxydizer
    Simulation.add_specie(138,1.13,60,fuel_cp);             //Fuel

    // Initial conditions
    double rho = p2/(thermo::r_mix_comp(init_comp)*T0);
    double u = md/(rho*Simulation.msh.A[0]);
    double momentum = rho*u;
    double energy = rho*thermo::enthalpy(T0,init_comp) -p2 + rho*u*u/2;

    // Initial conditions
    Simulation.initial_conditions(std::vector<double>{rho,0,rho,momentum,energy});

    // Lagrangian particles manager and specie init
    Simulation.init_particles(1e6,1e3,1000);
    // Simulation.add_lagrangian_mono_particles(2,20,700,1e-4,0,15,300,300,1e5,1e3);
    Simulation.add_lagrangian_unif_particles(2,20,700,1e-4,1.5e-4,0,15,300,300,1e5,1e-3);

    // Boundary functions and values
    Simulation.add_boundary_function(boundary::subsonic_inlet,std::vector<double>{md,300,0,0,1});
    Simulation.add_boundary_function(boundary::subsonic_outlet,std::vector<double>{p2});

    // Solving
    Simulation.solve(0.1,1000,0.9);

    // Exporting all lagrangian particles
    Simulation.export_particles(Simulation.par_man.particles);

    return 0;
}
