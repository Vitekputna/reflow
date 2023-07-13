#include <iostream>
#include <vector>
#include <memory>

#include "../src/reflow.hpp"
#include "../src/boundary_cond.hpp"
#include "../src/initial_cond.hpp"
#include "../src/params.hpp"

const auto init_comp = std::vector<double>{0.99,0,0.01};

double T0 = 300;

double p0 = 101325;
double p2 = 101325;
double m_A = 100;
double m_W = 0.01;

int main(int argc, char** argv)
{
    // Sim object
    reflow Simulation;

    // Mesh generation
    Simulation.refine_mesh(std::vector<std::vector<double>>{{0,5,2000}});
    Simulation.msh.constant_area(0.002);
    Simulation.msh.export_to_file();

    // Species
    Simulation.add_specie(287,1.3,29,air_cp,air_k,air_mu);                //Air
    Simulation.add_specie(138,1.13,60,fuel_cp,fuel_k,fuel_mu);            //Dummy
    Simulation.add_specie(461,1.4,18,steam_cp,steam_k,steam_mu);          //Water

    // set the fuel properties
    thermo::species[2].h_vap = 2257e3;
    thermo::species[2].T_ref = 373.15;
    thermo::species[2].p_ref = 101325;
    thermo::species[2].rho_liq = 1000;
    thermo::species[2].C = 4182;

    // Initial conditions
    double rho = p0/(thermo::r_mix_comp(init_comp)*T0);
    double u = m_A/(rho*Simulation.msh.A[0]);
    double momentum = rho*u;
    double energy = rho*thermo::enthalpy(T0,init_comp) -p0 + rho*u*u/2;

    // Initial conditions
    Simulation.initial_conditions(std::vector<double>{rho,0,0,momentum,energy});

    // Lagrangian particles manager and specie init
    Simulation.init_particles(1e6,1e3,1000);
    Simulation.add_lagrangian_mono_particles(2,m_W,1000,50.05e-6 ,-0.02,0.5*u,300,300,1e5,1e3);

    // Boundary functions and values
    Simulation.add_boundary_function(boundary::mass_flow_inlet,std::vector<double>{m_A,300,0.99,0,0.01}); 
    Simulation.add_boundary_function(boundary::subsonic_outlet,std::vector<double>{p2});

    // Solving
    Simulation.solve(1,1000,0.1);

    // Exporting all lagrangian particles
    Simulation.var.export_to_file(Simulation.msh,Simulation.par_man.particles);
    Simulation.export_particles(Simulation.par_man.particles);

    return 0;
}
