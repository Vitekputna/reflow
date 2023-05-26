#include <iostream>
#include <vector>
#include <memory>

#include "../src/reflow.hpp"
#include "../src/boundary_cond.hpp"
#include "../src/initial_cond.hpp"
#include "../src/params.hpp"

const auto init_comp = std::vector<double>{0,1,0};
double T0 = 300;
double p2 = 101325;
double md = 100;

// double OF = 1;
// double m_F = md/(OF+1);
// double m_OX = md-m_F;

double OF = 5;
double m_OX = 100;
double m_F = m_OX/OF;

int main(int argc, char** argv)
{
    // Sim object
    reflow Simulation;

    // Mesh generation
    // Simulation.refine_mesh(std::vector<std::vector<double>>{{0,0.1,200},{0.1,5,1000}});
    Simulation.refine_mesh(std::vector<std::vector<double>>{{0,5,1000}});
    Simulation.msh.constant_area(0.002);

    // Simulation.load_old_data("out/",3,0,true,true);

    Simulation.msh.export_to_file();

    // Species
    Simulation.add_specie(311.39,1.30,26.7,prod_cp,prod_k,prod_mu);          //Products
    Simulation.add_specie(188,1.31,44,oxi_cp,oxi_k,oxi_mu);                  //Oxydizer
    Simulation.add_specie(138,1.13,60,fuel_cp,fuel_k,fuel_mu);               //Fuel

    // set the fuel properties
    thermo::species[2].h_vap = 300e3;
    thermo::species[2].T_ref = 350;
    thermo::species[2].p_ref = 101325;
    thermo::species[2].rho_liq = 700;
    thermo::species[2].C = 2680;

    // Initial conditions
    double rho = p2/(thermo::r_mix_comp(init_comp)*T0);
    double u = m_OX/(rho*Simulation.msh.A[0]);
    double momentum = rho*u;
    double energy = rho*thermo::enthalpy(T0,init_comp) -p2 + rho*u*u/2;

    // Initial conditions
    Simulation.initial_conditions(std::vector<double>{rho,rho,0,momentum,energy});

    // Lagrangian particles manager and specie init
    Simulation.init_particles(1e6,1e3,1000);
    Simulation.add_lagrangian_mono_particles(2,m_F,700,30e-6 ,-0.0004,0.5*u,300,300,1e5,1e3);
    // Simulation.add_lagrangian_unif_particles(2,20,700,1e-4,1e-6,0,15,300,300,1e5,1e-3);
    // Simulation.add_lagrangian_norm_particles(2,m_F,700,30e-6,1e-6,0,40,300,300,1e5,1e-3);

    // Boundary functions and values
    Simulation.add_boundary_function(boundary::mass_flow_inlet,std::vector<double>{m_OX,300,0,1,0});
    Simulation.add_boundary_function(boundary::supersonic_outlet,std::vector<double>{p2});

    // Solving
    Simulation.solve(2,1000,0.2);

    // Exporting all lagrangian particles
    Simulation.var.export_to_file(Simulation.msh,Simulation.par_man.particles);
    Simulation.export_particles(Simulation.par_man.particles);

    return 0;
}
