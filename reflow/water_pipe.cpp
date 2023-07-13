#include <iostream>
#include <vector>
#include <memory>

#include "../src/reflow.hpp"
#include "../src/boundary_cond.hpp"
#include "../src/initial_cond.hpp"
#include "../src/params.hpp"
#include "../src/euler_droplets.hpp"

const auto init_comp = std::vector<double>{0.99,0,0.01};
double p0 = 101325;
double T0 = 300;
double p2 = 101325;

int N_frac = 1;
int N_comp = 3;
bool droplet_momentum = true;
bool droplet_energy = true;

int N_var = N_comp+2*N_frac+droplet_momentum*N_frac+droplet_energy*N_frac+2;

double m_A = 100;
double m_W = 0.01;

int main(int argc, char** argv)
{
    // výpočet motoru
    reflow S;
    S.refine_mesh(std::vector<std::vector<double>>{{0,5,2000}});

    S.msh.export_to_file();

    // Species
    S.add_specie(287,1.3,29,air_cp,air_k,air_mu);                //Air
    S.add_specie(138,1.13,60,fuel_cp,fuel_k,fuel_mu);            //Dummy
    S.add_specie(461,1.4,18,steam_cp,steam_k,steam_mu);              //Water
    

    // set the fuel properties
    thermo::species[2].h_vap = 2257e3;
    thermo::species[2].T_ref = 373.15;
    thermo::species[2].p_ref = 101325;
    thermo::species[2].rho_liq = 1000;
    thermo::species[2].C = 4182;

    // Initial conditions
    double rho = p0/(thermo::r_mix_comp(init_comp)*T0);
    double u = m_A/(rho*S.msh.A[0]);
    double momentum = rho*u;
    double energy = rho*thermo::enthalpy(T0,init_comp) -p0 + rho*u*u/2;

    // Initial conditions
    S.initial_conditions(N_frac,true,true,std::vector<double>{rho,0,0,0,0,0,0,momentum,energy});

    using namespace boundary;

    S.add_boundary_function(mass_flow_inlet,std::vector<double>{m_A,300,0.99,0,0.01});
    S.add_boundary_function(active_thermal_drop_inlet,active_thermal_droplets(normal_distribution,N_frac,m_W,1000,300,0.5*u,50e-6,1e-6));
    S.add_boundary_function(subsonic_outlet,std::vector<double>{p2});

    S.solve(1,1000,0.1);
    S.var.export_to_file(S.msh,S.par_man.particles);

    return 0;
}
