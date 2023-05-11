#include <iostream>
#include <vector>
#include <memory>

#include "../src/reflow.hpp"
#include "../src/boundary_cond.hpp"
#include "../src/initial_cond.hpp"
#include "../src/params.hpp"
#include "../src/euler_droplets.hpp"

const auto init_comp = std::vector<double>{0,1,0};
double p0 = 101325;
double T0 = 300;
double p2 = 101325;
double md = 100;
double OF = 6.6;

int N_frac = 5;
int N_comp = 3;
bool droplet_momentum = true;
bool droplet_energy = true;

int N_var = N_comp+2*N_frac+droplet_momentum*N_frac+droplet_energy*N_frac+2;

double m_F = md/(OF+1);
double m_OX = md-m_F;

int main(int argc, char** argv)
{
    // výpočet motoru
    reflow S;
    S.set_numThreads(2);
    S.refine_mesh(std::vector<std::vector<double>>{{0,5,1000}});
    S.msh.constant_area(0.002);

    // S.load_old_data("out/",N_comp,N_frac,true,true);

    S.msh.export_to_file();

    // Species
    S.add_specie(311.39,1.30,26.7,prod_cp,prod_k,prod_mu);          //Products
    S.add_specie(188,1.31,44,oxi_cp,oxi_k,oxi_mu);                  //Oxydizer
    S.add_specie(138,1.13,60,fuel_cp,fuel_k,fuel_mu);               //Fuel

    // set the fuel properties
    thermo::species[2].h_vap = 300e3;
    thermo::species[2].T_ref = 350;
    thermo::species[2].p_ref = 101325;
    thermo::species[2].rho_liq = 700;
    thermo::species[2].C = 2680;

    // Initial conditions
    double rho = p2/(thermo::r_mix_comp(init_comp)*T0);
    double u = m_OX/(rho*S.msh.A[0]);
    double momentum = rho*u;
    double energy = rho*thermo::enthalpy(T0,init_comp) -p2 + rho*u*u/2;

    // Initial conditions
    S.initial_conditions(N_frac,true,true,std::vector<double>{rho,rho,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,momentum,energy});

    std::cout << "Fuel: " << m_F << ", Oxydizer: " << m_OX << "\n";

    using namespace boundary;

    S.add_boundary_function(mass_flow_inlet,std::vector<double>{m_OX,300,0,1,0}); 
    // S.add_boundary_function(active_thermal_drop_inlet,active_thermal_droplets(normal_distribution,N_frac,m_F,700,300,40,30e-6,1e-6));
    S.add_boundary_function(subsonic_outlet,std::vector<double>{p2});

    S.solve(0.2,1000,0.2);
    S.var.export_to_file(S.msh,S.par_man.particles);

    return 0;
}
