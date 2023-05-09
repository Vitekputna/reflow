#include <iostream>
#include <vector>
#include <memory>

#include "../src/reflow.hpp"
#include "../src/boundary_cond.hpp"
#include "../src/initial_cond.hpp"
#include "../src/params.hpp"
#include "../src/euler_droplets.hpp"

const auto init_comp = std::vector<double>{0,1,0};
double p0 = 25e5;
double T0 = 3000;
double p2 = 101325;
double md = 1.34;
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
    // geometrie
    std::vector<std::vector<std::vector<double>>> curves;
    std::vector<std::vector<double>> curve = {{0,4.418e-3},{0.15,4.418e-3},{1e-1,0},{1e-2,0}};
    curves.push_back(curve);
    curve = {{0.15,4.418e-3},{0.2,8.553e-4},{0.5e-1,0},{0.5e-1,0}};
    curves.push_back(curve);
    curve = {{0.2,8.553e-4},{0.319,3e-3},{0.5e-1,0},{2.1e-1,0}};
    curves.push_back(curve);

    // výpočet motoru
    reflow S;
    S.refine_mesh(std::vector<std::vector<double>>{{0,0.319,500}});
    S.spline_geometry(curves,100);

    // S.load_old_data("out/",3,N_frac,true,true);

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

    S.initial_conditions(N_frac,droplet_momentum,droplet_energy,init::nozzle(S.msh.N,N_var,md,400,p0,p2,0.15,init_comp,S.msh));

    std::cout << "Fuel: " << m_F << ", Oxydizer: " << m_OX << "\n";

    using namespace boundary;

    S.add_boundary_function(mass_flow_inlet,std::vector<double>{m_OX,400,0,1,0});
    S.add_boundary_function(active_thermal_drop_inlet,active_thermal_droplets(normal_distribution,N_frac,m_F,700,300,60,10e-6,1e-6));

    S.add_boundary_function(supersonic_outlet,std::vector<double>{p2});

    thermo::update(S.var.W);
    S.apply_boundary_conditions();

    S.solve(0.5,1000,0.5);
    S.var.export_to_file(S.msh,S.par_man.particles);

    return 0;
}
