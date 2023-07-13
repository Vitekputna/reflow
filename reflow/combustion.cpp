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
double T0 = 3000;
double p2 = 101325;
double md = 1.34;
double OF = 6.6;

int N_frac = 1;
int N_comp = 3;
bool droplet_momentum = true;
bool droplet_energy = true;

int N_var = N_comp+2*N_frac+droplet_momentum*N_frac+droplet_energy*N_frac+2;

double m_F = md/(OF+1);
double m_OX = md-m_F;

int main(int argc, char** argv)
{
    // // geometrie
    std::vector<std::vector<std::vector<double>>> curves;
    std::vector<std::vector<double>> curve = {{0,4.418e-3},{0.15,4.418e-3},{1e-1,0},{1e-2,0}};
    curves.push_back(curve);
    curve = {{0.15,4.418e-3},{0.2,8.553e-4},{0.5e-1,0},{0.5e-1,0}};
    curves.push_back(curve);
    curve = {{0.2,8.553e-4},{0.319,3e-3},{0.5e-1,0},{2.1e-1,0}};
    curves.push_back(curve);

    // výpočet motoru
    reflow S;

    // S.set_export_path("tests/euler500250/");
    S.set_numThreads(6);

    // S.refine_mesh(std::vector<std::vector<double>>{{0,0.1,500},{0.1,0.319,250}});  
    S.refine_mesh(std::vector<std::vector<double>>{{0,0.319,700}});
    S.spline_geometry(curves,100);

    S.export_mesh();

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
    double rho = p0/(thermo::r_mix_comp(init_comp)*T0);
    double u = m_OX/(rho*S.msh.A[0]);
    double momentum = rho*u;
    double energy = rho*thermo::enthalpy(T0,init_comp) -p0 + rho*u*u/2;

    // Initial conditions
    // S.initial_conditions(N_frac,true,true,std::vector<double>{rho,rho,0,0,0,0,0,momentum,energy});
    S.initial_conditions(N_frac,true,true,init::nozzle(S.msh.N,N_var,m_OX,T0,p0,p2,0.15,init_comp,S.msh));

    std::cout << "Fuel: " << m_F << ", Oxydizer: " << m_OX << "\n";

    using namespace boundary;

    S.add_boundary_function(mass_flow_inlet,std::vector<double>{m_OX,300,0,1,0}); 
    S.add_boundary_function(active_thermal_drop_inlet,active_thermal_droplets(normal_distribution,N_frac,m_F,700,300,5,100e-6,1e-6));
    S.add_boundary_function(supersonic_outlet,std::vector<double>{p2});

    S.solve_parallel(1,1000,0.3);
    // S.var.export_to_file(S.msh,S.par_man.particles);

    return 0;
}
