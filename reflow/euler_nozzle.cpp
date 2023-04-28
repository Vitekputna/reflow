#include <iostream>
#include <vector>
#include <memory>

#include "../src/reflow.hpp"
#include "../src/boundary_cond.hpp"
#include "../src/initial_cond.hpp"
#include "../src/params.hpp"
#include "../src/evaporation.hpp"

const auto init_comp = std::vector<double>{0,1,0};
double p0 = 25e5;
double T0 = 3000;
double p2 = 101325;
double md = 1.34;
double OF = 6.6;

int N_frac = 5;

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
    S.refine_mesh(std::vector<std::vector<double>>{{0,0.02,200},{0.02,0.319,500}});
    S.spline_geometry(curves,100);

    // S.load_old_data("out/",3,2*N_frac,0);

    S.msh.export_to_file();

    // Species
    S.add_specie(311.39,1.30,26.7,prod_cp,prod_k,prod_mu);          //Products
    S.add_specie(188,1.31,44,oxi_cp,oxi_k,oxi_mu);                  //Oxydizer
    S.add_specie(138,1.13,60,fuel_cp,fuel_k,fuel_mu);               //Fuel

    // set the fuel properties
    thermo::species[2].h_vap = 666e3;
    thermo::species[2].T_ref = 350;
    thermo::species[2].p_ref = 101325;

    // std::cout << evaporation::fuel_mass_fraction(101325,400) << "\n";

    S.initial_conditions(2*N_frac,0,init::nozzle(S.msh.N,2*N_frac+5,md,300,p0,p2,0.15,init_comp,S.msh));

    // S.apply_heat_source(1e6,0.01,0.05);

    std::cout << "Fuel: " << m_F << ", Oxydizer: " << m_OX << "\n";

    S.add_boundary_function(boundary::mass_flow_inlet_with_droplets,boundary::flow_with_droplets(m_OX,3200,init_comp,N_frac,m_F,700,0.5e-4,3e-6));

    S.add_boundary_function(boundary::supersonic_outlet,std::vector<double>{p2});

    S.solve(0.02,1000,0.3);

    return 0;
}
