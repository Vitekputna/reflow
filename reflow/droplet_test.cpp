#include <iostream>
#include <vector>
#include <memory>

#include "../src/reflow.hpp"
#include "../src/boundary_cond.hpp"
#include "../src/initial_cond.hpp"

std::vector<double> prod_cp = {9.95096576e+02,  4.87796972e-01, -1.33106793e-04, 1.32096919e-08,  3.86499580e-13, -9.79511577e-17};
std::vector<double> fuel_cp = {2.95260582e+02, 4.95204053e+00, -2.60871236e-03, 7.03837556e-07, -9.39772178e-11, 4.90497393e-15};
std::vector<double> oxi_cp = {6.27400878e+02, 1.09090162e+00, -6.21904849e-04, 1.77914259e-07, -2.46557076e-11, 1.31958533e-15};

const auto init_comp = std::vector<double>{0,1,0};
double p0 = 25e5;
double T0 = 3000;
double p2 = 101325;
double md = 1.34;
double OF = 6.6;

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
    S.refine_mesh(std::vector<std::vector<double>>{{0,0.02,200},{0.02,0.319,500}});
    S.spline_geometry(curves,100);

    S.msh.export_to_file();

    // Species
    S.add_specie(311.39,1.30,26.7,prod_cp);     //Products
    S.add_specie(188,1.31,44,oxi_cp);           //Oxydizer
    S.add_specie(138,1.13,60,fuel_cp);          //Fuel

    S.initial_conditions(10,0,init::nozzle(S.msh.N,15,md,300,p0,p2,0.15,init_comp,S.msh));

    std::cout << "Fuel: " << m_F << ", Oxydizer: " << m_OX << "\n";

    // S.add_boundary_function(boundary::mass_flow_inlet_with_droplets,std::vector<double>{m_OX,300,0,1,0,3,m_F/3,0.8e-3,m_F/3,1e-3,m_F/3,1.2e-3,700});

    S.add_boundary_function(boundary::mass_flow_inlet_with_droplets,std::vector<double>{m_OX,300,0,1,0,5,0.00021123295588010936,0.0009620007732977866,
                                                                                                            0.02624689043447874,0.0009825778335494264,
                                                                                                            0.12012475625169018,0.0010001171423589421,
                                                                                                            0.029151618922297836,0.0010175617214114906,
                                                                                                            0.0002653996573385169,0.0010380553687645798,700});

    S.add_boundary_function(boundary::supersonic_outlet,std::vector<double>{p2});

    thermo::update(S.var.W);
    S.var.export_to_file(S.msh);

    S.solve();

    return 0;
}
