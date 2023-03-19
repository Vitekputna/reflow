#include <iostream>
#include <vector>
#include <memory>

#include "reflow.hpp"
#include "boundary_cond.hpp"
#include "initial_cond.hpp"

#include "geometry.hpp"

std::vector<double> prod_cp = {0,  9.95096576e+02,  4.87796972e-01, -1.33106793e-04, 1.32096919e-08,  3.86499580e-13, -9.79511577e-17};
std::vector<double> fuel_cp = {0, 2.95260582e+02, 4.95204053e+00, -2.60871236e-03, 7.03837556e-07, -9.39772178e-11, 4.90497393e-15};
std::vector<double> oxi_cp = {0, 6.27400878e+02, 1.09090162e+00, -6.21904849e-04, 1.77914259e-07, -2.46557076e-11, 1.31958533e-15};

const auto init_comp = std::vector<double>{1,0,0};
double p_0 = 101325;
double T_0 = 300;

int main(int argc, char** argv)
{
    // // geometrie
    std::vector<std::vector<std::vector<double>>> curves;
    std::vector<std::vector<double>> curve = {{0,4.418e-3},{0.15,4.418e-3},{1e-1,0},{1e-2,0}};
    curves.push_back(curve);
    curve = {{0.15,4.418e-3},{0.2,8.553e-4},{0.5e-1,0},{0.5e-1,0}};
    curves.push_back(curve);
    // curve = {{0.2,8.553e-4},{0.319,4.185e-3},{0.5e-1,0},{2.1e-1,0}};
    curve = {{0.2,8.553e-4},{0.319,3e-3},{0.5e-1,0},{2.1e-1,0}};
    curves.push_back(curve);

    // std::vector<std::unique_ptr<geometry::curve>> geo;
    // geometry_vector geo;
    // geo.emplace_back(new geometry::line(0.0,0.0,2.0,1.0));
    // geo.emplace_back(new geometry::line(2.0,1.0,3.0,0.0));
    // geo.emplace_back(new geometry::arc(1.0,2.0,1.0,1.0,2.0,1.0));
    // geometry::test(geo);

    // výpočet motoru
    reflow S;
    S.refine_mesh(std::vector<std::vector<double>>{{0,0.319,500}});
    S.spline_geometry(curves,100);

    // Species
    S.add_specie(311.39,1.30,26.7,prod_cp); //Products
    S.add_specie(188,1.31,44,oxi_cp); //Oxydizer
    S.add_specie(138,1.13,60,fuel_cp); //Fuel

    // // Chemistry
    // // reaction R(std::vector<int>{1,2},std::vector<int>{},std::vector<double>{0.8684,0.1316},std::vector<double>{}, 5.719e6);
    // // S.add_reaction(R);

    // S.initial_conditions(init::flow(5,p_0,T_0,0,init_comp));
    S.initial_conditions(1,0,init::flow(6,p_0,T_0,0,init_comp));

    double md = 1.1943;
    double OF = 6.6;

    double m_F = md/(OF+1);
    double m_OX = md-m_F;

    std::cout << "Fuel: " << m_F << ", Oxydizer: " << m_OX << "\n";

    S.apply_mass_source(m_F,300,0.005,0.08,std::vector<double>{0,0,1});

    S.set_boundary(boundary::subsonic_inlet,std::vector<double>{m_OX,300,0,1,0}
                  ,boundary::subsonic_outlet,std::vector<double>{p_0});

    S.solve();



    // reflow S;
    // S.refine_mesh(std::vector<std::vector<double>>{{0,10,10000}});

    // S.add_specie(311.39,1.30,26.7,prod_cp);     //Products
    // S.add_specie(188,1.31,44,oxi_cp);           //Oxydizer
    // S.add_specie(138,1.13,60,fuel_cp);          //Fuel

    // // S.initial_conditions(init::flow(5,p_0,T_0,9.219,init_comp));
    // S.initial_conditions(1,1,init::flow(6,p_0,T_0,9.219,init_comp));

    // S.particle_distribution();

    // S.var.export_to_file(S.msh);

    // S.set_boundary(boundary::subsonic_inlet,std::vector<double>{10,300,1,0,0}
    //               ,boundary::subsonic_outlet,std::vector<double>{p_0});

    // S.solve();

    return 0;
}
