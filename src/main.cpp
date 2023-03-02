#include <iostream>
#include <vector>
#include <memory>

#include "reflow.hpp"
#include "boundary_cond.hpp"
#include "initial_cond.hpp"

#include "geometry.hpp"

std::vector<double> prod_cp = {92385.20332747149, 1307.9320148736522, 1.1843020600858043, -0.0004719070244224573, 1.018954129756949e-07, -1.1264527612965708e-11, 5.023488301548345e-16};
std::vector<double> fuel_cp = {-235519.57542488514, 1455.374204197742, 3.091173168072451, -0.001445727217091643, 3.428476619622498e-07, -4.0705333490788554e-11, 1.9221666948167e-15};
std::vector<double> oxi_cp = {-8000.032796543984, 629.1342945579523, 1.1223631596599488, -0.0006677168150395522, 1.948637060185014e-07, -2.7403372013652136e-11, 1.4815269244552882e-15};

// old thermo time 23s
// new thermo time 12s :)


int main(int argc, char** argv)
{
    // geometrie
    std::vector<std::vector<std::vector<double>>> curves;
    std::vector<std::vector<double>> curve = {{0,4.418e-3},{0.15,4.418e-3},{1e-1,0},{1e-2,0}};
    curves.push_back(curve);
    curve = {{0.15,4.418e-3},{0.2,8.553e-4},{0.5e-1,0},{0.5e-1,0}};
    curves.push_back(curve);
    curve = {{0.2,8.553e-4},{0.319,4.185e-3},{0.5e-1,0},{2.1e-1,0}};
    curves.push_back(curve);

    // std::vector<std::unique_ptr<geometry::curve>> geo;
    // geometry_vector geo;
    // geo.emplace_back(new geometry::line(0.0,0.0,2.0,1.0));
    // geo.emplace_back(new geometry::line(2.0,1.0,3.0,0.0));
    // geo.emplace_back(new geometry::arc(1.0,2.0,1.0,1.0,2.0,1.0));
    // geometry::test(geo);

    // výpočet motoru
    reflow S;
    S.refine_mesh(std::vector<std::vector<double>>{{0,0.319,1000}});
    S.spline_geometry(curves,100);

    S.msh.export_to_file();

    // Species
    S.add_specie(291,1.23,31,prod_cp); //Products
    S.add_specie(188,1.31,44,oxi_cp); //Oxydizer
    S.add_specie(138,1.13,60,fuel_cp); //Fuel

    // Chemistry
    // reaction R(std::vector<int>{1,2},std::vector<int>{},std::vector<double>{0.8684,0.1316},std::vector<double>{}, 5.719e6);
    // S.add_reaction(R);

    S.initial_conditions(init::flow(5,101325,300,0,std::vector<double>{1,0,0}));

    // S.init_particles(200000,10000,100);

    double md = 1.1943;
    double OF = 6.6;

    double m_F = md/(OF+1);
    double m_OX = md-m_F;

    std::cout << "Fuel: " << m_F << ", Oxydizer: " << m_OX << "\n";

    // S.apply_heat_source(6.837e6,0.005,0.08);
    // S.apply_mass_source(0.1573,300,0.005,0.08,std::vector<double>{0,0,1});

    S.apply_mass_source(m_F,300,0.005,0.08,std::vector<double>{0,0,1}); 

    S.set_boundary(boundary::subsonic_inlet,std::vector<double>{m_OX,300,0,1,0}
                  ,boundary::subsonic_outlet,std::vector<double>{101325});

    // S.set_boundary(boundary::subsonic_inlet,std::vector<double>{md,3225,1,0,0}
    //               ,boundary::zero_gradient_r,std::vector<double>{101325});

    S.solve();

    return 0;
}
