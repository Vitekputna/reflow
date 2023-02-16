#include <iostream>
#include <vector>

#include "reflow.hpp"
#include "boundary_cond.hpp"
#include "initial_cond.hpp"

std::vector<double> prod_cp = {72913.76001451393, 1119.6547177177858, 1.256662344680612, -0.000555166379245523, 1.3051468187215095e-07, -1.5588756032619083e-11, 7.455841969813714e-16};
std::vector<double> fuel_cp = {-235519.57542488514, 1455.374204197742, 3.091173168072451, -0.001445727217091643, 3.428476619622498e-07, -4.0705333490788554e-11, 1.9221666948167e-15};
std::vector<double> oxi_cp = {-8000.032796543984, 629.1342945579523, 1.1223631596599488, -0.0006677168150395522, 1.948637060185014e-07, -2.7403372013652136e-11, 1.4815269244552882e-15};


int main(int argc, char** argv)
{
    // geometrie
    std::vector<std::vector<std::vector<double>>> curves;
    std::vector<std::vector<double>> curve = {{0,4.418e-3},{0.15,4.418e-3},{1e-1,0},{1e-2,0}};
    curves.push_back(curve);
    curve = {{0.15,4.418e-3},{0.2,8.553e-4},{0.5e-1,0},{0.5e-1,0}};
    curves.push_back(curve);
    curve = {{0.2,8.553e-4},{0.25,3.455e-3},{0.5e-1,0},{1e-1,1e-3}};
    curves.push_back(curve);

    // výpočet motoru
    reflow S;

    // S.refine_mesh(std::vector<std::vector<double>>{{0,0.05,1000},{0.05,0.25,500}});
    S.refine_mesh(std::vector<std::vector<double>>{{0,0.25,500}});
    S.spline_geometry(curves,100);

    S.add_specie(291,1.23,31,prod_cp); //Products
    S.add_specie(188,1.31,44,oxi_cp); //Oxydizer
    S.add_specie(138,1.13,60,fuel_cp); //Fuel

    S.initial_conditions(init::flow(5,1e5,500,0,std::vector<double>{1,0,0}));

    // S.init_particles(200000,10000,100);
    // S.apply_heat_source(10e6,0.005,0.08);

    S.set_boundary(boundary::subsonic_inlet,std::vector<double>{1.18,3200,1,0,0}
                  ,boundary::zero_gradient_r,std::vector<double>{50000});

    S.solve();

    return 0;
}
