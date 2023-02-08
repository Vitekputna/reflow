#include <iostream>
#include <vector>

#include "reflow.hpp"
#include "boundary_cond.hpp"
#include "initial_cond.hpp"

#include "omp.h"

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

    // S.refine_mesh(std::vector<std::vector<double>>{{0,0.05,1000},{0.05,0.25,200}});
    S.refine_mesh(std::vector<std::vector<double>>{{0,0.25,1000}});
    S.spline_geometry(curves,100);

    S.add_specie(291,1.23); //Products
    S.add_specie(188,1.31); //Oxydizer
    S.add_specie(138,1.13); //Fuel

    S.initial_conditions(init::flow(5,1e5,300,0,std::vector<double>{0,1,0}));

    S.init_particles(200000,10000,100);
    // S.apply_heat_source(1e7,0.005,0.05);

    S.set_boundary(boundary::subsonic_inlet,std::vector<double>{1,300,0,1,0}
                  ,boundary::zero_gradient_r,std::vector<double>{150000});

    S.solve();

    return 0;
}
