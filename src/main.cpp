#include <iostream>
#include <vector>

#include "reflow.hpp"
#include "boundary_cond.hpp"

double kappa = 1.23;
double r = 451;

// double kappa = 1.4;
// double r = 287;

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
    reflow S(500,0,0.25,4,std::vector<double>{18.47,1,267,10.87e6});
    S.spline_geometry(curves,100);
    S.apply_heat_source(7e6,0,0.25);
    // S.init_particles(200000,100000,100);
    S.set_boundary(boundary::subsonic_inlet,std::vector<double>{1,600,r,kappa,1,0.8},boundary::subsonic_outlet,std::vector<double>{100000,kappa});
    S.solve();

    // test kapiček v trubce
    // reflow S(500,3,std::vector<double>{1.16,0,250000});
    // S.init_particles(200000,1000,100);
    // S.set_boundary(boundary::subsonic_inlet,std::vector<double>{100,300,r,kappa},boundary::zero_gradient_r,std::vector<double>{1e5,kappa});
    // S.solve();

    return 0;
}
