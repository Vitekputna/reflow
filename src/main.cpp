#include <iostream>
#include <vector>
#include "reflow.hpp"
#include "initial_cond.hpp"
#include "boundary_cond.hpp"

#include "mesh.hpp"
#include "lagrange_solver.hpp"
#include "particle.hpp"

#include <fstream>

// double kappa = 1.23;
// double r = 451;

double kappa = 1.4;
double r = 287;

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
    reflow S(500,0,0.25,3,std::vector<double>{18.47,267,10.87e6});
    S.spline_geometry(curves,100);
    S.apply_heat_source(7e6,0.03,0.05);
    S.set_boundary(boundary::subsonic_inlet,std::vector<double>{1.18,320,r,kappa},boundary::zero_gradient_r,std::vector<double>{1e5,kappa});
    S.solve();

    // test kapiček v trubce
    // reflow S(500,3,std::vector<double>{1.16,0,250000});
    // S.set_boundary(boundary::subsonic_inlet,std::vector<double>{100,300,r,kappa},boundary::zero_gradient_r,std::vector<double>{1e5,kappa});
    // S.solve();

    // std::vector<std::vector<double>> res{{}};

    // particle P(1,0,1e-6,70,1000,300);
    // std::vector<particle> particles;
    // particles.push_back(P);

    // auto stream = std::ofstream("up.txt");
    // double dt, t = 0;

    // while(particles[0].x < 0.25)
    // {
    //     stream << t << " " << particles[0].x << " " << particles[0].u << " " << particles[0].T << "\n";
    //     dt = lagrange_solver::update_particles(1e-5,particles,S.var,S.msh,res);
    //     t += dt;    
    // }

    // stream.close();

    for(int i = 0; i < 100; i++)
    {
        std::cout << S.particles[i].r << "\n";
    }

    return 0;
}
