#pragma once
#include "particle.hpp"
#include "variables.hpp"
#include "mesh.hpp"
#include <vector>

namespace lagrange_solver
{
    double integrate_particle(double dt, particle& P, std::vector<double>& W, std::vector<double>& res);
    void update_particles(double dt, std::vector<particle>& particles, variables& var, mesh const& msh, std::vector<std::vector<double>>& res);

    inline double acceleration(double Cd, double r, double du);
    inline double heat_flux(double alfa, double C, double dT);
    inline double radius_change(double D, double r, double rho, double dT);
};
