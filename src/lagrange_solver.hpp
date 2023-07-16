#pragma once
#include "particle.hpp"
#include "variables.hpp"
#include "mesh.hpp"
#include <vector>

namespace lagrange_solver
{
    double integrate_particle(double dt, double V, particle& P, std::vector<double>& W, std::vector<double>& res);
    void update_particles(double dt, std::vector<particle>& particles, variables& var, mesh const& msh, std::vector<std::vector<double>>& res);

    inline double acceleration(double r, double rho_l, double rho_g, double mu, double du);
    inline double heat_flux(double r, double T_gas, double T_drop, double rho_g, double u_l, double u_g, double mu, double cp, double k, double h_vap, double phi);
    inline double radius_change(double r, double Ys, double Yf, double rho_l, double rho_g, double u_l, double u_g, double mu, double D);
};
