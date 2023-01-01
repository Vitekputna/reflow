#include "particle.hpp"
#include <vector>
#include "thermodynamics.hpp"
#include <iostream>

particle::particle(){}

particle::particle(unsigned int _N, double _x, double _r, double _u, double _rho, double _T) : N{_N}, r{_r}, u{_u}, rho{_rho}, T{_T}, x{_x}
{
    x_0 = x;
    r_0 = r;
    u_0 = u;
    T_0 = T;
    rho_0 = rho;

    // std::cout << N << " " << r << "\n";
    m = 4/3*3.14159*r*r*r*rho; // hmotnost jedne kapičky
    M = N*m; // hmotnost grupy
}

void particle::reset()
{
    x = x_0;
    r = r_0;
    u = u_0;
    T = T_0;
    rho = rho_0;
    last_cell_idx = 0;
    m = 4/3*3.14159*r*r*r*rho; // hmotnost jedne kapičky
    M = N*m; // hmotnost grupy
}