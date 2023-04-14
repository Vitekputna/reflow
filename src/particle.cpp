#include "particle.hpp"
#include <vector>
#include "thermodynamics.hpp"
#include <iostream>
#include <random>
#include <exception>

particle::particle(){}

particle::particle(unsigned int _N, double _x, double _r, double _u, double _rho, double _T) : N{_N}, r{_r}, u{_u}, rho{_rho}, T{_T}, x{_x}
{
    in_use = true;

    x_0 = x;
    r_0 = r;
    u_0 = u;
    T_0 = T;
    rho_0 = rho;

    last_cell_idx = 0;

    // std::cout << N << " " << r << "\n";
    m = 4/3*3.14159*r*r*r*rho; // hmotnost jedne kapičky
    M = N*m; // hmotnost grupy
}

void particle::reset()
{
    in_use = false;

    x = x_0;
    r = r_0;
    u = u_0;
    T = T_0;
    rho = rho_0;
    last_cell_idx = 0;
    m = 4/3*3.14159*r*r*r*rho; // hmotnost jedne kapičky
    M = N*m; // hmotnost grupy
}

particle_manager::particle_manager(){}

particle_manager::particle_manager(int _N_max, int _max_per_group, int _N) : N_max{_N_max}, max_per_group{_max_per_group}, N{_N}
{
    particles.resize(N);
}

bool particle_manager::particle_inlet(double m, double r, double u, double x, double rho, double T)
{
    int N = m/(4*3.14159*r*r*r*rho/3);
    spawn_particles(1,N,r,u,x,rho,T);
    return true;
}

bool particle_manager::particle_inlet(double m, double r_from, double r_to, double u_from, double u_to, double x_from, double x_to, double rho, double T)
{
    double r = 0.5*(r_from+r_to);
    int N = m/(4*3.14159*r*r*r*rho/3);
    spawn_particles(1,N,r_from,r_to,u_from,u_to,x_from,x_to,rho,T);
    return true;
}

bool particle_manager::spawn_particles(int n_groups, int n_particles, double r, double u, double x, double rho, double T)
{
    int n = 0;
    int n_par_group = n_particles/n_groups;

    for(auto& par : particles)
    {
        if(!par.in_use)
        {
            par = particle(n_par_group,x,r,u,rho,T);
            n++;
        }

        if(n == n_groups) return true;
    }

    if(N+n_add > N_max) throw std::overflow_error("Number of particles is too big");
    particles.resize(N+n_add);
    N+=n_add;
    return false;
}

bool particle_manager::spawn_particles(int n_groups, int n_particles,double r_from, double r_to, double u_from, double u_to, double x_from, double x_to, double rho, double T)
{
    //randomizer for velocity
    std::random_device rnd;
    std::mt19937 gen(rnd());
    std::uniform_real_distribution<double> r_dist(r_from,r_to);
    std::uniform_real_distribution<double> u_dist(u_from,u_to);
    std::uniform_real_distribution<double> x_dist(x_from,x_to);

    int n = 0;
    int n_par_group = n_particles/n_groups;

    double dr, du, dx;

    for(auto& par : particles)
    {
        if(!par.in_use)
        {
            dr = r_dist(gen);
            du = u_dist(gen);
            dx = x_dist(gen);
            par = particle(n_par_group,dx,dr,du,rho,T);
            n++;
        }

        if(n == n_groups) return true;
    }

    if(N+n_add > N_max) throw std::overflow_error("Number of particles is too big");
    particles.resize(N+n_add);
    N+=n_add;
    return false;
}