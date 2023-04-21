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

    m = 4*M_PI*pow(r,3)*rho/3; // hmotnost jedne kapičky
    M = N*m; // hmotnost grupy
}

void particle::reset()
{
    in_use = false;
}

particle_manager::particle_manager(){}

particle_manager::particle_manager(int _N_max, int _max_per_group, int _N) : N_max{_N_max}, max_per_group{_max_per_group}, N{_N}
{
    particles.resize(N);
}

void particle_manager::apply_boundary(double dt)
{
    int i = 0;
    for(auto ptr : boundary_functions)
    {
        (this->*ptr)(dt,boundary_parameters[i]);
        i++;
    }
}

void particle_manager::add_monodispersion(std::vector<double> parameters)
{
    boundary_functions.push_back(&particle_manager::monodispersion_particle_inlet);
    boundary_parameters.push_back(parameters);
}

void particle_manager::add_uniform(std::vector<double> parameters)
{
    boundary_functions.push_back(&particle_manager::uniform_particle_inlet);
    boundary_parameters.push_back(parameters);
}

void particle_manager::add_normal(std::vector<double> parameters)
{
    boundary_functions.push_back(&particle_manager::normal_particle_inlet);
    boundary_parameters.push_back(parameters);
}

// mass, r, u, x, rho, T
void particle_manager::monodispersion_particle_inlet(double dt, std::vector<double>& parameters)
{
    const double m = parameters[0];
    const double r = parameters[1];
    const double rho = parameters[4];

    int N = dt*m/(4*M_PI*pow(r,3)*rho/3);

    spawn_particles_monodispersion(1,N,parameters);
}

// mass, r_mean, r_var, u, x, rho, T
void particle_manager::uniform_particle_inlet(double dt, std::vector<double>& parameters)
{
    const double m = parameters[0];
    const double r_mean = parameters[1];
    const double rho = parameters[5];

    int N = dt*m/(4*M_PI*pow(r_mean,3)*rho/3);

    spawn_particles_uniform(1,N,parameters);
}

void particle_manager::normal_particle_inlet(double dt, std::vector<double>& parameters)
{
    const double m = parameters[0];
    const double r_mean = parameters[1];
    const double rho = parameters[5];

    int N = dt*m/(4*M_PI*pow(r_mean,3)*rho/3);

    spawn_particles_normal(1,N,parameters);
}

bool particle_manager::spawn_particles_monodispersion(int n_groups, int n_particles, std::vector<double>& parameters)
{
    const double r = parameters[1];
    const double u = parameters[2];
    const double x = parameters[3];
    const double rho = parameters[4];
    const double T = parameters[5];

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

bool particle_manager::spawn_particles_uniform(int n_groups, int n_particles, std::vector<double>& parameters)
{
    const double r_mean = parameters[1];
    const double r_var = parameters[2];
    const double u = parameters[3];
    const double x = parameters[4];
    const double rho = parameters[5];
    const double T = parameters[6];

    //randomizer for radii
    std::random_device rnd;
    std::mt19937 gen(rnd());
    std::uniform_real_distribution<double> r_dist(r_mean-r_var,r_mean+r_var);

    int n = 0;
    int n_par_group = n_particles/n_groups;

    double dr;

    for(auto& par : particles)
    {
        if(!par.in_use)
        {
            dr = r_dist(gen);
            par = particle(n_par_group,x,dr,u,rho,T);
            n++;
        }

        if(n == n_groups) return true;
    }

    if(N+n_add > N_max) throw std::overflow_error("Number of particles is too big");
    particles.resize(N+n_add);
    N+=n_add;
    return false;
}

bool particle_manager::spawn_particles_normal(int n_groups, int n_particles, std::vector<double>& parameters)
{
    const double r_mean = parameters[1];
    const double r_var = parameters[2];
    const double u = parameters[3];
    const double x = parameters[4];
    const double rho = parameters[5];
    const double T = parameters[6];

    //randomizer for radii
    std::random_device rnd;
    std::mt19937 gen(rnd());
    std::normal_distribution<double> r_dist(r_mean,r_var);

    int n = 0;
    int n_par_group = n_particles/n_groups;

    double dr;

    for(auto& par : particles)
    {
        if(!par.in_use)
        {
            dr = r_dist(gen);
            par = particle(n_par_group,x,dr,u,rho,T);
            n++;
        }

        if(n == n_groups) return true;
    }

    if(N+n_add > N_max) throw std::overflow_error("Number of particles is too big");
    particles.resize(N+n_add);
    N+=n_add;
    return false;
}