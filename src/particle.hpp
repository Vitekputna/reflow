#pragma once
#include <vector>

struct particle
{
    double r, u, T, rho; // vlastnosti grupy kapiček
    double m, M; // hmotnost jedné kapičky, hmotnost celé grupy
    double x = 0; // poloha grupy
    unsigned int N ; // počet kapiček v grupě

    unsigned int specie_idx = 0;

    double x_0, r_0, u_0, T_0, rho_0; //počáteční hodnoty grupy

    int last_cell_idx = 0;
    bool in_use = false;

    particle();
    particle(unsigned int _N, double _x, double _r, double _u, double _rho, double _T);
    void reset();
};



class particle_manager
{
    public:
    int N_max, max_per_group, N = 0; // max. number of particles, max num, of real particles per group, current number of usable paricles
    int n_add = 1000; // number of particles to add if not enough

    public:
    std::vector<particle> particles; // holds particle objects

    // boundary/inlets
    using lagrange_func = void (particle_manager::*)(double, std::vector<double>&);

    std::vector<lagrange_func> boundary_functions;
    std::vector<std::vector<double>> boundary_parameters;

    particle_manager();
    particle_manager(int _N_max, int _max_per_group, int _N);
    
    void apply_boundary(double dt);

    void add_monodispersion(std::vector<double> parameters);
    void add_uniform(std::vector<double> parameters);
    void add_normal(std::vector<double> parameters);

    void monodispersion_particle_inlet(double dt, std::vector<double>& parameters);
    void uniform_particle_inlet(double dt, std::vector<double>& parameters);
    void normal_particle_inlet(double dt, std::vector<double>& parameters);

    bool spawn_particles_monodispersion(int n_groups, int n_particles, std::vector<double>& parameters);
    bool spawn_particles_uniform(int n_groups, int n_particles, std::vector<double>& parameters);
    bool spawn_particles_normal(int n_groups, int n_particles, std::vector<double>& parameters);
};

