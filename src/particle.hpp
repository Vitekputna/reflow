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
    int N_max, max_per_group, N; // max. number of particles, max num, of real particles per group, current number of usable paricles
    int n_add = 1000; // number of particles to add if not enough

    public:
    std::vector<particle> particles; // holds particle objects

    particle_manager();
    particle_manager(int _N_max, int _max_per_group, int _N);

    bool particle_inlet(double m, double r, double u, double x, double rho, double T);
    bool particle_inlet(double m, double r_from, double r_to, double u_from, double u_to, double x_from, double x_to, double rho, double T);

    bool spawn_particles(int n_groups, int n_particles, double r, double u, double x, double rho, double T);
    bool spawn_particles(int n_groups, int n_particles, double r_from, double r_to, double u_from, double u_to, double x_from, double x_to, double rho, double T);
};