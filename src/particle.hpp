#pragma once
#include <vector>

struct particle
{
    double r, u, T, rho; // vlastnosti grupy kapiček
    double m, M; // hmotnost jedné kapičky, hmotnost celé grupy
    double x = 0; // poloha grupy
    unsigned int N ; // počet kapiček v grupě

    double x_0, r_0, u_0, T_0, rho_0; //počáteční hodnoty grupy

    int last_cell_idx = 0;

    particle();
    particle(unsigned int _N, double _x, double _r, double _u, double _rho, double _T);
    void reset();
};