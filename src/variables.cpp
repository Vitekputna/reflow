#include "variables.hpp"
#include "thermodynamics.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

extern double kappa,r;

variables::variables(){}

variables::variables(int _N_var, int _N) : N_var{_N_var}, N{_N+2}, N_walls{_N+1}
{
    W = std::vector<std::vector<double>>(N,std::vector<double>(N_var,0.0));
    flux = std::vector<std::vector<double>>(N_walls,std::vector<double>(N_var,0.0));
    exact_flux = std::vector<std::vector<double>>(N,std::vector<double>(N_var,0.0));

    q = std::vector<double>(N,0.0);
    md = std::vector<double>(N,0.0);
}

variables::variables(int _N_var, int _N, std::vector<double> const& W_0) : variables(_N_var,_N)
{
    for(auto& wi : W)
    {
        for(int k = 0; k < N_var; k++)
        {
            wi[k] = W_0[k];
        }
    }
}

variables::variables(int _N_var, int _N, std::vector<std::vector<double>> const& W_0) : variables(_N_var, _N)
{
    for(int i = 0; i < N; i++)
    {
        for(int k = 0; k < N_var; k++)
        {
            W[i][k] = W_0[i][k];
        }
    }
}

void variables::apply_heat_source(double Q_tot, double x_from, double x_to, mesh const& msh)
{
    double V = 0;

    for(int i = 1; i < msh.N-1; i++)
    {
        if(msh.x[i] > x_from && msh.x[i] < x_to)
        {
            V += (msh.xf[i] - msh.xf[i-1])*msh.A[i];
        }
    }

    double Q_V = Q_tot/V;

    for(int i = 0; i < N; i++)
    {
        if(msh.x[i] > x_from && msh.x[i] < x_to)
        {
            q[i] = Q_V*msh.A[i];
        }
    }

    std::cout << "Source: " << Q_V << " " << V << "\n";
}

void variables::export_to_file(std::string path, mesh const& msh)
{
    auto stream = std::ofstream();

    for(int k = 0; k < N_var; k++)
    {
        stream =  std::ofstream("out/W" + std::to_string(k) + ".txt");

        for(int i = 0; i < N; i++)
        {
            stream << msh.x[i] << " " << W[i][k] << "\n";
        }

        stream << "\n";
        stream.close();
    }

    stream =  std::ofstream("out/p.txt");

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << thermo::pressure(W[i],kappa) << "\n";
    }

    stream << "\n";
    stream.close();

    stream =  std::ofstream("out/u.txt");

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << W[i][1]/W[i][0] << "\n";
    }

    stream << "\n";
    stream.close();

    stream =  std::ofstream("out/T.txt");

    double T_max = 0;

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << thermo::temperature(W[i],kappa,r) << "\n";
        T_max = std::max(T_max,thermo::temperature(W[i],kappa,r));
    }

    // std::cout << T_max << "\n";

    stream << "\n";
    stream.close();

    stream =  std::ofstream("out/a.txt");

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << thermo::speed_of_sound(W[i],kappa) << "\n";
    }

    stream << "\n";
    stream.close();
}

void variables::export_timestep(double t, mesh const& msh, std::vector<particle> const& particles)
{
    auto stream = std::ofstream("out/timesteps/data_" + std::to_string(t) + ".txt");
    
    stream << "x[m]\trho[kg/m3]\tu[m/s]\tp[Pa]\tT[K]\tP_r[m]\tP_u[m/s]\tP_T[K]\tP_N[-]\n";

    std::vector<std::vector<double>> particle_values(N,std::vector<double>(4,0.0));

    for(auto const& P : particles)
    {
        if(P.in_use) particle_values[P.last_cell_idx][3] += P.N;
    }

    int n;
    for(auto const& P : particles)
    {
        if(P.in_use)
        {
            n = particle_values[P.last_cell_idx][3];
            particle_values[P.last_cell_idx][0] += P.N*P.r/n;
            particle_values[P.last_cell_idx][1] += P.N*P.u/n;
            particle_values[P.last_cell_idx][2] += P.N*P.T/n;
        }
    }

    for(int i = 1; i < N-1; i++)
    {
        stream << msh.x[i] << "\t" << W[i][0] << "\t" << W[i][1]/W[i][0] << "\t" << thermo::pressure(W[i],kappa) << "\t" << thermo::temperature(W[i],kappa,r)
                           << "\t" << particle_values[i][0] << "\t" << particle_values[i][1] << "\t" << particle_values[i][2] << "\t" << particle_values[i][3] << "\n";
    }

}