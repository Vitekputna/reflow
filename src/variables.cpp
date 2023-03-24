#include "variables.hpp"
#include "thermodynamics.hpp"
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <exception>

std::vector<int> variables::drop_mom_idx;

int variables::mom_idx;
int variables::eng_idx;
int variables::N_comp;
int variables::N_drop_frac = 0;

variables::variables(){}

variables::variables(int _N_var, int _N) : N_var{_N_var}, N{_N}, N_walls{_N-1}
{
    variables::N_comp = _N_var-2;
    variables::mom_idx = N_comp;
    variables::eng_idx = N_comp+1;

    W = std::vector<std::vector<double>>(N,std::vector<double>(N_var,0.0));
    flux = std::vector<std::vector<double>>(N_walls,std::vector<double>(N_var,0.0));
    exact_flux = std::vector<std::vector<double>>(N,std::vector<double>(N_var,0.0));
    grad = std::vector<std::vector<double>>(N,std::vector<double>(N_var,0.0));

    q = std::vector<double>(N,0.0);
    md = std::vector<std::vector<double>>(N,std::vector<double>(N_comp,0.0));
}

variables::variables(int _N_var, int _N, std::vector<double> const& W_0) : variables(_N_var,_N)
{
    if((uint)N_var != W_0.size()) throw std::overflow_error("Size of initial condition is not the same as declared number of variables");

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
    if((uint)N_var != W_0.size()) throw std::overflow_error("Size of initial condition is not the same as declared number of variables");

    for(int i = 0; i < N; i++)
    {
        for(int k = 0; k < N_var; k++)
        {
            W[i][k] = W_0[i][k];
        }
    }
}

variables::variables(int _N_var, int _N, int _N_drop_frac, int _N_drop_mom, std::vector<double> const& W_0) : variables(_N_var, _N, W_0)
{
    variables::N_drop_frac = _N_drop_frac;

    variables::drop_mom_idx = std::vector<int>(_N_drop_frac,mom_idx);
    variables::N_comp = (N_var-2) - N_drop_frac;

    std::cout << "##########################################\n";
    std::cout << "Variables:\n";
    std::cout << "Number of compounds:\t\t" << N_comp << "\n";
    std::cout << "Number of droplet fractions:\t" << N_drop_frac << "\n";
    std::cout << "Number of droplet momenta:\t" << _N_drop_mom << "\n\n";
    std::cout << "##########################################\n";


    // if(_N_drop_mom == _N_drop_frac)
    // {
    //     for(int i = 0; i < N_drop_frac; i++)
    //     {
    //         drop_mom_idx[i] = i+N_comp+N_drop_frac;
    //         std::cout << drop_mom_idx[i] << "\n";
    //     }

    // }
    // else
    // {

    // }
}

variables::variables(int _N_var, int _N, int _N_drop_frac, int _N_drop_mom, std::vector<std::vector<double>> const& W_0) : variables(_N_var, _N, W_0)
{

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

    double Q_V = 2*Q_tot/V;

    for(int i = 0; i < N; i++)
    {
        if(msh.x[i] > x_from && msh.x[i] < x_to)
        {
            q[i] = Q_V/2*(1-cos(2*M_PI*(msh.x[i] - x_from)/(x_to-x_from)));
        }
    }

    std::cout << "heat source: " << Q_V << " " << V << "\n";

    
}

void variables::apply_mass_source(double M_tot, double T, double x_from, double x_to, mesh const& msh, std::vector<double> const& comp)
{
    double V = 0;

    for(int i = 1; i < msh.N-1; i++)
    {
        if(msh.x[i] > x_from && msh.x[i] < x_to)
        {
            V += (msh.xf[i] - msh.xf[i-1])*msh.A[i];
        }
    }

    std::vector<double> spec_frac = comp;

    spec_frac[0] = 1;

    double M_V = 2*M_tot/V;

    for(int i = 0; i < N; i++)
    {
        if(msh.x[i] > x_from && msh.x[i] < x_to)
        {
            for(auto k = 0; k < N_comp; k++)
            {
                md[i][k] = M_V/2*(1-cos(2*M_PI*(msh.x[i] - x_from)/(x_to-x_from)))*spec_frac[k];
                q[i] += thermo::cp_mix_comp(comp,T)*T*M_V/2*(1-cos(2*M_PI*(msh.x[i] - x_from)/(x_to-x_from)))*comp[k];
            }
        }
    }
}

void variables::export_to_file(mesh const& msh)
{
    double kappa,r;

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

    std::vector<double> comp(thermo::n_comp);
    stream =  std::ofstream("out/Y.txt");
    for(int i = 0; i < N; i++)
    {
        thermo::composition(comp,W[i]);

        stream << msh.x[i] << " ";
        for(auto const c : comp)
        {
            stream << c << " ";
        }
        stream << "\n";
    }
    stream.close();
    

    stream =  std::ofstream("out/p.txt");

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << thermo::p[i] << "\n";
    }

    stream << "\n";
    stream.close();

    stream =  std::ofstream("out/u.txt");

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << W[i][mom_idx]/W[i][0] << "\n";
    }

    stream << "\n";
    stream.close();

    stream =  std::ofstream("out/T.txt");

    double T_max = 0;

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << thermo::T[i] << "\n";
        T_max = std::max(T_max,thermo::temperature(W[i]));
    }

    stream << "\n";
    stream.close();

    stream =  std::ofstream("out/H.txt");

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << thermo::enthalpy(i,W[i]) << "\n";
    }

    stream << "\n";
    stream.close();

    stream =  std::ofstream("out/H0.txt");

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << thermo::enthalpy_stagnate(i,W[i]) << "\n";
    }

    stream << "\n";
    stream.close();

    stream =  std::ofstream("out/a.txt");

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << thermo::speed_of_sound(i,W[i]) << "\n";
    }

    stream << "\n";
    stream.close();

    stream =  std::ofstream("out/md.txt");

    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << W[i][mom_idx]*msh.A[i] << "\n";
    }

    stream << "\n";
    stream.close();


    for(int k = 0; k < N_drop_frac; k++)
    {
        stream =  std::ofstream("out/X" + std::to_string(k) + ".txt");
        for(int i = 0; i < N; i++)
        {
            stream << msh.x[i] << " " << W[i][3+k] << "\n";
        }
        stream << "\n";
        stream.close();
    }


    stream =  std::ofstream("out/grad.txt");
    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " ";
        for(int k = 0; k < N_var; k++)
        {
            stream << grad[i][k] << " ";
        }
        stream << "\n";
    }
    stream.close();

    stream =  std::ofstream("out/Q_add.txt");
    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " " << q[i] << "\n";
    }
    stream.close();

    stream =  std::ofstream("out/md_add.txt");
    for(int i = 0; i < N; i++)
    {
        stream << msh.x[i] << " ";
        for(int k = 0; k < N_comp; k++)
        {
            stream << md[i][k] << " ";
        }
        stream << "\n";
    }
    stream.close();
}

void variables::export_timestep(double t, mesh const& msh, std::vector<particle> const& particles)
{
    double kappa,r;

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
        stream << msh.x[i] << "\t" << W[i][0] << "\t" << W[i][N_comp]/W[i][0] << "\t" << thermo::p[i] << "\t" << thermo::T[i]
                           << "\t" << particle_values[i][0] << "\t" << particle_values[i][1] << "\t" << particle_values[i][2] << "\t" << particle_values[i][3] << "\n";
    }

}