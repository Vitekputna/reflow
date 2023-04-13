#include "solver.hpp"
#include "thermodynamics.hpp"
#include <cmath>
#include <iostream>
#include "omp.h"

void solver::compute_wall_flux(double dt, variables& var, mesh const& msh,
                               void(*flux)(variables&,mesh const&,parameters const&))
{
    flux(var, msh, parameters(msh.dx_min,dt,0.7));
}

void solver::compute_exact_flux(variables& var)
{
    for(int i = 0; i < var.N; i++)
    {
        Euler_flux(i,var.exact_flux[i],var.W[i]);
    }
}

void solver::compute_cell_res(std::vector<std::vector<double>>& res,variables& var, mesh const& msh)
{
    for(int i = 1; i < var.N-1; i++)
    {
        for(int k = 0; k < var.N_var; k++)
        {
            res[i][k] = -1/(msh.A[i]*(msh.xf[i] - msh.xf[i-1]))*(msh.Af[i]*var.flux[i][k] - msh.Af[i-1]*var.flux[i-1][k]);
        }
    }
}

void solver::apply_source_terms(std::vector<std::vector<double>>& res, variables& var, mesh const& msh)
{
    for(int i = 1; i < var.N-1; i++)
    {
        for(int k = 0; k < var.N_comp; k++)
        {
            // res[i][k] += var.md[i][k];
            var.md[i][k] = 0;
        }
        res[i][var.mom_idx] += ((msh.Af[i]-msh.Af[i-1])/(msh.xf[i]-msh.xf[i-1]))/msh.A[i]*thermo::p[i];
        res[i][var.eng_idx] += var.q[i];
    }
}

void solver::chemical_reactions(double dt,std::vector<std::vector<double>>& res, variables& var, mesh const& msh)
{
    double dm, V;
    for(int i = 1; i < var.N-1; i++)
    {
        dm = std::max(0.0,std::min(var.W[i][2],var.W[i][1]/6.6));

        res[i][1] += -6.6*dm/dt;    // Oxydizer
        res[i][2] += -dm/dt;        // Fuel

        res[i][var.eng_idx] += dm*33.326e6/dt;
        // res[i][var.eng_idx] += dm*20e6/dt;
    }
}

void solver::droplet_transport(std::vector<std::vector<double>>& res, variables& var, mesh const& msh)
{
    static double r;
    static double dm,dr;

    static int N_idx, Frac_idx;
    static int N_comp;
    
    N_comp = var.N_comp;

    for(int i = 1; i < var.N-1; i++)
    {
        if(msh.x[i] < 0.005) continue; // not solving for droplet evaporation near inlet boundary

        const double T_coeff = log(1 + 1e-3*std::max(0.0,thermo::T[i] - 200));

        var.md[i][2] = 0;

        for(int j = 0; j < var.N_drop_frac; j += 2)
        {
            N_idx = N_comp + j;
            Frac_idx = N_comp + j + 1;

            r = std::pow(3*var.W[i][Frac_idx]/(4*var.W[i][N_idx]*3.14159*700),0.3333);
            if(var.W[i][N_idx] == 0) r = 0;

            dm = std::max(0.0,3*var.W[i][N_idx]*r*T_coeff);

            res[i][Frac_idx] -= dm;
            res[i][2] += dm;
            var.md[i][2] += dm;
            res[i][var.eng_idx] += dm*thermo::enthalpy(thermo::T[i],std::vector<double>{0,0,1});
        }
    }
}

void solver::reconstruct(variables& var, mesh const& msh)
{
    double phi;
    for(auto i = 1; i < var.N-1; i++)
    {
        for(auto k = 0; k < var.N_var; k++)
        {
            // phi = minmod(var.W[i+1][k] - var.W[i][k] , var.W[i][k] - var.W[i-1][k]);
            phi = van_albada(var.W[i+1][k] - var.W[i][k] , var.W[i][k] - var.W[i-1][k]);
            // phi = van_leer(var.W[i+1][k] - var.W[i][k] , var.W[i][k] - var.W[i-1][k]);

            var.grad[i][k] = phi*(var.W[i+1][k] - var.W[i-1][k])/(msh.x[i+1] - msh.x[i-1]);
        }
    }
}   

inline double solver::minmod(double a, double b)
{
    double r = a/b;

    if(r > 0 && b != 0)
    {
        return std::min(2/(1+r),2*r/(1+r));
    }
    else
    {
        return 0.0; 
    }
}

inline double solver::van_albada(double a, double b)
{
    double r = a/b;

    if(r > 0 && b != 0)
    {
        return 2*r/(r*r+1);
    }
    else
    {
        return 0.0;
    }
}

inline double solver::van_leer(double a, double b)
{
    double r = a/b;

    if(a*b > 0)
    {
        return 4*r/(r+1)/(r+1);
    }
    else
    {
        return 0.0;
    }
}

void solver::Lax_Friedrichs_flux(variables& var, mesh const& msh, parameters const& par)
{
    compute_exact_flux(var);

    for(int i = 0; i < var.N_walls; i++)
    {
        for(int k = 0; k < var.N_var; k++)
        {
            var.flux[i][k] = 0.5*(var.exact_flux[i+1][k] + var.exact_flux[i][k]) - par.eps*0.5*par.dx/par.dt*(var.W[i+1][k] - var.W[i][k]);
        }
    }
}

void solver::HLL_flux(variables& var, mesh const& msh, parameters const& par)
{
    static double sr,sl;
    static double cr,cl,ur,ul;

    static bool right, left, center;

    compute_exact_flux(var);

    for(int i = 0; i < var.N_walls; i++)
    {
        ur = var.W[i+1][var.mom_idx]/var.W[i+1][0];
        ul = var.W[i][var.mom_idx]/var.W[i][0];

        cr = thermo::speed_of_sound(i+1,var.W[i+1]);
        cl = thermo::speed_of_sound(i,var.W[i]);

        sr = std::max(ur+cr,ul+cl);
        sl = std::min(ur-cr,ul-cl);

        if(sl >= 0)
        {
            left = true;
            right = false;
            center = false;
        }
        else if(sr < 0)
        {
            right = true;
            left = false;
            center = false;
        }
        else
        {
            center = true;
            right = false;
            left = false;
        }

        for(int k = 0; k < var.N_var; k++)
        {
            var.flux[i][k] = left*var.exact_flux[i][k] + center*((sr*var.exact_flux[i][k] - sl*var.exact_flux[i+1][k] + sr*sl*(var.W[i+1][k]-var.W[i][k]))/(sr-sl))
                            +right*var.exact_flux[i+1][k];
        }
    }

}

void solver::reconstructed_flux(int idx, std::vector<double>& flux, std::vector<double> W, std::vector<double> const& grad, double dx)
{
    for(auto i = 0; i < W.size(); i++)
    {
        W[i] += grad[i]*dx;
    }

    Euler_flux(idx, flux,W);
}

void solver::reconstructed_wave_speed(int i, std::vector<double>& a, std::vector<double> W, std::vector<double> const& grad, double dx)
{
    static int size;
    static double c;

    size = W.size();
    c = thermo::speed_of_sound(i,W);

    for(auto i = 0; i < size; i++)
    {
        W[i] += grad[i]*dx;
    }

    // c = thermo::speed_of_sound(i,W);

    a[0] = W[size-2]/W[0] - c;

    for(auto i = 1; i < size-2; i++)
    {
        a[i] = W[size-2]/W[0];
    }

    a[size-1] = W[size-2]/W[0] + c;
}

void solver::Kurganov_Tadmore(variables& var, mesh const& msh, parameters const& par)
{
    std::vector<double> fr(var.N_var,0.0);
    std::vector<double> fl(var.N_var,0.0);
    std::vector<double> ar(var.N_var,0.0);
    std::vector<double> al(var.N_var,0.0);

    double ul, ur, a;

    for(int i = 0; i < var.N_walls; i++)
    {
        reconstructed_flux(i,fl,var.W[i],var.grad[i],msh.xf[i]-msh.x[i]);
        reconstructed_flux(i+1,fr,var.W[i+1],var.grad[i+1],msh.xf[i]-msh.x[i+1]);

        reconstructed_wave_speed(i,al,var.W[i],var.grad[i],msh.xf[i]-msh.x[i]);
        reconstructed_wave_speed(i+1,ar,var.W[i+1],var.grad[i+1],msh.xf[i]-msh.x[i+1]);

        for(int k = 0; k < var.N_var; k++)
        {
            ul = var.W[i][k] + var.grad[i][k]*(msh.xf[i]-msh.x[i]);
            ur = var.W[i+1][k] + var.grad[i+1][k]*(msh.xf[i]-msh.x[i+1]);

            var.flux[i][k] = 0.5*(fr[k]+fl[k]) - std::max(std::abs(ar.back()),std::abs(al.back()))/2*(ur-ul);
        }
    }
}

inline void solver::Euler_flux(int i, std::vector<double>& flux, std::vector<double> const& W)
{
    static double p;
    static int n_var = flux.size()-2;

    p = thermo::p[i];

    flux[0] = W[n_var];

    for(auto idx = 1; idx < variables::N_comp; idx++)
    {
        flux[idx] = W[n_var]*W[idx]/W[0];
    }

    for(auto idx = 0; idx < variables::N_drop_frac; idx++)
    {
        flux[idx+variables::N_comp] = W[idx+variables::N_comp]*(W[variables::drop_mom_idx[idx]]/W[0]);
    }

    flux[n_var] = W[n_var]*W[n_var]/W[0] + p;
    flux[n_var+1] = (W[n_var+1] + p)*W[n_var]/W[0];
}

void solver::Explicit_Euler(variables& var, std::vector<std::vector<double>>& res, double dt)
{
    for(int i = 1; i < var.N-1; i++)
    {
        for(int k = 0; k < var.N_var; k++)
        {
            var.W[i][k] += dt*res[i][k];
        }
    }
}

void solver::Adams_Bashforth(variables& var,std::vector<std::vector<double>>& res_new,std::vector<std::vector<double>>& res_old, double dt)
{
    for(int i = 0; i < var.N; i++)
    {
        for(int k = 0; k < var.N_var; k++)
        {
            var.W[i][k] += dt*(1.5*res_new[i][k] - 0.5*res_old[i][k]);
        }
    }
}

double solver::time_step(variables const& var, mesh const& msh, double CFL)
{
    double dt = 1e10;

    for(int i = 0; i < msh.N; i++)
    {
        dt = std::min(dt, msh.dx_min/(std::abs(var.W[i][var.mom_idx]/var.W[i][0]) + thermo::speed_of_sound(i,var.W[i])));
    }

    return CFL*dt;
}

double solver::max_residual(std::vector<std::vector<double>> const& res, variables const& var, int res_idx) 
{
    double max_res = 0;

    for(unsigned int i = 1; i < res.size()-1; i++)
    {
        max_res = std::max(max_res, std::abs(res[i][res_idx]));
    }

    return max_res;
}