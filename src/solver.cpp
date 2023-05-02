#include "solver.hpp"
#include "thermodynamics.hpp"
#include <cmath>
#include <iostream>
#include "omp.h"
#include "evaporation.hpp"

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
    var.md[0][2] = 0;

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

        // res[i][var.eng_idx] += dm*33.326e6/dt;
        res[i][var.eng_idx] += dm*30e6/dt;
    }
}

void solver::droplet_transport(std::vector<std::vector<double>>& res, variables& var, mesh const& msh)
{
    static double r;
    static double dm,dr;

    static int N_idx, Frac_idx;
    static int N_comp;
    
    N_comp = var.N_comp;

    auto comp = std::vector<double>{0,0,0};

    for(int i = 1; i < var.N-1; i++)
    {
        if(msh.x[i] < 0.005) continue; // not solving for droplet evaporation near inlet boundary

        dm = evaporation::drop_combustion_steady(i,var.W[i],res[i]);
        var.md[i][2] += dm;
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

double solver::AUSM_wall_mach_number(double M_left, double M_right)
{
    // M+
    double M_wall = 0;

    if(M_left < -1) M_wall += 0;
    else if(M_left > 1) M_wall += M_left;
    else M_wall += 0.25*pow(M_left+1,2);

    //M-
    if(M_right < -1) M_wall += M_right;
    else if(M_right > 1) M_wall += 0;
    else M_wall += -0.25*pow(M_left-1,2);

    return M_wall;
}

double solver::AUSM_wall_pressure(double M_left, double M_right, double p_left, double p_right)
{
    // P+
    double P_wall = 0;

    if(M_left < -1) P_wall += 0;
    else if(M_left > 1) P_wall += p_left;
    else P_wall += 0.25*pow(M_left+1,2)*(2-M_left)*p_left;

    //P-
    if(M_right < -1) P_wall += p_right;
    else if(M_right > 1) P_wall += 0;
    else P_wall += 0.25*pow(M_right-1,2)*(2+M_right)*p_right;

    return P_wall;
}

void solver::AUSM_flux(variables& var, mesh const& msh, parameters const& par)
{
    double M_right, M_left;
    double p_right, p_left;
    double M_wall, p_wall;
    double c;

    M_left = thermo::mach_number(0,var.W[0]);
    p_left = thermo::p[0];

    int cell_idx;

    for(int i = 0; i < var.N_walls; i++)
    {
        M_right = thermo::mach_number(i+1,var.W[i+1]);
        p_right = thermo::p[i+1];

        M_wall = AUSM_wall_mach_number(M_left,M_right);
        p_wall = AUSM_wall_pressure(M_left,M_right,p_left,p_right);

        if(M_wall > 0) cell_idx = i;
        else cell_idx = i+1;

        c = thermo::speed_of_sound(cell_idx,var.W[cell_idx]);

        for(int k = 0; k < variables::N_comp + variables::N_drop_eq; k++)
        {
            var.flux[i][k] = M_wall*c*var.W[cell_idx][k];
        }

        var.flux[i][var.mom_idx] = M_wall*c*var.W[cell_idx][var.mom_idx] + p_wall;
        var.flux[i][var.eng_idx] = M_wall*c*(var.W[cell_idx][var.eng_idx] + thermo::p[cell_idx]);

        M_left = M_right;       // old right is the new left
        p_left = p_right;
    }
}

inline void solver::Euler_flux(int i, std::vector<double>& flux, std::vector<double> const& W)
{
    static double p;
    // static int n_var = flux.size()-2;
    const int n_var = variables::mom_idx;

    p = thermo::p[i];

    flux[0] = W[n_var];

    for(auto idx = 1; idx < variables::N_comp; idx++)
    {
        flux[idx] = W[n_var]*W[idx]/W[0];
    }

    int frac_idx, num_idx, mom_idx;

    for(auto idx = 0; idx < variables::N_drop_frac; idx++)
    {
        frac_idx = idx*2+1;
        num_idx = idx*2;

        flux[num_idx+variables::N_comp] = W[num_idx+variables::N_comp]*(W[variables::drop_mom_idx[idx]]/W[0]);
        flux[frac_idx+variables::N_comp] = W[frac_idx+variables::N_comp]*(W[variables::drop_mom_idx[idx]]/W[0]);
    }

    for(auto idx = 0; idx < variables::N_drop_mom_eq; idx++)
    {
        // mom_idx = variables::mom_idx;
        // mom_idx = idx + variables::N_comp + 2*variables::N_drop_frac;
        mom_idx = variables::drop_mom_idx[idx];
        frac_idx = idx*2+1+variables::N_comp;

        flux[mom_idx] = W[mom_idx]*W[mom_idx]/(W[frac_idx] + 1e-12);
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