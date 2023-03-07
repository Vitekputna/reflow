#include "reflow.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include "lagrange_solver.hpp"
#include "particle.hpp"

extern double kappa;

reflow::reflow(variables& _var, mesh& _msh) : var{_var}, msh{_msh} 
{
    thermo::init(N);    
}

reflow::reflow(int _N, int _N_var) : N{_N}, N_var{_N_var}
{
    msh = mesh(N);
    thermo::init(N);
}

reflow::reflow()
{
    msh = mesh();
    var = variables();
}

reflow::reflow(int _N, int _N_var, std::vector<double> const& init) : N{_N}, N_var{_N_var}
{
    msh = mesh(N);
    var = variables(N_var,N,init);
}

reflow::reflow(int _N, int _N_var, std::vector<std::vector<double>> const& init) : N{_N}, N_var{_N_var}
{
    msh = mesh(N);
    var = variables(N_var,N,init);
}

reflow::reflow(int _N, int _N_var, double from, double to, std::vector<double> const& init) : N{_N}, N_var{_N_var}
{
    msh = mesh(N,from,to);
    var = variables(N_var,N,init);
}

reflow::reflow(int _N, int _N_var, double from, double to) : N{_N}, N_var{_N_var}
{
    msh = mesh(N,from,to);
}

void reflow::initial_conditions(std::vector<double> const& init)
{
    N_var = init.size();
    var = variables(N_var,N,init);
}

void reflow::apply_heat_source(double Q, double x_from, double x_to)
{
    var.apply_heat_source(Q,x_from,x_to,msh);
}

void reflow::apply_mass_source(double M, double T, double x_from, double x_to, std::vector<double> const& comp)
{
    var.apply_mass_source(M, T, x_from, x_to, msh, comp);
}

void reflow::set_boundary(void(*left)(variables&,mesh&,std::vector<double>&), void(*right)(variables&,mesh&,std::vector<double>&))
{
    right_boundary = right;
    left_boundary = left;
}

void reflow::set_boundary(void(*left)(variables&,mesh&,std::vector<double>&), std::vector<double> _left_values,
                    void(*right)(variables&,mesh&,std::vector<double>&), std::vector<double> _right_values)
{
    right_boundary = right;
    left_boundary = left;

    right_values = _right_values;
    left_values = _left_values;
}    

void reflow::spline_geometry(std::vector<std::vector<std::vector<double>>> curves, int n)
{
    msh.cubic(curves,n);
}

void reflow::bump_geometry()
{
    msh.bump();
}

void reflow::refine_mesh(std::vector<std::vector<double>> ref)
{
    thermo::init(N);
    msh.refine(ref);
    N = msh.N-2;
    from = msh.x_from;
    to = msh.x_to;
}

void reflow::add_specie(double r, double kappa, double Mm, std::vector<double> cp_coeff)
{
    specie spec(r,kappa,Mm,cp_coeff);
    thermo_manager.load_specie(spec);
}

void reflow::add_reaction(reaction& R)
{
    chemistry.add_reaction(R);
}

void reflow::export_particles(std::vector<particle>& particles)
{
    auto p_stream = std::ofstream("out/particles.txt");

    for(auto& P : particles)
    {
        if(P.in_use) p_stream << P.x << " " << P.u << " " << P.T << " " << P.r << "\n";
    }
}

void reflow::init_particles(int N_max, int N_particles, int N_per_group)
{
    par_man = particle_manager(N_max,N_per_group,N_particles);
    run_w_particles = true;
}

void reflow::solve()
{
    std::vector<std::vector<double>> res(var.N+2,std::vector<double>(var.N_var,0.0));

    int n = 1;
    double t = 0;
    double dt = 2e-8;
    double t_end = 0.5;
    double residual = 2*max_res;
    double CFL = 0.2;

    auto stream = std::ofstream("out/res.txt");
    stream << "Time [s]\tResidual[...]\n";
    stream.close();

    do
    {
        // update pressure and temperature
        thermo::update(var.W);

        // flow field part  
        solver::reconstruct(var,msh);
        // solver::compute_wall_flux(dt,var,msh,solver::Lax_Friedrichs_flux);
        solver::compute_wall_flux(dt,var,msh,solver::Kurganov_Tadmore);
        solver::compute_cell_res(res,var,msh);
        solver::apply_source_terms(res,var,msh);
        solver::chemical_reactions(dt,res,var,msh);
        // chemistry.solve(dt,res,var,msh);

        // lagrangian particles part
        if(run_w_particles)
        {
            // par_man.particle_inlet(dt*0.18,1e-4,100,0,1000,300);
            if(!(n % 5))  par_man.particle_inlet(5*dt*0.18,1e-4,1.1e-4,60,70,0,0,1000,300);
            // par_man.particle_inlet(dt*0.18,1e-4,1.1e-4,60,70,0,0,1000,300);
            lagrange_solver::update_particles(dt,par_man.particles,var,msh,res);
        }

        // time integration
        solver::Explicit_Euler(var,res,dt);

        // Runtime stuff
        if(!(n % n_res)) 
        {
            stream = std::ofstream("out/res.txt",std::ios_base::app);
            residual = solver::max_residual(res,var,var.eng_idx);
            std::cout << t << " " << dt << " " << residual << "              \r" << std::flush; 
            stream << t << "\t" << solver::max_residual(res,var,0) << "\t" << solver::max_residual(res,var,1) << "\t" << residual << "\n";
            stream.close();
        }

        // Runtime export
        if(!(n % n_exp))
        {
            var.export_to_file(msh);
            // var.export_timestep(t,msh,par_man.particles);
        }

        dt = solver::time_step(var,msh,CFL);

        // boundary
        left_boundary(var,msh,left_values);
        right_boundary(var,msh,right_values);

        t += dt;
        n++;

    // } while (n < 5);
    } while (t < t_end && residual > max_res);

    std::cout << "\r" << std::flush;
    std::cout << "Computation done...                    \n";
    std::cout << "///////////////////////////////////////\n";
    std::cout << "total steps[/] \t|end time[s] \t|final residual[/]\n";
    std::cout << n << "\t\t " << t << "\t\t " << residual << "\n";
    std::cout << "///////////////////////////////////////\n";
    std::cout << "Number of particles: " << par_man.particles.size() << "\n";

    solver::reconstruct(var,msh);

    reflow::export_particles(par_man.particles);
    var.export_to_file(msh);
    msh.export_to_file();
}