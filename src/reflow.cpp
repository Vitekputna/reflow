#include "reflow.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <omp.h>
#include "lagrange_solver.hpp"
#include "particle.hpp"

extern double kappa;

reflow::reflow(variables& _var, mesh& _msh) : var{_var}, msh{_msh} 
{
    thermo::init(N);  
    variable_init();  
}

reflow::reflow(int _N, int _N_var) : N{_N}, N_var{_N_var}
{
    msh = mesh(N);
    thermo::init(N);
    variable_init();
}

reflow::reflow()
{
    msh = mesh();
    var = variables();
    variable_init();
}

reflow::reflow(int _N, int _N_var, std::vector<double> const& init) : N{_N}, N_var{_N_var}
{
    msh = mesh(N);
    var = variables(N_var,N,init);
    variable_init();
}

reflow::reflow(int _N, int _N_var, std::vector<std::vector<double>> const& init) : N{_N}, N_var{_N_var}
{
    msh = mesh(N);
    var = variables(N_var,N,init);
    variable_init();
}

reflow::reflow(int _N, int _N_var, double from, double to, std::vector<double> const& init) : N{_N}, N_var{_N_var}
{
    msh = mesh(N,from,to);
    var = variables(N_var,N,init);
    variable_init();
}

reflow::reflow(int _N, int _N_var, double from, double to) : N{_N}, N_var{_N_var}
{
    msh = mesh(N,from,to);
    variable_init();
}

void reflow::variable_init()
{
    n_dt = 2;
    n_res = 500;
    n_exp = 2000;

    max_res = 1000;
    t_end = 1;
    CFL = 0.1;
}

void reflow::set_flux_func_vectors()
{
    // fluid flux index
    for(int i = 0; i < variables::N_comp; i++)
    {
        fluid_flux_idx.push_back(i);
    }
    fluid_flux_idx.push_back(variables::mom_idx);
    fluid_flux_idx.push_back(variables::eng_idx);

    // dispersed flux idx
    for(int i = 0; i < variables::N_drop_eq; i++)
    {
        dispersed_flux_idx.push_back(variables::N_comp+i);
    }
}

void reflow::initial_conditions(std::vector<double> const& init)
{
    N_var = init.size();
    var = variables(N_var,N,init);
    set_flux_func_vectors();
}

void reflow::initial_conditions(std::vector<std::vector<double>> const& init)
{
    N_var = init[0].size();
    var = variables(N_var,N,init);
    set_flux_func_vectors();
}

void reflow::initial_conditions(int N_drop, bool drop_momenta, std::vector<double> const& init)
{
    N_var = init.size();
    var = variables(N_var,N,N_drop,drop_momenta,init);
    euler_particles = true;
    set_flux_func_vectors();
}

void reflow::initial_conditions(int N_drop, bool drop_momenta, std::vector<std::vector<double>> const& init)
{
    N_var = init[0].size();
    var = variables(N_var,N,N_drop,drop_momenta,init);
    euler_particles = true;
    set_flux_func_vectors();
}

void reflow::initial_conditions(int N_drop, bool drop_momenta, bool drop_energy, std::vector<double> const& init)
{
    N_var = init.size();
    var = variables(N_var,N,N_drop,drop_momenta,drop_energy,init);
    euler_particles = true;
    set_flux_func_vectors();
}

void reflow::initial_conditions(int N_drop, bool drop_momenta, bool drop_energy, std::vector<std::vector<double>> const& init)
{
    N_var = init[0].size();
    var = variables(N_var,N,N_drop,drop_momenta,drop_energy,init);
    euler_particles = true;
    set_flux_func_vectors();
}

void reflow::apply_heat_source(double Q, double x_from, double x_to)
{
    var.apply_heat_source(Q,x_from,x_to,msh);
}

void reflow::apply_mass_source(double M, double T, double x_from, double x_to, std::vector<double> const& comp)
{
    std::cout << "Aplying mass source:\n";
    std::cout << "Total flux:\t\t" << M << "\n";
    std::cout << "Temperature:\t\t" << T << "\n";
    std::cout << "Source from:\t\t" << x_from << "\n";
    std::cout << "Source to:\t\t" << x_to << "\n";
    std::cout << "##########################################\n";

    var.apply_mass_source(M, T, x_from, x_to, msh, comp);
}

void reflow::add_boundary_function(void(*func)(variables&,mesh&,std::vector<double>&),std::vector<double> values)
{
    boundary_func_vec.push_back(func);
    boundary_values_vec.push_back(values);
}

void reflow::apply_boundary_conditions()
{
    for(auto i = 0; i < boundary_func_vec.size(); i++)
    {
        boundary_func_vec[i](var,msh,boundary_values_vec[i]);
    }
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
    msh.refine(ref);
    N = msh.N;

    thermo::init(N);
    from = msh.x_from;
    to = msh.x_to;
}

void reflow::add_specie(double r, double kappa, double Mm, std::vector<double> cp_coeff)
{
    specie spec(r,kappa,Mm,cp_coeff);
    thermo_manager.load_specie(spec);
    // n_comp++;
}

void reflow::add_specie(double r, double kappa, double Mm, std::vector<double> cp_coeff, std::vector<double> k_coeff, std::vector<double> mu_coeff)
{
    specie spec(r,kappa,Mm,cp_coeff,k_coeff,mu_coeff);
    thermo_manager.load_specie(spec);
    // n_comp++;
}

void reflow::add_reaction(reaction& R)
{
    chemistry.add_reaction(R);
}

void reflow::export_data()
{
    var.export_to_file(export_path,msh,par_man.particles);
}

void reflow::export_mesh()
{
    msh.export_to_file(export_path);
}

void reflow::export_particles()
{
    auto p_stream = std::ofstream("out/Lagrangian.txt");

    auto particles = this->par_man.particles;

    for(auto& P : particles)
    {
        if(P.in_use) p_stream << P.x << " " << P.u << " " << P.T << " " << P.r << " " << P.N << " " << P.M << " " << P.rho << "\n";
    }
}

void reflow::init_particles(int N_max, int N_particles, int N_per_group)
{
    par_man = particle_manager(N_max,N_per_group,N_particles);
    lagrange_particles = true;
}

void reflow::add_lagrangian_mono_particles(double specie_idx, double mass_flux, double rho, double r, double x,
                                           double u, double T, double T_boil, double vap_heat, double C)
{
    auto parameters = std::vector<double>{mass_flux,r,u,x,rho,T,T_boil,vap_heat,C,specie_idx};
    par_man.add_monodispersion(parameters);
}

void reflow::add_lagrangian_unif_particles(double specie_idx, double mass_flux, double rho, double r_mean, double r_var, 
                                           double x, double u, double T, double T_boil, double vap_heat, double C)
{
    auto parameters = std::vector<double>{mass_flux,r_mean,r_var,u,x,rho,T,T_boil,vap_heat,C};
    par_man.add_uniform(parameters);
}

void reflow::add_lagrangian_norm_particles(double specie_idx, double mass_flux, double rho, double r_mean, double r_var, 
                                           double x, double u, double T, double T_boil, double vap_heat, double C)
{
    auto parameters = std::vector<double>{mass_flux,r_mean,r_var,u,x,rho,T,T_boil,vap_heat,C};
    par_man.add_normal(parameters);
}

void reflow::apply_lagrangian_particle_inlet(double dt)
{
    par_man.apply_boundary(dt);
}

std::vector<std::vector<double>> reflow::read_files(std::string path)
{
    std::cout << "Reading external data...\n";

    std::vector<double> x_data;
    std::vector<std::vector<double>> y_data;

    for (int i = 0; true; ++i) {
        std::string filename = path + "W" + std::to_string(i) + ".txt";

        std::ifstream file(filename);
        if (!file.is_open()) {
            break; // exit loop if file doesn't exist
        }

        y_data.push_back(std::vector<double>{});

        double x, y;
        while (file >> x >> y) {
            
            if(i == 0) x_data.push_back(x);             // only read first column for the first file
            y_data[i].push_back(y);                     // read data
        }
    }

    int N_var_data = y_data.size();
    int N_data = x_data.size();

    std::cout << "Read " << N_var_data << " datafiles of " << N_data << " elements.\n";

    if(N != N_data)
    {
        std::cout << "Loaded data is not compatible with defined mesh.\n";
        std::cout << "exiting...\n";

        return std::vector<std::vector<double>>{{}};
    }

    N_var = N_var_data;

    // transpose data
    std::vector<std::vector<double>> data(N,std::vector<double>(N_var,0.0));

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N_var; j++)
        {
            data[i][j] = y_data[j][i];
        }
    }

    return data;
}

void reflow::load_old_data(std::string path, int _n_comp)
{
    std::cout << "#########################################\n";
    std::cout << "Loading old data for initialization...\n";
    std::cout << "Path to data: " << path << "\n";
    std::cout << "Declared number of components: " << _n_comp << "\n";

    auto init_data = read_files(path);

    var = variables(N_var,N,init_data);
    set_flux_func_vectors();

    // n_comp = _n_comp;
    std::cout << "##########################################\n";
}

void reflow::load_old_data(std::string path, int _n_comp, int _n_drop_frac, bool droplet_momenta)
{
    auto init_data = read_files(path);

    int n_comp = _n_comp;
    int n_drop_frac = _n_drop_frac;
    int n_drop_mom = n_drop_frac;

    std::cout << N_var << " " << N << " " << n_drop_frac << " " << n_drop_mom << "\n";

    var = variables(N_var,N,n_drop_frac,true,init_data);
    set_flux_func_vectors();
}

void reflow::load_old_data(std::string path, int _n_comp, int _n_drop_frac, bool droplet_momenta, bool droplet_energy)
{
    auto init_data = read_files(path);

    int n_comp = _n_comp;
    int n_drop_frac = _n_drop_frac;
    int n_drop_mom = n_drop_frac;
    int n_drop_eng = n_drop_frac;

    std::cout << N_var << " " << N << " " << n_drop_frac << " " << n_drop_mom << "\n";

    var = variables(N_var,N,n_drop_frac,true,true,init_data);
    set_flux_func_vectors();
}

bool reflow::maximum_time(double T, double res)
{
    return T < t_end;
}

void reflow::set_numThreads(const int _numThreads)
{
    numThreads = _numThreads;
}

void reflow::set_export_path(std::string path)
{
    std::cout << "##########################################\n";
    std::cout << "Setting export path to: " + path +"\n";
    export_path = path;
}

void reflow::divide_data_parallel()
{
    std::cout << "##########################################\n";
    std::cout << "Running on: " << numThreads << " threads.\n";


    // divide data into threads
    cell_startIndices = std::vector<int>(numThreads);
    cell_endIndices = std::vector<int>(numThreads);

    cell_startIndices_full = std::vector<int>(numThreads);
    cell_endIndices_full = std::vector<int>(numThreads);

    wall_startIndices = std::vector<int>(numThreads);
    wall_endIndices = std::vector<int>(numThreads);

    int cell_chunkSize = this->var.N / numThreads;
    int cell_remaining = this->var.N % numThreads;

    int wall_chunkSize = this->var.N_walls / numThreads;
    int wall_remaining = this->var.N_walls % numThreads;

    // Calculate the start and end indices for each thread
    for (int i = 0; i < numThreads; ++i) {
        cell_startIndices_full[i] = i * cell_chunkSize;
        cell_endIndices_full[i] = (i == numThreads - 1) ? (cell_startIndices_full[i] + cell_chunkSize + cell_remaining - 1) : (cell_startIndices_full[i] + cell_chunkSize - 1);
        std::cout << cell_startIndices_full[i] << " " << cell_endIndices_full[i] << "\n";
    }

    for (int i = 0; i < numThreads; ++i) {
        wall_startIndices[i] = i * wall_chunkSize;
        wall_endIndices[i] = (i == numThreads - 1) ? (wall_startIndices[i] + wall_chunkSize + wall_remaining - 1) : (wall_startIndices[i] + wall_chunkSize - 1);
        std::cout << wall_startIndices[i] << " " << wall_endIndices[i] << "\n";
    }

    cell_startIndices = cell_startIndices_full;
    cell_endIndices = cell_endIndices_full;
    cell_startIndices[0] += 1;
    cell_endIndices.back() -= 1;

    std::cout << cell_startIndices[0] << " " << cell_endIndices.back() << "\n";

    std::cout << "##########################################\n";
}

void reflow::solve_parallel(double _t_end, double _max_residual, double _CFL)
{
    t_end = _t_end;
    max_res = _max_residual;
    CFL = _CFL;

    std::vector<std::vector<double>> res(var.N+2,std::vector<double>(var.N_var,0.0));

    int n = 1;
    double t = 0;
    double dt = 2e-8;
    double residual = 2*max_res;
    bool RUN_FLAG = true;

    auto stream = std::ofstream("out/res.txt");
    stream << "Time [s]\tResidual[...]\n";
    stream.close();

    std::cout << "Starting simulation...\n";
    std::cout << "##########################################\n";
    std::cout << "time[s]\ttime step[s]\tresidual[]\n";
    std::cout << "##########################################\n";

    const int cell_from = 0;
    const int cell_to = var.N-1;

    const int wall_from = 0;
    const int wall_to = var.N_walls-1;

    divide_data_parallel();
    #pragma omp parallel shared(dt,var,msh,res,t,n,residual,stream,t_end) private(RUN_FLAG) num_threads(numThreads)
    {
        int threadID = omp_get_thread_num();
        do
        {
            // update pressure and temperature
            // thermo::update(var.W,cell_startIndices_full[threadID],cell_endIndices_full[threadID]);
            if(threadID == 0)
            {
                thermo::update(var.W);

                if(reconstruct)
                {
                    solver::reconstruct(var,msh,cell_from+1,cell_to-1);
                    solver::reconstruct_pressure(var,msh,cell_from+1,cell_to-1);    
                }
            }

            // #pragma omp barrier
            // if(reconstruct)
            // {
            //     solver::reconstruct(var,msh,cell_startIndices[threadID],cell_endIndices[threadID]);
            //     solver::reconstruct_pressure(var,msh,cell_startIndices[threadID],cell_endIndices[threadID]);    
            // }

            solver::compute_exact_flux(var,cell_startIndices_full[threadID],cell_endIndices_full[threadID]);

            #pragma omp barrier
            
            solver::compute_wall_flux(dt,var,msh,fluid_flux,wall_startIndices[threadID],wall_endIndices[threadID],fluid_flux_idx);
            solver::compute_wall_flux(dt,var,msh,dispersed_flux,wall_startIndices[threadID],wall_endIndices[threadID],dispersed_flux_idx);

            #pragma omp barrier

            solver::compute_cell_res(res,var,msh,cell_startIndices[threadID],cell_endIndices[threadID]);
            solver::apply_source_terms(res,var,msh,cell_startIndices[threadID],cell_endIndices[threadID]);
            // solver::chemical_reactions(dt,res,var,msh,cell_startIndices[threadID],cell_endIndices[threadID]);
            solver::droplet_transport(res,var,msh,cell_startIndices[threadID],cell_endIndices[threadID]);

            // lagrangian particles part
            if(lagrange_particles)
            {
                if(!(n % 5)) apply_lagrangian_particle_inlet(5*dt);
                lagrange_solver::update_particles(dt,par_man.particles,var,msh,res);
            }

            // time integration
            #pragma omp barrier
            solver::Explicit_Euler(var,res,dt,cell_startIndices[threadID],cell_endIndices[threadID]);

            #pragma omp barrier

            if(threadID == 0)
            {
                // Runtime stuff
                if(!(n % n_res)) 
                {
                    stream = std::ofstream(export_path+"res.txt",std::ios_base::app);
                    residual = solver::max_residual(res,var,var.eng_idx);

                    std::cout << "\b\b\r";
                    std::cout << t << " " << dt << " " << residual << "\t" << par_man.N << std::flush; 

                    stream << t << "\t" << solver::max_residual(res,var,0) << "\t" << solver::max_residual(res,var,var.mom_idx) << "\t" << residual << "\n";
                    stream.close();

                    
                }

                // Runtime export
                if(!(n % n_exp))
                {
                    var.export_to_file(export_path, msh, par_man.particles);
                    if(lagrange_particles) export_particles();
                    // var.export_timestep(t,msh,par_man.particles);
                }
            
                dt = solver::time_step(var,msh,CFL);

                // boundary
                apply_boundary_conditions();

                t += dt;
                n++;
            }

            RUN_FLAG = maximum_time(t,residual);

            if(!RUN_FLAG)
            {
                std::cout << threadID << " done\n";
            }

        } while (RUN_FLAG);
        // } while (n < 100);  
    }

    std::cout << "\r\n" << std::flush;
    std::cout << "Simulation done...\n";
    std::cout << "##########################################\n";
    std::cout << "total steps[/] \t|end time[s] \t|final residual[/]\n";
    std::cout << n << "\t\t " << t << "\t\t " << residual << "\n";
    std::cout << "##########################################\n";
    std::cout << "Number of particles: " << par_man.particles.size() << "\n";

    var.export_to_file(export_path,msh,par_man.particles);
    msh.export_to_file(export_path);
}

void reflow::solve(double _t_end, double _max_residual, double _CFL)
{
    t_end = _t_end;
    max_res = _max_residual;
    CFL = _CFL;

    std::vector<std::vector<double>> res(var.N+2,std::vector<double>(var.N_var,0.0));

    int n = 1;
    double t = 0;
    double dt = 2e-8;
    double residual = 2*max_res;
    bool RUN_FLAG = true;

    auto stream = std::ofstream("out/res.txt");
    stream << "Time [s]\tResidual[...]\n";
    stream.close();

    std::cout << "Starting simulation...\n";
    std::cout << "##########################################\n";
    std::cout << "time[s]\ttime step[s]\tresidual[]\n";
    std::cout << "##########################################\n";
    
    int threadID = 0;

    const int cell_from = 0;
    const int cell_to = var.N-1;

    const int wall_from = 0;
    const int wall_to = var.N_walls-1;

    do
    {
        // update pressure and temperature
        thermo::update(var.W);
        
        if(reconstruct)
        {
            solver::reconstruct(var,msh,cell_from+1,cell_to-1);
            solver::reconstruct_pressure(var,msh,cell_from+1,cell_to-1);    
        }

        solver::compute_exact_flux(var,0,cell_to);
        
        solver::compute_wall_flux(dt,var,msh,fluid_flux,wall_from,wall_to,fluid_flux_idx);
        if(euler_particles) solver::compute_wall_flux(dt,var,msh,dispersed_flux,wall_from,wall_to,dispersed_flux_idx);

        solver::compute_cell_res(res,var,msh,cell_from+1,cell_to-1);
        solver::apply_source_terms(res,var,msh,cell_from+1,cell_to-1);
        solver::chemical_reactions(dt,res,var,msh,cell_from+1,cell_to-1);
        solver::droplet_transport(res,var,msh,cell_from+1,cell_to-1);

        // lagrangian particles part
        if(lagrange_particles)
        {
            // apply_lagrangian_particle_inlet(dt);
            if(!(n % 5)) apply_lagrangian_particle_inlet(5*dt);
            
            lagrange_solver::update_particles(dt,par_man.particles,var,msh,res);
        }

        // time integration
        solver::Explicit_Euler(var,res,dt,cell_from+1,cell_to-1);

        // Runtime stuff
        if(!(n % n_res)) 
        {
            stream = std::ofstream(export_path+"res.txt",std::ios_base::app);
            residual = solver::max_residual(res,var,var.eng_idx);

            std::cout << "\b\b\r";
            std::cout << t << " " << dt << " " << residual << "\t" << par_man.N << std::flush; 

            stream << t << "\t" << solver::max_residual(res,var,0) << "\t" << solver::max_residual(res,var,var.mom_idx) << "\t" << residual << "\n";
            stream.close();
        }

        // Runtime export
        if(!(n % n_exp))
        {
            var.export_to_file(export_path,msh,par_man.particles);
            if(lagrange_particles) export_particles();
            // var.export_timestep(t,msh,par_man.particles);
        }
    
        dt = solver::time_step(var,msh,CFL);

        // boundary
        apply_boundary_conditions();

        t += dt;
        n++;

        RUN_FLAG = maximum_time(t,residual);

    } while (RUN_FLAG);
    // } while (n < 100);  

    std::cout << "\r\n" << std::flush;
    std::cout << "Simulation done...\n";
    std::cout << "##########################################\n";
    std::cout << "total steps[/] \t|end time[s] \t|final residual[/]\n";
    std::cout << n << "\t\t " << t << "\t\t " << residual << "\n";
    std::cout << "##########################################\n";
    std::cout << "Number of particles: " << par_man.particles.size() << "\n";

    var.export_to_file(export_path,msh,par_man.particles);
    msh.export_to_file(export_path);
}