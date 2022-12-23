#include "reflow.hpp"
#include <vector>

reflow::reflow(variables& _var, mesh& _msh) : var{_var}, msh{_msh} {}

reflow::reflow(int N, int N_var, std::vector<double> init)
{
    msh = mesh(N);
    var = variables(N_var,N,init);
}

void reflow::solve()
{
    double t = 0;
    double dt = 1e-4;

    double t_end = 10*dt;

    std::vector<std::vector<double>> res(var.N+2,std::vector<double>(var.N_var,0.0));

    do
    {
        solver::compute_wall_flux(var,msh,solver::Lax_Friedrichs_flux);
        solver::compute_cell_res(res,var,msh);
        solver::Explicit_Euler(var,res,dt);

        t += dt;
    } while (t < t_end);

    var.export_to_file("out/test.txt");
}