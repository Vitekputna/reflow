#pragma once
#include <vector>
#include "mesh.hpp"

namespace init
{
    std::vector<std::vector<double>> step(int N, int N_var, std::vector<double> const& L, std::vector<double> const& R);

    std::vector<double> flow(int N_var, double p, double T, double u, std::vector<double> const& comp);

    std::vector<double> flow_droplets(int N_var, double p, double T, double u, std::vector<double> const& comp,
                                         std::vector<double> drp_frac,
                                         std::vector<double> drp_count,
                                         std::vector<double> drp_mom);

    std::vector<std::vector<double>> nozzle(int N, int N_var,double md, double T0, double p0, double p2, double L_chamber,
                                         std::vector<double> const comp, mesh const& msh);
                                         
    std::vector<std::vector<double>> nozzle_droplets(int N, int N_var,double md, double T0, double p0, double p2, double L_chamber,
                                         std::vector<double> const comp, mesh const& msh,
                                         std::vector<double> drp_frac,
                                         std::vector<double> drp_count,
                                         std::vector<double> drp_mom);
}