#pragma once
#include "variables.hpp"
#include "mesh.hpp"

namespace boundary
{
    // fundamental
    void _zero_gradient_r(variables& var, int idx, double value);
    void _set_gradient_r(variables& var, int idx, double value);
    void _set_value_r(variables& var, int idx, double value);

    void _zero_gradient_l(variables& var, int idx, double value);
    void _set_gradient_l(variables& var, int idx, double value);
    void _set_value_l(variables& var, int idx, double value);

    // abstractions
    void zero_gradient_l(variables& var, mesh& msh, std::vector<double>& values);
    void zero_gradient_r(variables& var, mesh& msh, std::vector<double>& values);
    void set_value_l(variables& var, mesh& msh, std::vector<double>& values);
    void mass_flow_inlet(variables& var, mesh& msh, std::vector<double>& values);
    void subsonic_inlet(variables& var, mesh& msh, std::vector<double>& values);
    void subsonic_outlet(variables& var, mesh& msh, std::vector<double>& values);
    void supersonic_outlet(variables& var, mesh& msh, std::vector<double>& values);

    void quiscent_dropplet_inlet(variables& var, mesh& msh, std::vector<double>& values);
}