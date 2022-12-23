#include "variables.hpp"
#include <vector>
#include <string>
#include <fstream>

#include <iostream>

variables::variables(){}

variables::variables(int _N_var, int _N) : N_var{_N_var}, N{_N+2}, N_walls{_N+1}
{
    W = std::vector<std::vector<double>>(N,std::vector<double>(N_var,0.0));
    flux = std::vector<std::vector<double>>(N_walls,std::vector<double>(N_var,0.0));
    exact_flux = std::vector<std::vector<double>>(N,std::vector<double>(N_var,0.0));
}

variables::variables(int _N_var, int _N, std::vector<double> W_0) : variables(_N_var,_N)
{
    for(auto& wi : W)
    {
        for(int k = 0; k < N_var; k++)
        {
            wi[k] = W_0[k];
        }
    }
}

void variables::export_to_file(std::string path)
{
    auto stream = std::ofstream();

    for(int k = 0; k < N_var; k++)
    {
        stream =  std::ofstream("out/W" + std::to_string(k) + ".txt");

        for(int i = 0; i < N; i++)
        {
            stream << W[i][k] << "\n";
        }

        stream << "\n";
        stream.close();
    }
}