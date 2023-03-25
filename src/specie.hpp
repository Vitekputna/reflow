#pragma once
#include <cmath>
struct specie
{
    double r;
    double kappa;
    double Mm;
    double a,b,c,d,e,f;

    specie(double _r, double _kappa, double _Mm, std::vector<double> cp_coeff) : r{_r}, kappa{_kappa}, Mm{_Mm} 
    {
        a = cp_coeff[0];
        b = cp_coeff[1];
        c = cp_coeff[2];
        d = cp_coeff[3];
        e = cp_coeff[4];
        f = cp_coeff[5];
    };

    double cp(double T)
    {
        T = std::min(T,4000.0);
        return a + b*T + c*pow(T,2) + d*pow(T,3) + e*pow(T,4) + f*pow(T,5);
    }

    double h(double T)
    {
        T = std::min(T,4000.0);
        return a*T + b*pow(T,2)/2 + c*pow(T,3)/3 + d*pow(T,4)/4 + e*pow(T,5)/5 + f*pow(T,6)/6;
    }
};