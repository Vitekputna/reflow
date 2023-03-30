#pragma once
#include <cmath>
#include <vector>
struct specie
{
    double r;
    double kappa;
    double Mm;
    double a,b,c,d,e,f;
    double a1,b1,c1,d1,e1,f1;

    specie(double _r, double _kappa, double _Mm, std::vector<double> cp_coeff) : r{_r}, kappa{_kappa}, Mm{_Mm} 
    {
        a = cp_coeff[0];
        b = cp_coeff[1];
        c = cp_coeff[2];
        d = cp_coeff[3];
        e = cp_coeff[4];
        f = cp_coeff[5];

        a1 = a;
        b1 = b/2;
        c1 = c/3;
        d1 = d/4;
        e1 = e/5;
        f1 = f/6;
    };

    inline double cp(double T)
    {
        T = std::min(T,5000.0);
        // return a + b*T + c*pow(T,2) + d*pow(T,3) + e*pow(T,4) + f*pow(T,5);

        double result = f*T;
        result = (result + e)*T;
        result = (result + d)*T;
        result = (result + c)*T;
        result = (result + b)*T;
        result += a;
        return result;
    }

    inline double h(double T)
    {
        T = std::min(T,5000.0);
        // return a*T + b*pow(T,2)/2 + c*pow(T,3)/3 + d*pow(T,4)/4 + e*pow(T,5)/5 + f*pow(T,6)/6;

        double result = f1*T;
        result = (result + e1)*T;
        result = (result + d1)*T;
        result = (result + c1)*T;
        result = (result + b1)*T;
        result = (result + a1)*T;
        return result;
    }
};