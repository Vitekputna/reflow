#pragma once
#include <cmath>
#include <vector>
struct specie
{
    double r;
    double kappa;
    double Mm;
    double h_vap;

    double a,b,c,d,e,f;         //cp coeffs
    double a1,b1,c1,d1,e1,f1;   //cp coeffs divided

    double a2,b2,c2,d2,e2,f2;   //k coeffs
    double a3,b3,c3,d3,e3,f3;   //mu coeffs

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

    specie(double _r, double _kappa, double _Mm, std::vector<double> cp_coeff, std::vector<double> k_coeff) : specie(_r,_kappa,_Mm,cp_coeff)
    {
        a2 = k_coeff[0];
        b2 = k_coeff[1];
        c2 = k_coeff[2];
        d2 = k_coeff[3];
        e2 = k_coeff[4];
        f2 = k_coeff[5];
    };

    specie(double _r, double _kappa, double _Mm, std::vector<double> cp_coeff, std::vector<double> k_coeff, std::vector<double> mu_coeff) : specie(_r,_kappa,_Mm,cp_coeff,k_coeff)
    {
        a3 = mu_coeff[0];
        b3 = mu_coeff[1];
        c3 = mu_coeff[2];
        d3 = mu_coeff[3];
        e3 = mu_coeff[4];
        f3 = mu_coeff[5];
    };

    inline double cp(double T)
    {
        T = std::min(T,5000.0);

        double result = f*T;
        result = (result + e)*T;
        result = (result + d)*T;
        result = (result + c)*T;
        result = (result + b)*T;
        result += a;
        return result;
    }

    inline double k(double T)
    {
        T = std::min(T,5000.0);

        double result = f2*T;
        result = (result + e2)*T;
        result = (result + d2)*T;
        result = (result + c2)*T;
        result = (result + b2)*T;
        result += a2;
        return result;
    }

    inline double h(double T)
    {
        T = std::min(T,5000.0);

        double result = f1*T;
        result = (result + e1)*T;
        result = (result + d1)*T;
        result = (result + c1)*T;
        result = (result + b1)*T;
        result = (result + a1)*T;
        return result;
    }

    inline double mu(double T)
    {
        T = std::min(T,5000.0);

        double result = f3*T;
        result = (result + e3)*T;
        result = (result + d3)*T;
        result = (result + c3)*T;
        result = (result + b3)*T;
        result += a3;
        return result;
    }
};