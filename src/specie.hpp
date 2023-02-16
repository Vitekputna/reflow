#pragma once

struct specie
{
    double r;
    double kappa;
    double Mm;
    double a,b,c,d,e,f,g;

    specie(double _r, double _kappa, double _Mm, std::vector<double> cp_coeff) : r{_r}, kappa{_kappa}, Mm{_Mm} 
    {
        a = cp_coeff[0];
        b = cp_coeff[1];
        c = cp_coeff[2];
        d = cp_coeff[3];
        e = cp_coeff[4];
        f = cp_coeff[5];
        g = cp_coeff[6];
    };

    double cp(double T)
    {
        return a/T + b + c*T + d*T*T + e*T*T*T + f*T*T*T*T + g*T*T*T*T*T;
    }
};