#pragma once

struct specie
{
    double r;
    double kappa;
    double Mm;

    specie(double _r, double _kappa, double _Mm) : r{_r}, kappa{_kappa}, Mm{_Mm} {};
};