#include "geometry.hpp"
#include <iostream>
#include <memory>
#include <cmath>

geometry::curve::curve(){}
std::vector<double> geometry::curve::interpolate(double t)
{
    return std::vector<double>{};
}

geometry::line::line(double _x0, double _y0, double _x1, double _y1)
{
    x0 = _x0;
    y0 = _y0;
    x1 = _x1;
    y1 = _y1;
}

std::vector<double> geometry::line::interpolate(double t)
{
    double x = (x1-x0)*t + x0;
    double y = (y1-y0)*t + y0;

    return std::vector<double>{x,y};
}

geometry::arc::arc(double _x0, double _y0, double _rx, double _ry, double _x1, double _y1)
{
    x0 = _x0;
    y0 = _y0;
    x1 = _x1;
    y1 = _y1;
    rx = _rx;
    ry = _ry;

    alfa = std::atan((y0-ry)/(x0-rx));
    phi = std::atan((y1-ry)/(x1-rx)) - alfa;
    r = std::sqrt(( std::pow(x0-rx,2) + std::pow(y0 - ry,2) ));

    std::cout << r << " " << alfa << " " << phi << "\n";
}

std::vector<double> geometry::arc::interpolate(double t)
{
    double x = r*std::cos(phi*t + alfa) + rx;
    double y = r*std::sin(phi*t + alfa) + ry;

    return std::vector<double>{x,y};
}

void geometry::test(std::vector<std::unique_ptr<curve>> const& vec)
{
    std::vector<double> coor(2,0.0);

    for(auto const& c : vec)
    {
        for(double t : {0.0,0.2,0.4,0.6,0.8,1.0})
        {
            coor = c->interpolate(t);

            std::cout << coor[0] << " " << coor[1] << "\n";
        }
    }
}