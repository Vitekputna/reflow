#pragma once
#include <vector>
#include <memory>

namespace geometry
{

    class curve
    {
        public:
        double x0,y0;

        curve();
        virtual std::vector<double> interpolate(double t);
    };


    class line : public curve
    {
        public:
        double x1,y1;

        line(double _x0, double _y0, double _x1, double _y1);
        std::vector<double> interpolate(double t);
    };

    class arc : public curve
    {
        public:
        double x1,y1,rx,ry;

        double phi, alfa, r;

        arc(double _x0, double _y0, double _rx, double _ry, double _x1, double _y1);
        std::vector<double> interpolate(double t);
    };

    void test(std::vector<std::unique_ptr<curve>> const& vec);
};

typedef std::vector<std::unique_ptr<geometry::curve>> geometry_vector;