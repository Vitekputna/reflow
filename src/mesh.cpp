#include "mesh.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include "thermodynamics.hpp"

mesh::mesh(){}

mesh::mesh(int _N) : N{_N+2}
{
    // Resize fields
    x.resize(N);
    xf.resize(N-1);
    A.resize(N);
    Af.resize(N-1);

    // fill x field
    double step = (x_to-x_from)/(N-3);
    dx_min = step;

    for(int i = 1; i < N-1; i++)
    {
        x[i] = (i-1)*step + x_from;
    }

    x[0] = x[1] - (x[2] - x[1]);
    x[N-1] = x[N-2] + (x[N-2] - x[N-3]);

    // fill xf field

    for(int i = 0; i < N-1; i++)
    {
        xf[i] = (x[i+1] + x[i])/2;
    }

    // fill A field
    step = (A_to-A_from)/(N-3);

    for(int i = 1; i < N-1; i++)
    {
        A[i] = i*step + A_from;
    }

    A[0] = A[1];
    A[N-1] = A[N-2];

    // fill Af field

    for(int i = 0; i < N-1; i++)
    {
        Af[i] = (A[i+1] + A[i])/2;
    }
}

mesh::mesh(int _N, double from, double to) : N{_N+2}, x_from{from}, x_to{to}
{
    // Resize fields
    x.resize(N);
    xf.resize(N-1);
    A.resize(N);
    Af.resize(N-1);

    // fill x field
    double step = (x_to-x_from)/(N-3);
    dx_min = step;

    for(int i = 1; i < N-1; i++)
    {
        x[i] = (i-1)*step + x_from;
    }

    x[0] = x[1] - (x[2] - x[1]);
    x[N-1] = x[N-2] + (x[N-2] - x[N-3]);

    // fill xf field

    for(int i = 0; i < N-1; i++)
    {
        xf[i] = (x[i+1] + x[i])/2;
    }

    // fill A field
    step = (A_to-A_from)/(N-3);

    for(int i = 1; i < N-1; i++)
    {
        A[i] = i*step + A_from;
    }

    A[0] = A[1];
    A[N-1] = A[N-2];

    // fill Af field

    for(int i = 0; i < N-1; i++)
    {
        Af[i] = (A[i+1] + A[i])/2;
    }
}

void mesh::bump()
{
    for(int i = 0; i < N; i++)
    {
        if(x[i] > 0.3 && x[i] < 0.7)
        {
            A[i] -= 0.2*sin((x[i] - 0.3)/0.4*3.14159);
        }
    }

    for(int i = 0; i < N-1; i++)
    {
        Af[i] = (A[i+1] + A[i])/2;
    }
}

void mesh::cubic(std::vector<std::vector<std::vector<double>>> curves, int n)
{
    int n_curves = curves.size();

    double px, py;
    double t;

    std::vector<double> inter_x;
    std::vector<double> inter_y;

    for(int i = 0; i < n_curves; i++)
    {
        for(int j = 0; j < n; j++)
        {
            t = j*1/(double(n-1));
            
            px = (2*t*t*t - 3*t*t + 1)*curves[i][0][0] + (t*t*t -2*t*t + t)*curves[i][2][0]
               + (-2*t*t*t + 3*t*t)*curves[i][1][0] + (t*t*t -t*t)*curves[i][3][0];

            py = (2*t*t*t - 3*t*t + 1)*curves[i][0][1] + (t*t*t -2*t*t + t)*curves[i][2][1]
               + (-2*t*t*t + 3*t*t)*curves[i][1][1] + (t*t*t -t*t)*curves[i][3][1];
            
            inter_x.push_back(px);
            inter_y.push_back(py);
        }
    }

    A[0] = inter_y[0];

    for(int i = 1; i < N-1; i++)
    {
        for(unsigned int j = 0; j < inter_x.size() - 1; j++)
        {
            if(x[i] >= inter_x[j] && x[i] <= inter_x[j+1])
            {
                A[i] = (inter_y[j+1] - inter_y[j])/(inter_x[j+1] - inter_x[j])*(x[i] - inter_x[j]) + inter_y[j];
                break;
            }
        }
    }

    A[N-1] = inter_y.back();

    // fill Af field

    for(int i = 0; i < N-1; i++)
    {
        Af[i] = (A[i+1] + A[i])/2;
    }
}

void mesh::export_to_file()
{
    auto stream =  std::ofstream("out/A.txt");

    for(int i = 0; i < N; i++)
    {
        stream << A[i] << "\n";
    }

    stream << "\n";
    stream.close();
}