#include "mesh.hpp"

#include <iostream>

mesh::mesh(){}

mesh::mesh(int _N) : N{_N+2}
{
    // Resize fields
    x.resize(N);
    xf.resize(N-1);
    A.resize(N);

    // fill x field
    double step = (x_to-x_from)/(N-3);

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
    step = (A_to-A_from)/(N-1);

    for(int i = 1; i < N-1; i++)
    {
        A[i] = i*step + A_from;
    }

    A[0] = A[1];
    A[N-1] = A[N-2];

}