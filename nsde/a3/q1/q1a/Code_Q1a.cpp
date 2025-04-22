#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include "assert.h"

// default name space
using namespace std;

// necessary declarations
typedef double real;

// allocation
void allocate1D(real *&x, int n)
{
    x = new real[n];
}

// printing output in a file
void print_array(real x[], int n, const char *filename)
{
    cout << " ...printing " << filename << endl;
    string name = filename;
    ofstream fout;
    fout.open(name);
    for (int i = 0; i < n; i++)
        fout << x[i] << scientific << endl;
    fout.close();
}

void readFile(vector<double> &values)
{
    ifstream fin;
    fin.open("./input", ios::in);
    assert(!fin.fail());

    double temp = 0.0;

    while (fin >> temp)
    {
        values.push_back(temp);
    }
    fin.close();
}

// begin the main program
int main()
{

    vector<double> inputs;
    readFile(inputs);
    // variables
    int N;
    real dx, dt, r, L, pi = 2 * acos(0.0), tend, Uavg = 0.0;
    real *U, *Unew, *Uana, *Uerr;

    cout << "Solving an ODE: " << endl;

    N = inputs[0];
    L = 1.0;
    tend = 0.075;
    dx = L / N;
    dt = 0.0004;
    r = dt / dx;
    N = N + 1;

    // allocate arrays
    allocate1D(U, N);
    allocate1D(Unew, N);
    allocate1D(Uana, N);
    allocate1D(Uerr, N);

    // initial condition
    for (int i = 0; i < N; i++)
    {
        U[i] = sin(4 * pi * i * dx) + sin(6 * pi * i * dx) + sin(10 * pi * i * dx);
    }

    // calculation of U
    for (double k = dt; k <= tend; k += dt)
    {
        // left corner
        if (U[0] >= 0)
            Unew[0] = U[0] - r * U[0] * (U[1] - U[0]);
        else
            Unew[0] = U[0] - r * U[0] * (U[0] - U[N - 2]);

        // interior points
        for (int i = 1; i < N - 1; i++)
            if (U[i] < 0)
                Unew[i] = U[i] - r * U[i] * (U[i + 1] - U[i]);
            else
                Unew[i] = U[i] - r * U[i] * (U[i] - U[i - 1]);

        // right corner
        if (U[N - 1] >= 0)
            Unew[N - 1] = U[N - 1] - r * U[N - 1] * (U[N - 1] - U[N - 2]);
        else
            Unew[N - 1] = U[N - 1] - r * U[N - 1] * (U[0] - U[N - 1]);

        U = Unew;
    }

    // error calculation
    for (int i = 0; i < N; i++)
    {
        Uana[i] = exp(-4 * tend) * sin(2 * (i - 1) * dx);
        Uerr[i] = abs(Uana[i] - U[i]);
        Uavg = Uavg + Uerr[i];
    }

    cout << "average error = " << Uavg / (N) << endl;
    cout << "solved" << endl;

    // write data
    print_array(U, N, "U_64.txt");

    return 0;
}
