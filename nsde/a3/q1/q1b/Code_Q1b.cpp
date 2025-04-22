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
    real dx, dt, L, alpha, pi = 2 * acos(0.0), tend, Uavg = 0.0;
    real *U, *Unew, *Uana, *Uerr;

    cout << "Solving an ODE: " << endl;

    N = inputs[0];
    L = 1.0;
    tend = 0.075;
    dx = L / N;
    dt = 0.0004;
    alpha = 0.001;
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
        Unew[0] = U[0] - (dt / (2 * dx)) * U[0] * (U[1] - U[N - 2]) + ((alpha * dt) / (dx * dx)) * (U[1] - 2 * U[0] + U[N - 2]);

        // interior points
        for (int i = 1; i < N - 1; i++)
            Unew[i] = U[i] - (dt / (2 * dx)) * U[i] * (U[i + 1] - U[i - 1]) + ((alpha * dt) / (dx * dx)) * (U[i + 1] - 2 * U[i] + U[i - 1]);

        // right corner
        Unew[N - 1] = U[N - 1] - (dt / (2 * dx)) * U[N - 1] * (U[0] - U[N - 2]) + ((alpha * dt) / (dx * dx)) * (U[0] - 2 * U[N - 1] + U[N - 2]);

        U = Unew;
    }

    cout << "solved" << endl;

    // write data
    print_array(U, N, "U_1024.txt");

    return 0;
}
