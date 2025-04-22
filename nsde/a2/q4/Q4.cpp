#include <iostream>
#include <fstream>
#include <cmath>
#include "assert.h"
#include <vector>

using namespace std;

const double pi = 2 * acos(0.0);

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

double analytical_solution(double x, double t)
{
    return (sin(4 * x) * exp(-24 * t) + sin(x) * exp(-1.5 * t));
}

double compute_error(double *u, double dx, int N, double t)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        double x = i * dx;
        double error = fabs(u[i] - analytical_solution(x, t));
        sum += error;
    }
    double avg_error = sum / N;
    cout << "Average absolute error at t = " << t << " : " << avg_error << endl;
    return avg_error;
}

int main()
{
    vector<double> inputs;
    readFile(inputs);
    const double rd_array[2] = {inputs[0], inputs[1]};
    const int N_array[4] = {int(inputs[2]), int(inputs[3]), int(inputs[4]), int(inputs[5])};
    // Creating .dat file for the error computed
    ofstream file("solution.dat");
    for (int r = 0; r < 2; r++)
    {
        for (int n = 0; n < 4; n++)
        {

            double rd = rd_array[r];
            int N = N_array[n];

            double dx = 2 * pi / N;
            double dt = (rd * dx * dx) / 1.5;
            int K = (int)(0.4 / dt);
            dt = 0.4 / K;
            const double tend = 0.4;

            // Forming arrays
            double u[N], un[N], x[N];
            double error[N];

            // Setting initial condition as given in the problem
            for (int i = 0; i < N; i++)
            {
                x[i] = 2 * pi * i / N;
                u[i] = sin(4 * x[i]) + sin(x[i]); // Initial condition for the PDE at time t=0 and spatial location x[i].
                un[i] = u[i];                     // un[i] is the value of the numerical solution at time t + dt and spatial location x[i]. Sets the value of un[i] to be equal to u[i]. This initializes the numerical solution at the same value as the initial condition.
                error[i] = 0;                     // Initializing an array that will be used to store the error between the numerical solution and the true solution.
            }

            // Evolving the solution in time
            double t;
            for (t = 0; t < tend; t += dt)
            {
                // Updating the solution at the interior points
                for (int i = 1; i < N-1; i++)
                {
                    un[i] = u[i] + rd * (u[i + 1] - 2 * u[i] + u[i - 1]); // Updating the solution at each interior grid point.
                }
                // Updating the solution at the boundary points (periodic boundary condition)
                un[0] = u[0] + rd * (u[1] - 2 * u[0] + u[N - 2]);
                un[N-1] = un[0];
                // Updating the solution for the next time step
                for (int i = 0; i < N; i++)
                {
                    u[i] = un[i];
                }
            }

            // Computing the error at t = 0.4
            double E = compute_error(u, dx, N, 0.4);

            // Write to file
            file << rd << " " << N << " " << E << endl;
        }
    }

    file.close();
    return 0;
}
