#include <iostream>
#include <fstream>
#include <cmath>
#include "assert.h"
#include <vector>

using namespace std;
using Real = double;

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

int main()
{

    vector<double> inputs;
    readFile(inputs);
    const Real rd_array[2] = {inputs[0], inputs[1]};
    const int N_array[4] = {int(inputs[2]), int(inputs[3]), int(inputs[4]), int(inputs[5])};

    // Creating .dat file for the error computed
    ofstream file("solution.dat");

    for (int r = 0; r < 2; r++)
    {
        for (int n = 0; n < 4; n++)
        {
            double rd = rd_array[r];
            int N = N_array[n];
            const Real L = 2 * M_PI;
            const Real tend = 0.4;
            const Real dx = L / N;        // Spatial grid size
            const Real dt = rd * dx * dx; // Time step size

            // Declare and allocate memory for arrays
            Real u[N];
            Real u_new[N];
            Real u_err[N];
            Real u_ana[N];
            Real u_tol[N];

            // Initializing the values of the array u, which has size N, using a sine function
            for (int i = 0; i < N; i++)
            {
                u[i] = sin(4.0 * dx * i) + sin(dx * i);
            }

            // Apply periodic boundary conditions
            u[0] = u[N - 2];
            u[N - 1] = u[1];

            Real tolerance = 1.0;
            while (tolerance > 0.0001)
            {
                {
                    // Jacobi method to solve equation
                    for (int i = 1; i < N - 1; i++)
                    u_new[i] = (u[i] + rd * u_new[i - 1] + rd * u_new[i + 1]) / (1 + 2.0 * rd);
                    u_new[0] = (u[0] + rd * u_new[1] + rd * u_new[N - 2]) / (1 + 2.0 * rd);
                    u_new[N - 1] = (u[N - 1] + rd * u_new[N - 2] + rd * u_new[1]) / (1 + 2.0 * rd);

                    // calculate tolerance
                    for (int i = 0; i < N; i++)
                        u_tol[i] = std::abs(u_new[i] - u[i]);

                    tolerance = u_tol[0];
                    for (int i = 1; i < N; i++)
                        if (u_tol[i] > tolerance)
                            tolerance = u_tol[i];

                    // update u
                    for (int i = 0; i < N; i++)
                        u[i] = u_new[i];
                }

                {

                }
            }
          

            // calculate error
            Real u_avg = 0.0;
            for (int i = 0; i < N; i++)
            {
                u_ana[i] = (sin(4 * (i - 1) * dx) * exp(-24 * tend) + sin((i - 1) * dx) * exp(-1.5 * tend));
                u_err[i] = fabs(u_ana[i] - u[i]);
                u_avg += u_err[i];
            }
            u_avg = u_avg / N;

            std::cout << "average error = " << u_avg << '\n';
            std::cout << "solved\n";

            file << rd << " " << N << " " << u_avg << endl;
        }
    }

    return 0;
}
