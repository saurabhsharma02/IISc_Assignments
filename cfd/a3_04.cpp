#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Define constants
const double PI = 3.14159265358979323846;
const double L = 1.0;
const double T = 1.0;

// Define functions for the exact solution and the source term
double exact(double x, double y, double t) {
    return sin(2 * PI * x) * sin(2 * PI * y) * sin(2 * PI * t);
}

double source(double x, double y, double t) {
    return sin(2 * PI * x) * sin(2 * PI * y) * sin(2 * PI * t);
}

// Define the matrix diagonalization method for solving the linear system of equations
void matrixDiagonalization(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d) {
    int n = d.size();

    // Forward sweep
    for (int i = 1; i < n; i++) {
        double m = a[i] / b[i - 1];
        b[i] = b[i] - m * c[i - 1];
        d[i] = d[i] - m * d[i - 1];
    }

    // Backward substitution
    d[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        d[i] = (d[i] - c[i] * d[i + 1]) / b[i];
    }
}

int main() {
    // Set up the parameters and mesh
    double dt = 0.0001;
    int N_min = 10, N_max = 100;
    vector<int> N_values;
    vector<double> L2_errors, Linf_errors;

    // Loop over different mesh sizes
    for (int N = N_min; N <= N_max; N += 10) {
        // Set up the mesh
        double dx = L / (N - 1);
        double dy = dx;
        double alpha = dt / (2 * dx * dx);

        int Nx = N, Ny = N;
        vector<vector<double>> u(Nx, vector<double>(Ny, 0.0));

        // Initialize the solution
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double x = i * dx;
                double y = j * dy;
                u[i][j] = exact(x, y, 0);
            }
        }

        // Apply the Crank-Nicolson method
        auto start = high_resolution_clock::now();
        for (double t = dt; t <= T; t += dt) {
            vector<double> a(Nx, -alpha);
            vector<double> b(Nx, 1 + 2 * alpha);
            vector<double> c(Nx, -alpha);
            vector<double> d(Nx, 0);

            // Set the boundary conditions
            for (int i = 0; i < Nx; i++) {
                u[i][0] = 0;
                u[i][Ny - 1] = 0;
            }
            for (int j =
                    0; j < Ny; j++) {
            u[0][j] = 0;
            u[Nx - 1][j] = 0;
        }

        // Set up the right-hand side vector
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                double x = i * dx;
                double y = j * dy;
                d[i] = u[i][j] + alpha * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4 * u[i][j]) + dt * source(x, y, t + dt / 2);
            }
        }

        // Solve the linear system of equations using the matrix diagonalization method
        matrixDiagonalization(a, b, c, d);

        // Update the solution
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                u[i][j] = d[i];
            }
        }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    // Compute the L2 and Linfinity errors
    double L2_error = 0, Linf_error = 0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double x = i * dx;
            double y = j * dy;
            double u_exact = exact(x, y, T);
            double u_approx = u[i][j];
                L2_error += pow(u_exact - u_approx, 2);
                Linf_error = max(Linf_error, abs(u_exact - u_approx));
            }
        }
        L2_error = sqrt(L2_error / (Nx * Ny));
        L2_errors.push_back(L2_error);
        Linf_errors.push_back(Linf_error);
        N_values.push_back(N);

        // Print the results
        cout << "N = " << N << ", L2 error = " << L2_error << ", Linfinity error = " << Linf_error << ", Time elapsed = " << duration.count() << " microseconds" << endl;
    }

    // Write the L2 and Linfinity errors to CSV files
    ofstream fout_L2("L2_errors.csv");
    ofstream fout_Linf("Linf_errors.csv");
    fout_L2 << "N,L2 error" << endl;
    fout_Linf << "N,Linfinity error" << endl;
    for (int i = 0; i < N_values.size(); i++) {
        fout_L2 << N_values[i] << "," << L2_errors[i] << endl;
        fout_Linf << N_values[i] << "," << Linf_errors[i] << endl;
    }
    fout_L2.close();
    fout_Linf.close();

    return 0;
    }
