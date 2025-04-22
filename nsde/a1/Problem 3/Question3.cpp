#include <iostream>
#include <fstream>
#include <cmath>
#include "assert.h"
#include <vector>

using namespace std;

// Function to calculate the derivative of v with respect to t
double f(double t, double x, double v, double c, double k, double m) {
    return (-c*v - k*x) / m;
}

void rungeKutta(double t0, double x0, double v0, double c, double k, double m, double h, int n, string filename) {
    double t = t0, x = x0, v = v0;
    ofstream file;
    file.open(filename);

    file << t << " " << x << endl;

    for (int i = 0; i < n; i++) {
        double k1x = h*v;
        double k1v = h*f(t, x, v, c, k, m);
        double k2x = h*(v + 0.5*k1v);
        double k2v = h*f(t + 0.5*h, x + 0.5*k1x, v + 0.5*k1v, c, k, m);
        double k3x = h*(v + 0.5*k2v);
        double k3v = h*f(t + 0.5*h, x + 0.5*k2x, v + 0.5*k2v, c, k, m);
        double k4x = h*(v + k3v);
        double k4v = h*f(t + h, x + k3x, v + k3v, c, k, m);

        x = x + (k1x + 2*k2x + 2*k3x + k4x) / 6;
        v = v + (k1v + 2*k2v + 2*k3v + k4v) / 6;
        t = t + h;

        file << t << " " << x << endl;
    }

    file.close();
}

void readFile(vector<double> &values)
{
    ifstream fin;
    fin.open("./input", ios::in);
    assert(!fin.fail());

    double temp = 0.0;
    
    while(fin >> temp)
    {
        values.push_back(temp);        
    }
    fin.close();   
}

int main() {
    vector<double> inputs;
    readFile(inputs);

    
    const double m = inputs[0]; // mass in kg
    const double k = inputs[1]; // spring constant in N/m
    const double x0 = inputs[2]; // initial displacement in m
    const double v0 = inputs[3]; // initial velocity in m/s
    const double t0 = inputs[4]; // initial time in s
    const double tmax = inputs[5]; // maximum time in s

    int n = 1500; // Number of steps
    double h = 15.0 / n; // Step size
    string filename = "";

    // Under-damped condition
    double c = 5.0;
    filename = "underdamped.txt";
    rungeKutta(t0, x0, v0, c, k, m, h, n, filename);

    // Critically damped condition
    c = 40.0;
    filename = "criticallydamped.txt";
    rungeKutta(t0, x0, v0, c, k, m, h, n, filename);

    // Over-damped condition
    c = 200.0;
    filename = "overdamped.txt";
    rungeKutta(t0, x0, v0, c, k, m, h, n, filename);

    return 0;
}
