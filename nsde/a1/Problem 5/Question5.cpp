#include <iostream>
#include <fstream>
#include <cmath>
#include "assert.h"
#include <vector>

using namespace std;
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

int main()
{
    vector<double> inputs;
    readFile(inputs);


    const int N = 1024;
    const double h = 1.0 / (N - 1);
    const double* S = &inputs[0];
    double T[N] = {0};
    double a[N] = {0};
    double b[N] = {0};
    double c[N] = {0};
    long double d[N] = {0};
    double r[N] = {0};
    double s[N] = {0};
    double t[N] = {0};

    
  
    for (int j = 0; j < inputs.size(); j++)
    {   
        d[1]= 4.0 / (3.0* h * h);
        a[1]= -4.0 / (3.0* h * h);
        c[1]= S[j];

        // set up coefficients for Thomas Algorithm
        for (int i = 2; i < N - 2; i++)
        {
            b[i] = - 1.0 / (h * h) + 1.0 / (2.0 * h * i*h);
            d[i] = +2.0 / (h * h) ;
            a[i] = -1.0 / (h * h) - 1.0 / (2.0 * h * i*h);
            c[i] = S[j];
        }
    b[N-2] = 1.0 / (2.0 * (N - 2) * h * h) -1.0 / (h * h);
    d[N-2] = 2.0 / (h * h);
    c[N-2] = S[j] + 1.0 / (2.0 * (N - 2) * h * h) + 1.0 / (h * h);

    // Thomas
    for (int i = 2; i < N - 1; i++)
    {   
        d[i] = d[i] - b[i] * a[i-1] / d[i-1];
        c[i]= c[i] - c[i-1] * b[i] / d[i-1];
    }

    T[N-1] = 1;
    T[N-2]= c[N-2] / d[N-2];
    for (int i = N - 3; i > 0; i--)
    {
        T[i] = (c[i] - a[i] * T[i+1]) / d[i];
    }
    T[0] = (4 * T[1] - T[2]) / 3;
    

        // write temperature values to csv file
        string filename = "Source_value_" + to_string(int(S[j])) + ".csv";
        ofstream outfile(filename);
        for (int i = 0; i < N; i++)
        {
            outfile << i * h << "," << T[i] << endl;
        }
        outfile.close();
    }

    return 0;
}
