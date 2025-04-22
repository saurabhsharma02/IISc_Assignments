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

    
int main() {
    vector<double> inputs;
    readFile(inputs);

    // Initial conditions
    double x = inputs[0]; 
    double y = inputs[1];

    //Step size
    const double h = inputs[2]; 

    // Number of time steps
    const int N = inputs[3]; 
    double y_new;

    ofstream outfile("EulerQ4.csv");

    outfile << "x,y" << endl; // header row

    for (int i = 0; i <= N; i++) {
        outfile << x << "," << y << endl;
        y_new = (y + h*(2e5*exp(-x) - exp(-x)))/(1 + 2e5*h);
        y = y_new;
        x = x + h;
    }

    outfile.close();

    return 0;
}
