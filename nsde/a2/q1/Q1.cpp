#include<fstream>
#include<vector>
#include<cassert>
#include<iostream>
#include<Eigen/Sparse>

typedef Eigen::Triplet<double> tri;
typedef double Real;
using namespace std;

void readFile(vector<double> &values)
{
    ifstream fin;
    fin.open("./input.txt", ios::in);
    assert(!fin.fail());
    double temp = 0.0;
    while(fin >> temp)
    {
        values.push_back(temp);        
    }
    fin.close(); 
}




int main(){
    vector<double> inputs;
    readFile(inputs);
    int Nx = inputs[0];
    int Ny = inputs[1];
    Real x_start = inputs[2];
    Real x_end = inputs[3];
    Real y_start = inputs[4];
    Real y_end = inputs[5];
    double h = (x_end- x_start)/Nx;
    double k = (y_end- y_start)/Ny;
    double a = (h*h*(k*k+1))/(2*h*h+2*k*k-2*k*k*k*h*h);
    double b = ((k*k))/(2*h*h+2*k*k-2*k*k*k*h*h);
    double c = ((k*k)/(2*h*h+2*k*k));
    double d = ((h*h)/(2*h*h+2*k*k));

    Eigen::SparseMatrix <double,Eigen::RowMajor> A(((Nx-2)*Ny),((Nx-2)*Ny));
    std::vector<Eigen::Triplet<Real>> v ;
    A.setZero();
    A.setFromTriplets(v.begin(),v.end());

    //Bottom surface

    for (int i=0;i <Nx-2; i++)
    {
        v.push_back(tri(i,(i+Nx-2), (-h*h/(h*h+k*k))));
        v.push_back(tri(i,i,1));
        if(i!=0 && i!=Nx-3){
        v.push_back(tri(i,(i-1), (-h*h/(h*h+k*k))));   
        v.push_back(tri(i,(i+1), (-k*k/(h*h+k*k))));
        }
        else if (i==0) {
        v.push_back(tri(i,(i+1), (-k*k/(h*h+k*k))));
        }
        else 
        v.push_back(tri(i,(i-1), (-h*h/(h*h+k*k)))); 
    }

    //Top surface

    for (int i=((Nx-2)*(Ny-1)); i<((Nx-2)*(Ny)); i++ )
    {
        v.push_back(tri(i,i,1));
        v.push_back(tri(i,(i-Nx-2), -a));

        if(i!=((Nx-2)*(Ny-1)) && i!=((Nx-2)*(Ny)-1)){
        v.push_back(tri(i,(i-1), -b));   
        v.push_back(tri(i,(i+1), -b));
        }
        else if (i==(Nx-2)*(Ny-1)) {
        v.push_back(tri(i,(i+1), -b));
        }
        else 
        v.push_back(tri(i,(i-1), -b)); 
    }


    //Interior
    int k1 = (Nx-2);
    for (int j=0 ; j<(Ny-2); j++)
    {  
        for (int i=0; i<(Nx-2); i++)
        {
        v.push_back(tri(k1,k1,1));
        v.push_back(tri(k1,(k1-Nx+2), -d));
        v.push_back(tri(k1,(k1+Nx-2), -d));

        if(i!=0&& i!=(Nx-3)){
        v.push_back(tri(k1,(k1-1), -c));   
        v.push_back(tri(k1,(k1+1), -c));
        }
        else if (i==0) {
        v.push_back(tri(k1,(k1+1), -c));
        }
        else 
        v.push_back(tri(k1,(k1-1), -c)); 
        }

    }
    

    //RHS
    Eigen::VectorXd e((Nx-2)*(Ny));
    e.setZero();
    for (int i=0; i<Ny-1;i++){
        e(i*(Nx-2))= -30;
        e(i*(Nx-2)+(Nx-3))=-60;
    }
    e((Nx-2)*(Ny-1))=120*h*h/k-30;
    e((Nx-2)*Ny-1)=120*h*h/k-60;
    for(int i=(Nx-2)*(Ny-1)+1; i<(Nx-2)*Ny-1; i++){
        e(i)=120*h*h/k;
    }

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> Solver;
    Solver.analyzePattern(A);
    Solver.factorize(A);
	Eigen::VectorXd x((Nx-2)*(Ny));
	x = Solver.solve(e);
    std::vector<double> Solution(x.data(), x.data() + x.size());
    
    ofstream file("temperature_distribution.txt");
    file << x << endl;

    file.close();
    return 0;


 }



       
     







