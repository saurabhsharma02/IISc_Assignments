//************************************************//
// ME257 Homework 2
// Submitted by : Saurabh Sharma
// IDE : Visual Studio Code
//************************************************//

// Includes
#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include "Eigen/Sparse"   // Import Eigen library

//************************************************//
// User defined functions
//************************************************//

// Function to calculate the nodal coordinates
void  getCoordinates(double NodalCoordinates[], double h, int nNodes)
{
    for(int i =0; i<nNodes; i++)
        NodalCoordinates[i] = static_cast<double> (i) * h;
}

// Function to create the connectivity list
void getConnectivity(int Connectivity[], int nElements)
{
    for(int e = 0; e<nElements; e++)
    {
        Connectivity[2*e] = e;
        Connectivity[2*e + 1] = e+1;
    }
}

// Function to compute element mass matrix
void ComputeElementMassMatrix (double Xa, double Xb, double ElemMassMatrix[][2])
{
            ElemMassMatrix [0][0] = (Xb-Xa)/3.0;
            ElemMassMatrix [0][1] = (Xb-Xa)/6.0;
            ElemMassMatrix [1][0] = (Xb-Xa)/6.0;
            ElemMassMatrix [1][1] = (Xb-Xa)/3.0;
}

//Function to compute the global mass matrix
void ComputeGlobalMassMatrix (std::vector<Eigen::Triplet<double> > &matTriplets, Eigen::SparseMatrix<double, Eigen::RowMajor, int> &GlobalMassMatrix, double NodalCoordinates [], int Connectivity[], int nElements)
{
    double ElemMassMatrix [2][2];
    for (int nElm = 0; nElm<nElements; nElm++)
    {
        const int a = Connectivity[2*nElm];
        const int b = Connectivity[2*nElm + 1];
        const double Xa = NodalCoordinates[a];
        const double Xb = NodalCoordinates[b];

        ComputeElementMassMatrix(Xa, Xb, ElemMassMatrix);  // Function call to obstain element stiffness matrix

        matTriplets.push_back(Eigen::Triplet<double> (a,a,ElemMassMatrix[0][0]));
        matTriplets.push_back(Eigen::Triplet<double> (a,b,ElemMassMatrix[0][1]));
        matTriplets.push_back(Eigen::Triplet<double> (b,a,ElemMassMatrix[1][0]));
        matTriplets.push_back(Eigen::Triplet<double> (b,b,ElemMassMatrix[1][1]));
    }
    matTriplets.shrink_to_fit();
    GlobalMassMatrix.setFromTriplets(matTriplets.begin(), matTriplets.end());
}

// Function to compute element force vector
void ComputeElementForceVector(int FunctionIdentifier, double Xa, double Xb, double ElemForceVector[])
{
    // f(x) = sin(pi*x)
    if (FunctionIdentifier == 1)
    {
        ElemForceVector[0] = (M_PI*(Xa-Xb)*std::cos(M_PI*Xa) - std::sin(M_PI*Xa) + std::sin(M_PI*Xb))/(std::pow(M_PI,2)*(Xa-Xb));
        ElemForceVector[1] = (M_PI*(-Xa+Xb)*std::cos(M_PI*Xb) + std::sin(M_PI*Xa) - std::sin(M_PI*Xb))/(std::pow(M_PI,2)*(Xa-Xb));
    }
    // f(x) = x^2
    else if (FunctionIdentifier == 2)
    {
        ElemForceVector[0] = - (1.0/12.0) * (Xa - Xb) * (3*pow(Xa,2) + 2*Xa*Xb + pow(Xb,2));
        ElemForceVector[1] = - (1.0/12.0) * (Xa - Xb) * (pow(Xa,2) + 2*Xa*Xb + 3*pow(Xb,2));
    }
    // f(x) = e^x
    else if (FunctionIdentifier == 3)
    {
        ElemForceVector[0] = (-exp(Xb) + exp(Xa)*(1.0 - Xa + Xb))/(Xa-Xb);
        ElemForceVector[1] = (-exp(Xa) + exp(Xb)*(1.0 + Xa - Xb))/(Xa-Xb);
    }

    // f(x) = (sin(3x))^2
    else if (FunctionIdentifier == 4)
    {
        ElemForceVector[0] = (std::cos(6*Xa) - std::cos(6*Xb) + 6*(Xa-Xb)*(-3*Xa + 3*Xb + std::sin(6*Xa))) / (72.0*(Xa-Xb));
        ElemForceVector[1] = (-std::cos(6*Xa) + std::cos(6*Xb) - 6*(Xa-Xb)*(3*Xa - 3*Xb + std::sin(6*Xb))) / (72.0*(Xa-Xb));
    }

    // f(x) = cos(3*pi*x)
    else if (FunctionIdentifier == 5)
    {
        ElemForceVector[0] = (-std::cos(3 * M_PI * Xa)+std::cos(3 * M_PI * Xb) + 3.0*M_PI*(-Xa+Xb)*(std::sin(3*M_PI*Xa))) / (9.0*std::pow(M_PI,2)*(Xa-Xb));
        ElemForceVector[1] = (-std::cos(3 * M_PI * Xa)+std::cos(3 * M_PI * Xb) + 3.0*M_PI*(-Xa+Xb)*(std::sin(3*M_PI*Xa))) / (9.0*std::pow(M_PI,2)*(Xa-Xb));
    }
}

// Function to compute global mass matrix
void ComputeGlobalForceVector(int FunctionIdentifier, Eigen::VectorXd &ForceVec, int nElements, double NodalCoordinates[], int Connectivity [])
{
    double ElementForceVector [2];
    for (int nElm = 0; nElm<nElements; nElm++)
    {

        const int a = Connectivity[2*nElm];
        const int b = Connectivity[2*nElm + 1];
        const double Xa = NodalCoordinates[a];
        const double Xb = NodalCoordinates[b];

        ComputeElementForceVector(FunctionIdentifier, Xa,Xb, ElementForceVector);
        
        ForceVec.coeffRef(a) += ElementForceVector[0];
        ForceVec.coeffRef(b) += ElementForceVector[1];
    }
}

// Fucntion to compute the value of f(x) at any point in [0,1]
double getFunctionValueAtX(double X,int FunctionIdentifier)
{
    double fx;
    if (FunctionIdentifier == 1)
        fx = static_cast<double> (std::sin(M_PI*X));

    else if (FunctionIdentifier == 2)
        fx = static_cast<double> (std::pow(X,2));
    
    else if (FunctionIdentifier == 3)
        fx = static_cast<double> (exp(X));
        
    else if (FunctionIdentifier == 4)
        fx = static_cast<double> (std::pow(std::sin(X),2));
        
    else if (FunctionIdentifier == 5)
        fx = static_cast<double> (std::cos(3*M_PI*X));

    return fx;

}

// Function to compute the Max norm of the solution
double getMaxNorm(int nNodes, double NodalCoordinates [], Eigen::VectorXd DOFSolution, int FunctionIdentifier, Eigen::VectorXd &Error)
{
    double maxError = 0.0;
    double fx;
    double X;
    for (int i = 0 ; i<nNodes; i++)
    {
        X = NodalCoordinates[i];
        fx = getFunctionValueAtX(X,FunctionIdentifier);
        Error.coeffRef(i) = abs(fx - DOFSolution(i)); 

        if(maxError < Error.coeffRef(i))
            maxError = Error.coeffRef(i);
    }
    return maxError;
}

//************************************************//
// Main function
//************************************************//

int main()
{
    // variables declaration
    const int X0 = 0;                   // Start of the interval
    const int X1 = 1;                   // End of the interval
    const int nNodes = 10;             // Number of nodes
    const int nElements = nNodes - 1;   // Number of elements
    const double h = (X1 - X0) / static_cast<double> (nElements);   // Interval size
    double NodalCoordinates [nNodes];   // Array of nodal coodinates
    int Connectivity [2*nElements];     // Array have the connectivity indices
    double ElemMassMatrix[2][2];        // Elemental Mass matrix
    int FunctionIdentifier = 3;         // An indentifier to select the appropriate function

    //FunctionIdentifier=1 --> f(x) = sin(pi*x)
    //FunctionIdentifier=2 --> f(x) = x^2
    //FunctionIdentifier=3 --> f(x) = exp(x)
    //FunctionIdentifier=4 --> f(x) = sin^2(3*pi*x)
    //FunctionIdentifier=5 --> f(x) = cos(3*pi*x)
    
    //Declaration of Triplet, Sparse matrix and the force vector
    std::vector<Eigen::Triplet<double> > matTriplets;
    Eigen::SparseMatrix<double, Eigen::RowMajor, int> GlobalMassMatrix (nNodes, nNodes); // Mass matrix for the system
    Eigen::VectorXd GlobalForceVec(nNodes); // Global force vector

    // Solver creation and DOF vector creation
    Eigen::SparseLU<Eigen::SparseMatrix< double, Eigen::RowMajor> > Solver; 
    Eigen::VectorXd DOFSolution (nNodes); // vector to store the solution for the DOFs

    getCoordinates(NodalCoordinates, h, nNodes); // get the cooridnates of the nodes

    getConnectivity(Connectivity, nElements); // get the element connectivity table

    matTriplets.reserve(4*nElements);  // Triplets to store non-zero elements of the sparse matrix
    ComputeGlobalMassMatrix(matTriplets, GlobalMassMatrix, NodalCoordinates, Connectivity, nElements); // Assemble and get Global Mass matrix

    GlobalForceVec.setZero();
    ComputeGlobalForceVector(FunctionIdentifier,GlobalForceVec,nElements,NodalCoordinates,Connectivity); // Assemble and get Global Force vector

    // Solving the linear system of equations
    Solver.analyzePattern(GlobalMassMatrix);
    Solver.factorize(GlobalMassMatrix);
    DOFSolution = Solver.solve(GlobalForceVec);

    //Error calculation
    double MaxNorm = 0.0;
    Eigen::VectorXd Error(nNodes);
    Error.setZero();
    MaxNorm = getMaxNorm(nNodes, NodalCoordinates, DOFSolution, FunctionIdentifier, Error); // Fucntion call to calculate the max norm

    // Write nodal DOF to file
    std::fstream fout1;
    fout1.open((char*)"solution.dat", std::ios::out);
    for (int i= 0; i<nNodes; i++)
    fout1<< NodalCoordinates[i] << " " << DOFSolution(i) << "\n";
    fout1.close();

    //write error to file
    std::fstream fout2;
    fout2.open((char*)"error.dat", std::ios::out);
    for (int i= 0; i<nNodes; i++)
    fout2<< i+1 << " " << Error(i) << "\n";
    fout2.close();

    std::cout<<"Max Norm is "<<MaxNorm<<"\n";
    std::cout<<"Program terminated successfully";

    return 0;
}
//************************************************//
//                  End of Code
//************************************************//
