//************************************************//
// ME257 Homework 5
// Submitted by : Saurabh Sharma
// IDE : Visual Studio Code
//************************************************//

// Includes
#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include "Eigen/Dense"                          // Import Eigen library
#include "Eigen/Sparse"                         // Import Eigen library

//************************************************//
// Structure data tyoe for the mesh
//************************************************//
struct Mesh
{
    int nNodes;
    int nElements;
    int nBdNodes;
    std::vector<std::vector<double> > Coordinates;
    std::vector<std::vector<int> > Connectivity;
    std::vector<int> bdNodes;
};

//************************************************//
// User defined functions
//************************************************//

//************************************************//
// Function to read the mesh data
void ReadMesh(std::string coordFileName, std::string connFileName, std::string bdFileName, Mesh & M)
{
    // Reading the nodal coordinate data
    std::ifstream FileOpen;
    FileOpen.open(coordFileName);

    if(!FileOpen.is_open())
        {
            std::cout<<"File not found \n Terminating the program";
            exit(0);
        }

    FileOpen >> M.nNodes;                               // reading the number of nodes
    double nodeCoordinates;
        for(int i = 0; i<M.nNodes; i++)
        {
            std::vector<double> temp;
            for(int j=0; j<2;j++)
            {  \
                FileOpen >> nodeCoordinates;
                temp.push_back(nodeCoordinates);
            }   
            M.Coordinates.push_back(temp);              // storing the coordinates
        }
    FileOpen.close();

    // Reading the connectivity data
    FileOpen.open(connFileName);

    if(!FileOpen.is_open())
        {
            std::cout<<"File not found \n Terminating the program";
            exit(0);
        }
    FileOpen >> M.nElements;                            // reading the number of Elements
    int Connectivity;
        for(int i = 0; i<M.nElements; i++)
        {
            std::vector<int> temp;
            for(int j=0; j<3;j++)
            {  
                FileOpen >> Connectivity;
                temp.push_back(Connectivity);
            }   
            M.Connectivity.push_back(temp);             // storing the element connectivity
        }
    FileOpen.close();

    // Reading the boundary nodes
    FileOpen.open(bdFileName);

    if(!FileOpen.is_open())
        {
            std::cout<<"File not found \n Terminating the program";
            exit(0);
        }
    FileOpen >> M.nBdNodes;                             // reading the number of boundary nodes
    int BoundaryNodes;
        for(int i = 0; i<M.nBdNodes; i++)
        {
            FileOpen >> BoundaryNodes;
            M.bdNodes.push_back(BoundaryNodes);         // storing the boundary nodes
        }
    FileOpen.close();                                   
    std::cout<<"\nAll files have been successfully read";
}
//************************************************//

//************************************************//
// Function to define the local to global map
int Local2GlobalMap(const Mesh &M, int elem_num, int loc_node_num)
{
    if (elem_num > M.nElements+1 || loc_node_num > 4)
       { 
        std::cout<<"Index out of bound\n"<<"Program is being terminated";
        exit(0);
        }
    else
         return (M.Connectivity[elem_num-1][loc_node_num-1]);
}

//************************************************//
// Function to get the basis fucntion at specified coordinates
double getBasisFunction(double Xi, double Eta, int i)
{
    double val, N1,N2,N3;
    N1 = Xi;
    N2 = Eta;
    N3 = 1-Xi-Eta;                                      // Partition of unity check

    if (i==1)
        val= N1;
    
    else if (i==2)
        val= N2;
    
    else if (i==3)
        val= N3;
    
    return val;
}

//************************************************//
// Function to define the local to global map
void ComputeElementStiffnessMatrix(const Mesh &M, int elemNum, Eigen::Matrix<double,3,3> &Ke)
{
    double Xi,Eta;
    Xi = 1.0;
    Eta = 0.0;
    int a,b,c;
    double x1,y1,x2,y2,x3,y3;
    Eigen::Matrix<double, 3, 2> DerivativeOfNhats;
    Eigen::Vector2d NiHatPrime,NjHatPrime;
    DerivativeOfNhats << 1.0,0.0,0.0,1.0,-1.0,-1.0; // Initialization of the derivative of the basis functions

    Eigen::Matrix2d gradPhi;
    double Jacobian;
    double N1_hat, N2_hat, N3_hat;

    //Getting the global number of node with respected to local
    a = Local2GlobalMap(M, elemNum,1);
    b = Local2GlobalMap(M, elemNum,2);
    c = Local2GlobalMap(M, elemNum,3);

    // Getting the coordinates of nodes
    x1 = M.Coordinates[a][0];
    y1 = M.Coordinates[a][1];
    x2 = M.Coordinates[b][0];
    y2 = M.Coordinates[b][1];
    x3 = M.Coordinates[c][0];
    y3 = M.Coordinates[c][1];

    // Computing grad Phi
    gradPhi(0,0) = x1-x3;
    gradPhi(0,1) = x2-x3;
    gradPhi(1,0) = y1-y3;
    gradPhi(1,1) = y2-y3;

    Eigen::Matrix2d gradPhiInvTranspose = (gradPhi.transpose()).inverse(); // Finding transpose of grad phi and inverting it
    Jacobian = gradPhi.determinant();                                      // Jacobian of the transformation

    // Check for the Jacobian value
    if (Jacobian <=0)
        {
            std::cout<<"Linear Mapping failed : Jacobian is negative\n"<<"Terminating the program";
            exit(0);
        }

    // Estimating and verifying the value of the computed Jacobian
    double PhysicalElementArea, ParametricElementArea,AreaRatio;
    Eigen::Matrix<double, 3, 3> CoodinateMatrix;
    CoodinateMatrix << x1,x2,x3,y1,y2,y3,1.0,1.0,1.0;
    
    PhysicalElementArea = CoodinateMatrix.determinant() * 0.5;                  // Determinant gives the area of the triangle
    ParametricElementArea = 0.5;

    AreaRatio = PhysicalElementArea/ParametricElementArea;

    // Checking the correctness of the Jacobian value
    if (Jacobian - AreaRatio > 1.0e-12)
    {
        std::cout<<"Warning: Jacobian value not accurate: Check linear mapping Phi and derivative of basis functions";
    }

    // Computing the element stiffness matrix
    for (int i = 0; i<3;i++)
    {
        NiHatPrime (0) = DerivativeOfNhats(i,0);        // Assigning the values for NiPrime
        NiHatPrime (1) = DerivativeOfNhats(i,1);        // Assigning the values for NiPrime
        
        for (int j = 0; j<3; j++)
        {   
            NjHatPrime (0) = DerivativeOfNhats(j,0);    // Assigning the values for NjPrime
            NjHatPrime (1) = DerivativeOfNhats(j,1);    // Assigning the values for NjPrime
            Ke(i,j) = ((gradPhiInvTranspose*NiHatPrime).dot(gradPhiInvTranspose*NjHatPrime)) * Jacobian;        // Populating the element stiffness matrix
        }
    }
}

//************************************************//
// Function to assemble the global stiffness matrix
void ComputeGlobalStiffnessMatrix(std::vector<Eigen::Triplet<double> > &Triplets, Eigen::SparseMatrix<double, Eigen::RowMajor, int> &GlobalStiffnessMatrix, const Mesh &M)
{
   int a,b,c;
   Eigen::Matrix<double,3,3> Ke;                                         // Element stiffness matrix

   for (int i = 1; i<M.nElements+1; i++)
   { 
        a = Local2GlobalMap(M, i,1);
        b = Local2GlobalMap(M, i,2);
        c = Local2GlobalMap(M, i,3);

        Ke.setZero();
        ComputeElementStiffnessMatrix(M,i,Ke);                           // Function call to obtain the element stiffness matrix

        // Generating the triplets for the sparse matrix
        Triplets.push_back(Eigen::Triplet<double> (a,a,Ke(0,0)));
        Triplets.push_back(Eigen::Triplet<double> (a,b,Ke(0,1)));
        Triplets.push_back(Eigen::Triplet<double> (a,c,Ke(0,2)));

        Triplets.push_back(Eigen::Triplet<double> (b,a,Ke(1,0)));
        Triplets.push_back(Eigen::Triplet<double> (b,b,Ke(1,1)));
        Triplets.push_back(Eigen::Triplet<double> (b,c,Ke(1,2)));

        Triplets.push_back(Eigen::Triplet<double> (c,a,Ke(2,0)));
        Triplets.push_back(Eigen::Triplet<double> (c,b,Ke(2,1)));
        Triplets.push_back(Eigen::Triplet<double> (c,c,Ke(2,2)));
   }

    Triplets.shrink_to_fit();
    GlobalStiffnessMatrix.setFromTriplets(Triplets.begin(), Triplets.end());    // Constructing the stiffness matrix from triplets

    std::cout<<"\nGlobal stiffness matrix computed successfully";
}

//************************************************//
// Function to compute the element force vector
void ComputeElementForceVector(int FunctionIdentifier, Eigen::Vector3d &Fe, const Mesh &M, int elemNum)
{
    double Xi,Eta;
    int a,b,c;
    double x1,y1,x2,y2,x3,y3;
    Eigen::Matrix2d gradPhi;
    double Jacobian;
    double N1_hat, N2_hat, N3_hat;
    double phiX, phiY;

    //Getting the global number of node with respected to local
    a = Local2GlobalMap(M, elemNum,1);
    b = Local2GlobalMap(M, elemNum,2);
    c = Local2GlobalMap(M, elemNum,3);

    // Getting the coordinates of nodes
    x1 = M.Coordinates[a][0];
    y1 = M.Coordinates[a][1];
    x2 = M.Coordinates[b][0];
    y2 = M.Coordinates[b][1];
    x3 = M.Coordinates[c][0];
    y3 = M.Coordinates[c][1];

    // Computing grad Phi
    gradPhi(0,0) = x1-x3;
    gradPhi(0,1) = x2-x3;
    gradPhi(1,0) = y1-y3;
    gradPhi(1,1) = y2-y3;
    Jacobian = gradPhi.determinant();                                      // Jacobian of the transformation

    // Gauss quadrature points for a three point rule on a parametric triangle
    double Xi_Gauss [] = {1.0/6.0 , 2.0/3.0, 1.0/6.0};
    double Eta_Gauss [] = {1.0/6.0, 1.0/6.0, 2.0/3.0};
    double weights_Gauss [] = {1.0/6.0, 1.0/6.0, 1.0/6.0};

    // Element force vector for forcing function f = 0
    if (FunctionIdentifier == 0)
        {
            Fe(0) = 0.0;
            Fe(1) = 0.0;
            Fe(2) = 0.0;
        }
    // Element force vector for forcing function f = 2
    else if (FunctionIdentifier == 1)
        {
            for (int i = 0; i<3; i++)
            {
                for (int j=0; j<3; j++)
                {
                    phiX = Xi_Gauss[j]*(x1-x3) + Eta_Gauss[j]*(x2-x3) + x3;
                    phiY = Xi_Gauss[j]*(y1-y3) + Eta_Gauss[j]*(y2-y3) + y3;
                    Fe(i)+= 2.0*weights_Gauss[j]*getBasisFunction(Xi_Gauss[j], Eta_Gauss[j],(i+1));        // Integral evaluated using 3 point guass quadrature
                }   
            } 
                Fe(0) = 2.0*Fe(0)*Jacobian;
                Fe(1) = 2.0*Fe(1)*Jacobian;      
                Fe(2) = 2.0*Fe(2)*Jacobian;              
        }

    // Element force vector for forcing function f = sin(3*pi*x)
    else if (FunctionIdentifier == 2)
        {
            for (int i = 0; i<3; i++)
            {
                for (int j=0; j<3; j++)
                {
                    phiX = Xi_Gauss[j]*(x1-x3) + Eta_Gauss[j]*(x2-x3) + x3;
                    phiY = Xi_Gauss[j]*(y1-y3) + Eta_Gauss[j]*(y2-y3) + y3;
                    Fe(i)+=weights_Gauss[j]*(std::sin(3*M_PI*phiX)*getBasisFunction(Xi_Gauss[j], Eta_Gauss[j],(i+1)));
                }  
            }  
                Fe(0) = 2.0*Fe(0)*Jacobian;
                Fe(1) = 2.0*Fe(1)*Jacobian;      
                Fe(2) = 2.0*Fe(2)*Jacobian;              
        }

    // Element force vector for forcing function f = sin(3*pi*x) + sin(3*pi*y)
    else if (FunctionIdentifier == 3)
        {
            for (int i = 0; i<3; i++)
            {
                for (int j=0; j<3; j++)
                {
                    phiX = Xi_Gauss[j]*(x1-x3) + Eta_Gauss[j]*(x2-x3) + x3;
                    phiY = Xi_Gauss[j]*(y1-y3) + Eta_Gauss[j]*(y2-y3) + y3;
                    Fe(i)+= weights_Gauss[j]*(std::sin(3*M_PI*phiX)*getBasisFunction(Xi_Gauss[j], Eta_Gauss[j],(i+1))) + weights_Gauss[j]*(std::sin(3*M_PI*phiY)*getBasisFunction(Xi_Gauss[j], Eta_Gauss[j],(i+1)));
                }  
            }  
                Fe(0) = Fe(0)*Jacobian;
                Fe(1) = Fe(1)*Jacobian;      
                Fe(2) = Fe(2)*Jacobian;              
        }

    // Element force vector for forcing function f = cos(7*pi*x) + cos(7*pi*y)
    else if (FunctionIdentifier == 4)
        {
            for (int i = 0; i<3; i++)
            {
                for (int j=0; j<3; j++)
                {
                    phiX = Xi_Gauss[j]*(x1-x3) + Eta_Gauss[j]*(x2-x3) + x3;
                    phiY = Xi_Gauss[j]*(y1-y3) + Eta_Gauss[j]*(y2-y3) + y3;
                    Fe(i)+= weights_Gauss[j]*(std::cos(7*M_PI*phiX)*getBasisFunction(Xi_Gauss[j], Eta_Gauss[j],(i+1))) + weights_Gauss[j]*(std::cos(7*M_PI*phiY)*getBasisFunction(Xi_Gauss[j], Eta_Gauss[j],(i+1)));            
                }   
            } 
                Fe(0) = Fe(0)*Jacobian;
                Fe(1) = Fe(1)*Jacobian;      
                Fe(2) = Fe(2)*Jacobian;              
        }
}

//************************************************//
// Function to compute the global force vector
void ComputeGlobalForceVector(int FunctionIdentifier, Eigen::VectorXd &GlobalForceVec, const Mesh &M)
{
    Eigen::Vector3d Fe;
    int a,b,c;
    for (int i = 1; i<M.nElements+1; i++)
        {
            Fe.setZero();
            a = Local2GlobalMap(M, i,1);
            b = Local2GlobalMap(M, i,2);
            c = Local2GlobalMap(M, i,3);

            ComputeElementForceVector(FunctionIdentifier, Fe, M, i);

            // Assembling the global force vector
            GlobalForceVec.coeffRef(a) += Fe(0);
            GlobalForceVec.coeffRef(b) += Fe(1);
            GlobalForceVec.coeffRef(c) += Fe(2);
        }
    std::cout<<"\nGlobal force vector computed successfully";
}

//************************************************//
// Function to set the Dirichlet boundary conditions
void SetDirichletBcs(Eigen::SparseMatrix<double, Eigen::RowMajor, int> &GlobalStiffnessMatrix, Eigen::VectorXd &GlobalForceVec, const Mesh &M)
{
    int node;
    for (int i =0; i<M.nBdNodes; i++)
        {
            node = M.bdNodes[i];
            GlobalForceVec.coeffRef(node) = 0.0;                                    // Setting the corresponding rows of global force vector to zero
            for(int j=0;j<GlobalStiffnessMatrix.rows(); j++)
                GlobalStiffnessMatrix.coeffRef(node,j) = 0.0;                       // Setting the corresponding rows of global stiffness matrix to zero
            
            GlobalStiffnessMatrix.coeffRef(node,node) = 1.0;                        // Setting the corresponding diagonal entry of stiffness matrix to one
        }
    std::cout<<"\nDirichlet boundary conditions setup successfully";
}

//************************************************//
// Function to print nodal solutions to dat file
void PrintSolutions(std::string solutionfilename, const Mesh &M, const Eigen::VectorXd DOFSolution)
{
    // Write nodal DOF solutions to file
    std::fstream fout;
    fout.open((solutionfilename), std::ios::out);
    for (int i = 0; i < M.nNodes; i++)
        {
            for (int j = 0; j<2; j++) 
                fout << M.Coordinates[i][j]<<" ";

        fout << DOFSolution(i) << "\n";
        }
    fout.close();
}

//************************************************//
// Function to compute the analytical solution at a point
double getAnalyticalSolutionAtPoint(double x, double y, int numterms, int FunctionIdentifier)
{
    double Ua, seriessum;
    if (FunctionIdentifier == 1)
        {
            seriessum = 0.0;
        for (int k=1; k<numterms; k++)
        {
            seriessum += ((std::sinh((2.0*k-1)*M_PI*(1.0-y)) + std::sinh((2.0*k-1)*M_PI*y))/(std::sinh((2.0*k-1.0)*M_PI))) * ((std::sin((2.0*k-1.0)*M_PI*x))/(std::pow((2.0*k - 1.0),3)));
        }
        Ua = (1.0 - x)*x - (8.0/std::pow(M_PI,3)) * seriessum;
        }

    return Ua;
}

//************************************************//
// Function to compute the analytical solution vector
void getAnalyticalSolution(Eigen::VectorXd &AnalyticalSolution, const Mesh &M, int numterms, int FunctionIdentifier)
{  
    double x, y;
    for(int i = 0; i<M.nNodes; i++)
        {
            x = M.Coordinates[i][0];
            y = M.Coordinates[i][1];
            AnalyticalSolution.coeffRef(i) = getAnalyticalSolutionAtPoint(x,y,numterms, FunctionIdentifier);
        }
}

//************************************************//
// Function to compute the error
double getError( const Eigen::VectorXd AnalyticalSolution, Eigen::VectorXd DOFSolution, Eigen::VectorXd &ErrorPointWise)
{
    double NodeWiseErrorCummulative = 0;
    int nNodes = DOFSolution.size();
    double Ua, Uh, E;
    ErrorPointWise = AnalyticalSolution - DOFSolution;

    for (int i=0; i<DOFSolution.size(); i++)
    {
        Ua = AnalyticalSolution.coeff(i);
        Uh = DOFSolution.coeff(i);
        NodeWiseErrorCummulative += std::pow((Ua-Uh),2);
    }
    E = std::sqrt(NodeWiseErrorCummulative/nNodes);                     // Total error

    return E;
}

//************************************************//
// Main function
//************************************************//
int main()
{
    int FunctionIdentifier = 1;                         // Variable to choose the different forcing fucntions
                                                        // FunctionIdentifier = 0 --> f = 0
                                                        // FunctionIdentifier = 1 --> f = 2
                                                        // FunctionIdentifier = 2 --> f = sin(3*pi*x)
                                                        // FunctionIdentifier = 3 --> f = sin(3*pi*x) + sin(3*pi*y)
                                                        // FunctionIdentifier = 4 --> f = cos(7*pi*x) + cos(7*pi*y)

    std::string coordFileName = "coordinates.dat";      // File with details on the nodal coordinates
    std::string connFileName =  "connectivity.dat";     // File with details on the element connectivity
    std::string bdFileName = "boundaryNodes.dat";       // File with details on the boundary nodes
    std::string solutionfilename = "sol.dat";           // File with details of the approximated nodal solution
    Mesh M;         

    ReadMesh(coordFileName,connFileName,bdFileName,M);  // Reading the data files

    Eigen::SparseMatrix<double, Eigen::RowMajor, int> GlobalStiffnessMatrix (M.nNodes, M.nNodes); // Stiffness matrix for the system
    Eigen::VectorXd GlobalForceVec(M.nNodes); 
    std::vector<Eigen::Triplet<double> > Triplets;
    Triplets.reserve(9*M.nElements);                                                               // Reserving memory for the triplets

    ComputeGlobalStiffnessMatrix(Triplets,GlobalStiffnessMatrix,M);         // Function call to compute the global stiffness matrix
    
    GlobalForceVec.setZero();                                               // Initialization of the global force vector
    ComputeGlobalForceVector(FunctionIdentifier, GlobalForceVec, M);        // Function call to compute the global force vector
    
    SetDirichletBcs(GlobalStiffnessMatrix, GlobalForceVec, M);              // Function call to set the Dirichlet BCs

    // Solver creation and DOF vector creation
    Eigen::SparseLU<Eigen::SparseMatrix< double, Eigen::RowMajor> > Solver; // Object for the Sparse solver
    Eigen::VectorXd DOFSolution (M.nNodes);                                 // Vector to store the solution at the nodes

    // Solving the linear system of equations
    Solver.analyzePattern(GlobalStiffnessMatrix);
    Solver.factorize(GlobalStiffnessMatrix);
    DOFSolution = Solver.solve(GlobalForceVec);
    std::cout<<"\nSolution computed successfully";

    // Printing solution to file
    PrintSolutions(solutionfilename, M, DOFSolution);
    std::cout<<"\nSolution file created successfully";

    if (FunctionIdentifier ==1)                                                      //Computing analytical solution for forcing function = 1.0
    {
        Eigen::VectorXd AnalyticalSolution(M.nNodes); 
        Eigen::VectorXd ErrorPointWise(M.nNodes);
        AnalyticalSolution.setZero();
        ErrorPointWise.setZero();

        double Error = 0;
        int numterms = 50;                                                          //Number to terms to approximate the analytical solution.
        getAnalyticalSolution(AnalyticalSolution,M,numterms,FunctionIdentifier);

        Error = getError(AnalyticalSolution,DOFSolution,ErrorPointWise);
        std::cout<<"\nAnalytical solution exist and the error is computed below";
        std::cout<<"\nTotal Error = "<<Error;
    }                                                      
    std::cout <<"\nProgram terminated successfully";

    return 0;
}
//************************************************//
//                  End of Code
//************************************************//
