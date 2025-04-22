// MODIFIED CODE BY SAURABH SHARMA TO FIND L2 ERROR AND MAX ERROR IN .DAT FILE

#include<iostream>
#include<cstdlib>
#include<fstream>
#include<cstring>
#include<math.h>

using namespace std ;

double const PI = 4.0*atan(1.0); // Value of PI

void Allocate_2D_R(double**& m, int d1, int d2) {
        m=new double* [d1];
        for (int i=0; i<d1; ++i) {
                m[i]=new double [d2];
                for (int j=0; j<d2; ++j)
                        m[i][j]=0.0;
        }
}

int main() {

	cout.flags( ios::dec | ios::scientific );
	cout.precision(5);

	double *x, *y, **uExact, **P, **PINV, *PEig, **u, **RHS ;	 
	double **RHS_Tilde, **u_Tilde, **Temp ;
	double L2_Error = 0.0, Max_Error = 0.0, hx, hy  ;

	int Nx = 200, Ny = 200, i, j, k ;

	ofstream File1("Error.dat",ios::out);
	File1.flags( ios::dec | ios::scientific );
	File1.precision(16) ;
	if(!File1) {cerr<< "Error: Output file couldnot be opened.\n";}


    for(int p=0; p<10; ++p){

	x = new double[Nx+1] ;  y = new double[Ny+1] ;
	Allocate_2D_R(uExact, Nx+1, Ny+1) ; Allocate_2D_R(u, Nx+1, Ny+1) ; 
	Allocate_2D_R(P, Nx-1, Nx-1) ; Allocate_2D_R(PINV, Nx-1, Nx-1) ; PEig = new double[Nx-1] ;
	Allocate_2D_R(RHS, Nx-1, Ny-1) ; Allocate_2D_R(Temp, Nx-1, Ny-1) ; Allocate_2D_R(RHS_Tilde, Nx-1, Ny-1) ; Allocate_2D_R(u_Tilde, Nx-1, Ny-1) ;

	// set a uniform grid.
	for(i = 0 ; i <= Nx ; i++) x[i] = i/double(Nx) ;	
	for(i = 0 ; i <= Ny ; i++) y[i] = i/double(Ny) ;	
	hx = 1.0/double(Nx) ; 
	hy = 1.0/double(Ny) ;

	if( Nx != Ny) { cerr << "Code assumes Nx = Ny" << endl ; exit(1); } 

	// set initial solution and Exact solution for computing error.
	for(i = 0 ; i <= Nx ; i++) {
		for(j = 0 ; j <= Ny ; j++) {
			u[i][j] = 0.0 ; 
			uExact[i][j] = sin(10.0*PI*x[i])*sin(10.0*PI*y[j]) ;
		}
	}

	// Set up eigenvalues and matrices
	for(i = 1 ; i < Nx ; i++) {
		PEig[i-1] = -4.0*sin(PI*0.5*x[i])*sin(PI*0.5*x[i]) ;
		for(j = 1 ; j < Nx ; j++) {
			P[i-1][j-1] = sin(i*PI*x[j]) ; 
			PINV[i-1][j-1] = 2.0*hx*sin(j*PI*x[i]) ;
		}
	}

	// Set the right hand side
	for(i = 1 ; i < Nx ; i++) {
		for(j = 1 ; j < Ny ; j++) {
			RHS[i-1][j-1] = -200.0*PI*PI*sin(10.0*PI*x[i])*sin(10.0*PI*y[j]) ;
			RHS[i-1][j-1] += 1000.0*( sin(10.0*PI*x[i])*sin(10.0*PI*y[j]) ) ;
			RHS[i-1][j-1] *= hx*hx ;
		}
	}

	// compute P^{-1} F Q^{-T}
	for(i = 1 ; i < Nx ; i++) {
		for(j = 1 ; j < Nx ; j++) {
			Temp[i-1][j-1] = 0.0 ;
			for(k = 1 ; k < Nx ; k++) Temp[i-1][j-1] += PINV[i-1][k-1]*RHS[k-1][j-1] ;
		}
	}
	for(i = 1 ; i < Nx ; i++) {
		for(j = 1 ; j < Nx ; j++) {
			RHS_Tilde[i-1][j-1] = 0.0 ;
			for(k = 1 ; k < Nx ; k++) RHS_Tilde[i-1][j-1] += Temp[i-1][k-1]*P[k-1][j-1] ;
		}
	}
	
	// compute U_Tilde
	for(i = 1 ; i < Nx ; i++) {
		for(j = 1 ; j < Nx ; j++) {
			RHS_Tilde[i-1][j-1] /= (PEig[i-1] + PEig[j-1] + 1000.0*hx*hx) ;
		}
	}
	for(i = 1 ; i < Nx ; i++) {
		for(j = 1 ; j < Nx ; j++) {
			Temp[i-1][j-1] = 0.0 ;
			for(k = 1 ; k < Nx ; k++) Temp[i-1][j-1] += P[i-1][k-1]*RHS_Tilde[k-1][j-1] ;
		}
	}

	// compute U = P \bar{U} Q^T
	for(i = 1 ; i < Nx ; i++) {
		for(j = 1 ; j < Nx ; j++) {
			u[i][j] = 0.0 ;
			for(k = 1 ; k < Nx ; k++) u[i][j] += Temp[i-1][k-1]*PINV[k-1][j-1] ;
		}
	}

	


	// ofstream File("Output.dat", ios::out) ;
	// File.flags( ios::dec | ios::scientific );
	// File.precision(16) ;
	// if(!File) {cerr<< "Error: Output file couldnot be opened.\n";}
	
	Max_Error = L2_Error = 0.0 ;

	// File << "TITLE = Flow" << endl << "VARIABLES = X, Y, u, Exact " << endl;
	// File << "Zone T = psi I = " << Ny+1 << " J = " << Nx+1 << endl ;

	for(i = 0 ; i <= Nx ; i++) {
		for(j = 0 ; j <= Ny ; j++) {
			if( fabs(u[i][j] - uExact[i][j] ) > Max_Error) Max_Error = fabs( u[i][j] - uExact[i][j] );
			L2_Error += ( u[i][j] - uExact[i][j] )*( u[i][j] - uExact[i][j] )/ ( ( Nx+1.0 )*(Ny+1.0) ) ;

			// File << x[i] << "\t" << y[j] << "\t" << u[i][j] << "\t" << uExact[i][j] << endl ;
		}

	}

	L2_Error = sqrt(L2_Error) ;
	File1<< Nx<<"\t"<<Ny<<"\t"<<L2_Error<<"\t"<<Max_Error<<"\n";


    Nx=Nx+50;
	Ny=Ny+50;

	}


	// L2_Error = sqrt(L2_Error) ;
	// File.close() ;
	cout << "\n L2 : " << L2_Error << "\t Max : " << Max_Error <<  endl ;

	// Free all the memory and exit
	delete[] x ;  delete [] y ; 

	return 0;
}
