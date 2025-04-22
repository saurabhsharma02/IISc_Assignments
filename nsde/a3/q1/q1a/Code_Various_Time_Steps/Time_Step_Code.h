#include <iostream>
#include <vector>
#include <cmath>
#include "math.h"
#include <fstream>
#include "string.h"
#include "stdlib.h"
#include <sstream>
#include "assert.h"
#include "iomanip"
#include <algorithm>

// Default namespace
using namespace std;

// Custom datatypes
typedef double real;
typedef vector<real> Vector1D;
typedef vector<vector<real>> Vector2D;
typedef long unsigned int LUint;

// Constants
const real PI = acos(-1);

// Function declarations
void writeVector(Vector1D u, string filename); // write the arrays into a file
void print(Vector2D matrix);				   // print 2D vector to stdout
void print(Vector1D column);				   // print 1D vector to stdout
real rms(Vector1D column);

// ODE functions
class PDE
{
public:
	// Vector declarations
	Vector1D u;
	Vector1D unew;
	vector<LUint> savesteps;

	// Simulation extents
	real Length;
	LUint TNx;
	real Dx;
	real Dt;
	real t_begin;
	real t_end;
	real alpha;
	LUint nSteps;

	// Physical domain
	real PNx1;
	real PNx2;

	PDE() // Constructor; called while object creation
	{
		cout << "\t\t===== START of simulation =====\n";
	}

	~PDE() // Destructor; called when an object goes out of scope
	{
		cout << "\n\t\t===== END of simulation =====\n\n";
	}

	// Set Vector sizes
	void defineProblem(real lx, LUint gridx, real timestep, real t_i, real t_f)
	{
		Length = lx;
		TNx = gridx;
		Dt = timestep;
		t_begin = t_i;
		t_end = t_f;
	}

	void prepareSystem()
	{
		// Physical Space: Start and End
		PNx1 = 0;
		PNx2 = TNx - 1;

		// Stepsizes
		Dx = Length / (real)(TNx - 1); // Delta x
		nSteps = (LUint)((t_end - t_begin) / Dt) + 1;

		cout << "Dx = " << Dx << "; Dt = " << Dt << "; nSteps = " << nSteps << endl;

		// Revise t_end according to Dt
		cout << "Revising the end time according to timestep size...\n";
		t_end = nSteps * Dt;
		cout << "Revised t_end = " << t_end << endl;

		// Resizing vectors
		u.resize(TNx);
		unew.resize(TNx);

		// Setting the given initial condition
		cout << "Setting Initial condition...\t\t\t\t\t";
		real x;
		for (LUint i = PNx1; i <= PNx2; i++)
		{
			x = i * Dx;			
			unew[i] = sin(4 * PI * x) + sin(6 * PI * x) + sin(10 * PI * x);
			if (fabs(unew[i]) < 1e-9)
				unew[i] = 0;
		}
		cout << "Done\n";
	}

	// Solve the inviscid problem
	void computeInviscid()
	{
		prepareSystem();

		// Write the solution data to files
		ostringstream sobj_grid, sobj_time;
		ofstream fout_grid, fout_time;
		string filename_grid, filename_func;
		string gridstring, timestring;
		filename_grid = "../data/";
		sobj_grid << fixed << setprecision(0) << TNx;
		gridstring = sobj_grid.str();

		// Start the computation
		cout << "Computing solution...\t\t\t\t\t\t";
		for (LUint n = 0; n <= nSteps; n++)
		{
			// cout << "Current timestep = " << (n + 1) << "\t";

			// The if-else conditions have been used to keep the finite difference
			// representation always upwind.

			u = unew;
			// Left vertex
			if (u[PNx1] >= 0)
				unew[PNx1] = u[PNx1] * (1 - (Dt / Dx) * (u[PNx1] - u[PNx2 - 1])); // Upwind
			else
				unew[PNx1] = u[PNx1] * (1 - (Dt / Dx) * (u[PNx1 + 1] - u[PNx1])); // Upwind

			// Interior points
			for (LUint i = PNx1 + 1; i <= PNx2 - 1; i++)
			{
				if (u[i] >= 0)
					unew[i] = u[i] * (1 - (Dt / Dx) * (u[i] - u[i - 1])); // Upwind
				else
					unew[i] = u[i] * (1 - (Dt / Dx) * (u[i + 1] - u[i])); // Upwind
			}

			// Right vertex
			if (u[PNx2] >= 0)
				unew[PNx2] = u[PNx2] * (1 - (Dt / Dx) * (u[PNx2] - u[PNx2 - 1])); // Upwind
			else
				unew[PNx2] = u[PNx2] * (1 - (Dt / Dx) * (u[PNx1 + 1] - u[PNx2])); // Upwind

			// Continue the filesaving
			if (binary_search(savesteps.begin(), savesteps.end(), n))
			{
				filename_func = "../data/alpha_0.000_u_N_" + gridstring + "_nStep_";
				sobj_time.str(string()); // clears the string in sobj_time
				sobj_time << fixed << setprecision(0) << n;
				timestring = sobj_time.str();
				filename_func += timestring + ".txt";
				writeVector(unew, filename_func);
			}
		}
		cout << "Done\n";

		// Displaying at what all timesteps the solution has been saved
		cout << "The solution has been saved at timesteps: ";
		for(LUint i = 0; i < savesteps.size(); i++)
		{
			cout << savesteps[i] << " ";
		}
		cout << endl;

		// Write a file containing position vector
		cout << "Writing position vector...\t\t\t\t\t";
		filename_grid = "../data/position_N_" + gridstring + ".txt";
		fout_grid.open(filename_grid);
		for (LUint i = 0; i < TNx; i++)
		{
			fout_grid << i * Dx << " ";
		}
		fout_grid.close();
		cout << "Done\n";
	}

	// Solve the viscid problem
	void computeViscid()
	{
		prepareSystem();

		// Write the solution data to files
		ostringstream sobj_grid, sobj_time;
		ofstream fout_grid, fout_time;
		string filename_grid, filename_func;
		string gridstring, timestring;
		filename_grid = "../data/";
		sobj_grid << fixed << setprecision(0) << TNx;
		gridstring = sobj_grid.str();

		// Start the computation
		cout << "Computing solution...\t\t\t\t\t\t";
		for (LUint n = 0; n <= nSteps; n++)
		{
			// cout << "Current timestep = " << (n + 1) << "\t";

			u = unew;
			// Left vertex			
			unew[PNx1] = u[PNx1] - (u[PNx1] * Dt / (2 * Dx)) * (u[PNx1 + 1] - u[PNx2 - 1])
							+ (alpha * Dt / (Dx * Dx)) * (u[PNx1 + 1] - 2 * u[PNx1] + u[PNx2 - 1]);

			// Interior points
			for (LUint i = PNx1 + 1; i <= PNx2 - 1; i++)
			{
				unew[i] = u[i] - (u[i] * Dt / (2 * Dx)) * (u[i + 1] - u[i - 1])
							+ (alpha * Dt / (Dx * Dx)) * (u[i + 1] - 2 * u[i] + u[i - 1]);
			}

			// Right vertex
			unew[PNx2] = u[PNx2] - (u[PNx2] * Dt / (2 * Dx)) * (u[PNx1 + 1] - u[PNx2 - 1])
							+ (alpha * Dt / (Dx * Dx)) * (u[PNx1 + 1] - 2 * u[PNx2] + u[PNx2 - 1]);

			// Continue the filesaving
			if (binary_search(savesteps.begin(), savesteps.end(), n))
			{
				filename_func = "../data/alpha_0.001_u_N_" + gridstring + "_nStep_";
				sobj_time.str(string()); // clears the string in sobj_time
				sobj_time << fixed << setprecision(0) << n;
				timestring = sobj_time.str();
				filename_func += timestring + ".txt";
				writeVector(unew, filename_func);
			}
		}
		cout << "Done\n";

		// Displaying at what all timesteps the solution has been saved
		cout << "The solution has been saved at timesteps: ";
		for(LUint i = 0; i < savesteps.size(); i++)
		{
			cout << savesteps[i] << " ";
		}
		cout << endl;

		// Write a file containing position vector
		cout << "Writing position vector...\t\t\t\t\t";
		filename_grid = "../data/position_N_" + gridstring + ".txt";
		fout_grid.open(filename_grid);
		for (LUint i = 0; i < TNx; i++)
		{
			fout_grid << i * Dx << " ";
		}
		fout_grid.close();
		cout << "Done\n";
	}
};
