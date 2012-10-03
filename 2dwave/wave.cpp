#include <iostream>
#include <fstream>
#include <cmath>
#include "lib.cpp"
#define PI 3.1415
using namespace std;
double g(double x, double y)
{
	return sin(PI*x)*sin(2*PI*y);
}
void print_matrix(double **v,int n)
{
	ofstream file("test.d");
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++)
			file << v[i][j] << " ";
		file << endl;
	}
	file.close();
		
}
int main()
{
	//parameters -->
	double c = 1;
	double T = 1;
	double Lx = 1;
	double Ly = 1;
	double dt = 0.1;
	double dx = 0.01;
	double dy = 0.01;
	// <--
	int Nt = (int)floor(T/dt);
	int Nx = (int)floor(Lx/dx);
	int Ny = (int)floor(Ly/dy);
	double cc = c*c;
	double dtt = dt*dt;
	double dxx = dx*dx;
	double dyy = dy*dy;
	double **Unext = (double**)matrix(Nx,Ny,sizeof(double));
	double **U = (double**)matrix(Nx,Ny,sizeof(double));
	double **Uold = (double**)matrix(Nx,Ny,sizeof(double));
	// initial cond --> 
	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Ny; j++) {
			Uold[i][j] = g(i*dx,j*dy);
		}
	}
	U[0][0] = 0;
	U[0][Ny-1] = 0;
	U[Nx-1][0] = 0;
	U[Nx-1][Ny-1]=0;
	double U_1 = 0;
	for (int i=1; i<(Nx-1); i++) {
		for (int j=1; j<(Ny-1); j++) {
			U_1 = Uold[i][j] + (dt/(2*dxx))*(-4*Uold[i][j] + Uold[i+1][j]
				 + Uold[i-1][j] + Uold[i][j+1] + Uold[i][j-1]);
			U[i][j] = 2*Uold[i][j] - U_1 + (dtt/dxx)*(-4*Uold[i][j] + 
				Uold[i+1][j] + Uold[i-1][j] + Uold[i][j+1] + Uold[i][j-1]);
		}
		// boundary cond.
		Uold[i][0] = 0;
		Uold[i][Ny-1] = 0;
		Uold[0][i] = 0;
		Uold[Nx-1][i] = 0;
		U[i][0] = 0;
		U[i][Ny-1] = 0;
		U[0][i] = 0;
		U[Nx-1][i] = 0;		
	}
	for (int i=0;i<Nx;i++) {
		for (int j=0;j<Nx;j++) {
			if (U[i][j] > 5)
				cout << i << "," << j << endl;
		}
	}
	print_matrix(U,Nx);
	for (int n=0; n<(Nt-1); n++) {
		for (int i=1; i<(Nx-1); i++) {
			for (int j=1; j<(Ny-1); j++) {
	Unext[i][j] = 2*U[i][j] - Uold[i][j] + (dtt/dxx)*(-4*U[i][j] + U[i+1][j] + U[i-1][j] 
			+ U[i][j+1] + U[i][j-1]);
			}
		}
		for (int i=1; i<(Nx-1); i++) {
			for (int j=1; j<(Ny-1); j++) {
				Uold[i][j] = U[i][j];
				U[i][j] = Unext[i][j];
			}
		}
	}


	return 0;
}
