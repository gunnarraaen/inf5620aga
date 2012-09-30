#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "tools.cpp"
using namespace std;
class wavesolver {
public:
	wavesolver(double(*_c)(double x, double y), double (*_initc)(double x,double y),double** _u, double _dt = 0.0001, double _h = 0.001, double _T=1, double _Lx = 1, double _Ly = 1)
	{
		initc = _initc;
		c = _c; 
		Lx = _Lx; Ly = _Ly; dt = _dt; h = _h; T = _T;
		Nx = (int) ceil(Lx/h);
		Nt = (int) ceil(T/dt);
		u = _u;
		current_time_step = 0;
		init();
	}
	int nextstep()
	{
		wave(unext,u,ulast,1,1,1);
	// update
		for (int i=1; i<Nx-1; i++) {
			for (int j=1; j<Nx-1; j++) { 
				ulast[i][j] = u[i][j];
				u[i][j] = unext[i][j];
			}
		}
		return ++current_time_step;
	}
	void wave(double **un, double **uc, double **ul, double a, double b, double c) 
	{
		for (int i = 1; i < (Nx-1); i++) {
			for (int j=1; j<Nx-1; j++) {
				un[i][j] = 2*a*uc[i][j] - b*ul[i][j] + c*deltaDiffB(i+1,i,i-1,j+1,j,j-1,uc);
			}
		}
		// boundary
		for (int j = 2; j < (Nx-1); j++) {
			un[0][j] =    2*a*uc[0][j]    - b*ul[0][j]     + c*deltaDiffB(1,0,1		,j+1,j,j-1	,uc);
			un[Nx-1][j] = 2*a*uc[Nx-1][j] - b*ul[Nx-1][j]  + c*deltaDiffB(Nx-2,Nx-1,Nx-2	,j+1,j,j-1	,uc);
			un[j][0] =    2*a*uc[j][0]    - b*ul[j][0]     + c*deltaDiffB(j+1,j,j-1		,1,0,1		,uc);
			un[j][Nx-1] = 2*a*uc[j][Nx-1] - b*ul[j][Nx-1]  + c*deltaDiffB(j+1,j,j-1		,Nx-2,Nx-1,Nx-2	,uc);
		}
		// corner
		un[0][0] =       2*a*uc[0][0]       - b*ul[0][0]       + c*deltaDiffB(1,0,1		,1,0,1		,uc);
		un[0][Nx-1] =    2*a*uc[0][Nx-1]    - b*ul[0][Nx-1]    + c*deltaDiffB(1,0,1		,Nx-2,Nx-1,Nx-2	,uc);
		un[Nx-1][0] =    2*a*uc[Nx-1][0]    - b*ul[Nx-1][0]    + c*deltaDiffB(Nx-2,Nx-1,Nx-2	,1,0,1		,uc);
		un[Nx-1][Nx-1] = 2*a*uc[Nx-1][Nx-1] - b*ul[Nx-1][Nx-1] + c*deltaDiffB(Nx-2,Nx-1,Nx-2	,Nx-2,Nx-1,Nx-2	,uc);
		
	}
	double getTsteps() {return Nt;}
	double getXsteps() {return Nx;}
private:
	double T, Lx, Ly, dt,h;
	int Nx,Ny,Nt;
	double **u;
	double **unext;
	double **ulast;
	int current_time_step;
	double(*c)(double x, double y);
	double(*initc)(double x, double y);
/*	double deltaDiff(int i, int j, double ** uc) 
	{
		return (dt*dt)/(h*h)*(uc[i+1][j] - 4*uc[i][j] + uc[i-1][j] + uc[i][j+1] + uc[i][j-1]);
	}*/
	double deltaDiffB(int i1, int i, int i_1, int j1, int j, int j_1, double **uc) 
	{
	/*	double ui1,ui_1,uj1,uj_1;
		if (i==0) {
			ui1=uc[1][j];    ui_1 = uc[1][j];
			uj1=uc[i][j+1];  uj_1=uc[i][j-1];
		} else if(i==Nx-1)
		{
			ui1=uc[Nx-2][j]; ui_1 = uc[Nx-2][j];
			uj1=uc[i][j+1];  uj_1=uc[i][j-1];
		} else if(j==0)
		{
			ui1=uc[i+1][j];  ui_1 = uc[i-1][j];
			uj1=uc[i][1];    uj_1=uc[i][1];
		} else if(j==Nx-1)
		{
			ui1=uc[i+1][j];  ui_1 = uc[i-1][j];
			uj1=uc[i][Nx-2]; uj_1=uc[i][Nx-2];
		} else {cout << "error!";}
	*/
		return (dt*dt)/(h*h)*(uc[i1][j] - 4*uc[i][j] + uc[i_1][j] + uc[i][j1] + uc[i][j_1]);
	}
	void init() 
	{
		unext = (double**)matrix(Nx,Nx,sizeof(double));
		ulast = (double**)matrix(Nx,Nx,sizeof(double));
		for (int i=0; i<Nx; i++) {
			for (int j=0; j<Nx; j++) {
				ulast[i][j] = 0;
				unext[i][j] = 0;
				//u[i][j] = 0; 
				u[i][j] = (*initc)(h*i,h*j); 
			}
		}
/*		for (int i=1; i<Nx-1; i++) {
			for (int j=1; j<Nx-1; j++) {
				u[i][j] = (*initc)(h*i,h*j); 
			}
		}*/
		wave(ulast,u,ulast,0.5,0,0.5);
		

	}
};
const double pi = 3.1415;
double initcond_exp(double x,double y)
{
	//return sin(pi*x)*sin(2*pi*y);
	return 3*exp( -((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))*100 );
}
double initcond_sin(double x, double y)
{
	return 2*sin(pi/4*x)*sin(pi/4*y);
}
double initcond_sin_exact(double x, double y, double t, double c)
{
	return 2*sin(pi/4.*x)*sin(pi/4.*y)*cos(pi/2.*c*c*t);
}
double coeff(double x, double y)
{
	return 1;
}
int main()
{
	stringstream str;
	int write_delay = 10;
	double T = 3;
	double Lx = 1;
	double dt = 0.001;
	double h = 0.01;
	int Nx = (int) ceil(Lx/h);
	int Nt = (int) ceil(T/dt);
	double **u = (double**)matrix(Nx,Nx,sizeof(double));
	wavesolver wave1(&coeff,&initcond_exp,u, dt, h, T, Lx);
	for (int n=0; n<Nt; n++) {
		wave1.nextstep();
/*		for (int i=0; i<Nx; i++) {
		for (int j=0; j<Nx; j++) {
			u[i][j] = initcond_sin_exact(i*h,j*h,n*dt,1);
		}}*/
		if (n % write_delay == 0) {
			str.str(std::string(""));
			str << "test.d" << ZeroPadNumber(n/write_delay);
			print_matrix(Nx,u,str.str().c_str());
		}
	}	
	return 0;
}
/*	

	double **unext = (double**)matrix(Nx,Nx,sizeof(double));

	double **ulast = (double**)matrix(Nx,Nx,sizeof(double));

	for (int i=0; i<Nx; i++) {
		for (int j=0; j<Nx; j++) {
			u[i][j] = 0;
			unext[i][j] = 0;
			ulast[i][j] = initcond(h*i,h*j);
		}
	}
	// first timestep n = 0
	double temp;
	for (int i=1; i<Nx-1; i++) {
		for (int j=1; j<Nx-1; j++) {
			temp = ulast[i][j] + (dt)/(2*h*h)*(ulast[i+1][j] - 4*ulast[i][j] + ulast[i-1][j] + ulast[i][j+1] + ulast[i][j-1]);
			u[i][j] = 2*ulast[i][j] - temp + (dt*dt)/(h*h)*(ulast[i+1][j] - 4*ulast[i][j] + ulast[i-1][j] + ulast[i][j+1] + ulast[i][j-1]);
		}
	}
	for (int n=1; n<Nt; n++) {
		for (int i=1; i<Nx-1; i++) {
			for (int j=1; j<Nx-1; j++) {
				unext[i][j] = 2*u[i][j] - ulast[i][j] + (dt*dt)/(h*h)*(u[i+1][j] - 4*u[i][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]);
				}
		}
		for (int i=1; i<Nx-1; i++) {
			for (int j=1; j<Nx-1; j++) { 
				ulast[i][j] = u[i][j];
				u[i][j] = unext[i][j];
			}
		}
		if (n % 100 == 0) {
		str.str(std::string(""));
		str << "test.d" << ZeroPadNumber(n/100);
		print_matrix(Nx,u,str.str().c_str());
		}
	} // time
	*/

