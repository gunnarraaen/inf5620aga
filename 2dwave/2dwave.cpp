#include "bitmap/EasyBMP.cpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "tools.cpp"
#include <cstdlib>
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
		unext = (double**)matrix(Nx,Nx,sizeof(double));
		ulast = (double**)matrix(Nx,Nx,sizeof(double));
		for (int i=0; i<Nx; i++) {
			for (int j=0; j<Nx; j++) {
				ulast[i][j] = 0;
				unext[i][j] = 0;
				u[i][j] = (*initc)(h*i,h*j); 
			}
		}
		// artificial n=-1 step
		wave(ulast,u,ulast,0.5,0,0.5);

	}
	~wavesolver() {
		free_matrix((void**)unext);
		free_matrix((void**)ulast);
	}
	int nextstep()
	{
		wave(unext,u,ulast,1,1,1);
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
			//							      i+1,i,i-1		 j+1,j,j-1
			un[0][j] =    2*a*uc[0][j]    - b*ul[0][j]     + c*deltaDiffB(1,0,1		,j+1,j,j-1	,uc); // i = 0
			un[Nx-1][j] = 2*a*uc[Nx-1][j] - b*ul[Nx-1][j]  + c*deltaDiffB(Nx-2,Nx-1,Nx-2	,j+1,j,j-1	,uc); // i = Nx-1
			un[j][0] =    2*a*uc[j][0]    - b*ul[j][0]     + c*deltaDiffB(j+1,j,j-1		,1,0,1		,uc); // j = 0
			un[j][Nx-1] = 2*a*uc[j][Nx-1] - b*ul[j][Nx-1]  + c*deltaDiffB(j+1,j,j-1		,Nx-2,Nx-1,Nx-2	,uc); // j = Nx-1
		}
		// corners
		un[0][0] =       2*a*uc[0][0]       - b*ul[0][0]       + c*deltaDiffB(1,0,1		,1,0,1		,uc); // i = 0, j = 0
		un[0][Nx-1] =    2*a*uc[0][Nx-1]    - b*ul[0][Nx-1]    + c*deltaDiffB(1,0,1		,Nx-2,Nx-1,Nx-2	,uc); // i = 0, j = Nx-1
		un[Nx-1][0] =    2*a*uc[Nx-1][0]    - b*ul[Nx-1][0]    + c*deltaDiffB(Nx-2,Nx-1,Nx-2	,1,0,1		,uc); // i = Nx-1, j = 0
		un[Nx-1][Nx-1] = 2*a*uc[Nx-1][Nx-1] - b*ul[Nx-1][Nx-1] + c*deltaDiffB(Nx-2,Nx-1,Nx-2	,Nx-2,Nx-1,Nx-2	,uc); // i = Nx-1, j = Nx-1
		
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
	double deltaDiffB(int i1, int i, int i_1, int j1, int j, int j_1, double **uc) 
	{
		return 	(dt*dt)/(h*h)*( 
			(*c)((i+0.5)*h,j*h)*(uc[i1][j] - uc[i][j]) - (*c)((i-0.5)*h,j*h)*(uc[i][j] - uc[i_1][j]) + 
			(*c)(i*h,(j+0.5)*h)*(uc[i][j1] - uc[i][j]) - (*c)(i*h,(j-0.5)*h)*(uc[i][j] - uc[i][j_1])  
			);
	}
};
const double pi = 3.1415;
double initcond_exp(double x,double y)
{
//	return 4*exp( -((x-1.2)*(x-1.2) + (y-1.5)*(y-1.5))*90 ) + 4*exp( -((x-1.2)*(x-1.2) + (y-2.5)*(y-2.5))*90 );
	return 10*exp( -((x-2)*(x-2) + (y-2)*(y-2))*40 );
;
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
	if (x>1.6 && x<1.7) 
	{
		if ((y>1.5 && y<1.6) || (y>2.4 && y<2.5))
			return 1;
		else 
			return 0;
	}
	return 1;
}
double coeff2(double x, double y)
{
	return 1 - exp( - ((x-2)*(x-2) + (y-2)*(y-2))*20 ) - exp( -((x-3)*(x-3) + (y-3)*(y-3))*20 );
}
double coeff3(double x, double y)
{
	return 1.0;

	if (x<0 || y<0) return 0;
	int Lx = 4;
	static bool used_once = false;
	static BMP input;
	if (!used_once)
	{
		input.ReadFromFile("test4.bmp");
		used_once = true;
	}
	int relevant_x = x*input.TellHeight()/Lx;
	int relevant_y = y*input.TellWidth()/Lx;
	double val = (int) input(relevant_y,relevant_x)->Red + (int) input(relevant_y,relevant_x)->Green + (int) input(relevant_y,relevant_x)->Blue;
//	cout << "x:" << relevant_x << " y:" << relevant_y << endl;
	val /=(3*256.0);
	
	//cout << relevant_x << ", " << relevant_y << ":" << val << "\n";
	
	return val;
}
int main(int argc, char** argv)
{
	
/*	BMP input;
	input.ReadFromFile("test3.bmp");
	cout << "   (400,100)  (400,200)  (400,300)" << endl;
	cout << "Blue: " << (int) input(400,150)->Blue << " " << (int) input(400,300)->Blue << " " << (int) input(400,420)->Blue << endl;
	cout << "Red: "  << (int) input(400,150)->Red << " " << (int) input(400,250)->Red << " " << (int) input(400,420)->Red << endl;
	cout << "Green: "<< (int) input(400,150)->Green << " " << (int) input(400,250)->Green << " " << (int) input(400,420)->Green << endl;
	return 0;*/

	double T = 1;
	if (argc==2) {
		T = atof(argv[1]);
		cout << endl << argc;
	}
	stringstream str;
	int write_delay = 20;
	double Lx = 4;
	double dt = 0.0005;
	double h = 0.01;
	int Nx = (int) ceil(Lx/h);
	int Nt = (int) ceil(T/dt);
	double **u = (double**)matrix(Nx,Nx,sizeof(double));
	wavesolver wave1(&coeff3,&initcond_exp,u, dt, h, T, Lx);
	for (int n=0; n<Nt; n++) {
		wave1.nextstep();
/*		for (int i=0; i<Nx; i++) {
		for (int j=0; j<Nx; j++) {
			u[i][j] = initcond_sin_exact(i*h,j*h,n*dt,0.75);
		}}*/
		if (n % write_delay == 0) {
			str.str(std::string(""));
            str << "test.d" << padnumber(n/write_delay);
			print_matrix(Nx,u,str.str().c_str(),2);
		}
	}
	free_matrix((void**)u);
	return 0;
}
