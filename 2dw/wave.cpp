#include "various.cpp"
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstdlib>
using namespace std;
using namespace boost::numeric::ublas;

double T, Lx, Ly, dt;
int Nx,Ny,Nt;
double dx, dy;
double ddx, ddy, ddt;
double(*c)(double x, double y);
double(*source)(double x, double y,double t);
double(*initc)(double x, double y);
double deltaDiffB(int i1, int i, int i_1, int j1, int j, int j_1, int n, matrix<double> *uc);
void diff(matrix<double> *u,matrix<double> *uc,matrix<double> *ul, double a, double b, double c,int n);

int main(int argc, char** argv)
{
    stringstream str;
    int current_time_step = 0;
    int write_delay;

    // Input parameters ----->
    T = 1;
    write_delay = 10;
    Lx = 4;
    Ly = 4;
    dt = 0.001;
    dx = 0.04;
    dy = 0.04;
    // cmd line override
    if (argc==2) T = atof(argv[1]);

    // functions
    initc = &initcond_exp;  // init cond.
    c = &coeff2;            // coeff.
    source = &small_source;    // source
    // <----------------------
    ddx = dx*dx;
    ddy = dy*dy;
    ddt = dt*dt;
    Nx = (int) ceil(Lx/dx);
    Ny = (int) ceil(Ly/dy);
    Nt = (int) ceil(T/dt);
    matrix<double> u(Nx,Ny);
    matrix<double> unext(Nx,Ny);
    matrix<double> ulast(Nx,Ny);

    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {
            ulast(i,j) = 0;
            unext(i,j) = 0;
            u(i,j) = (*initc)(dx*i,dy*j);
        }
    }
    // artificial n=-1 step
    diff(&ulast,&u,&ulast,0.5,0,0.5,0);

    // time loop
    #pragma omp parallel
    for (int n=0; n<Nt; n++) {
        diff(&unext,&u,&ulast,1,1,1,dt*n);
        for (int i=1; i<Nx-1; i++) {
            for (int j=1; j<Ny-1; j++) {
                ulast(i,j) = u(i,j);
                u(i,j) = unext(i,j);
            }
        }
        current_time_step++;
        if (n % write_delay == 0) {
            str.str(std::string(""));
            str << "test.d" << padnumber(n/write_delay,'0');
            print_matrix(Nx,Nx,&u,str.str().c_str(),1);
        }
    }

    return 0;
}
void diff(matrix<double> *un,matrix<double> *uc,matrix<double> *ul, double a, double b, double c,int n)
{
    for (int i = 0; i < Nx; i++)
        for (int j=0; j < Ny; j++)
            (*un)(i,j) = 2*a*(*uc)(i,j) - b*(*ul)(i,j) + c*deltaDiffB(i+1,i,i-1,j+1,j,j-1,n,uc);
}
double deltaDiffB(int i1, int i, int i_1, int j1, int j, int j_1,int n, matrix<double> *uc)
{
    if (i_1 <0) i_1 = i1;
    if (j_1 <0) j_1 = j1;
    if (i1 > Nx-1) i1=i_1;
    if (j1 > Ny-1) j1=j_1;
    return ddt/ddx*((*c)((i+0.5)*dx,j*dy)*((*uc)(i1,j) - (*uc)(i,j)) - (*c)((i-0.5)*dx,j*dy)*((*uc)(i,j) - (*uc)(i_1,j))) +
           ddt/ddy*((*c)(i*dx,(j+0.5)*dy)*((*uc)(i,j1) - (*uc)(i,j)) - (*c)(i*dx,(j-0.5)*dy)*((*uc)(i,j) - (*uc)(i,j_1))) + dt*dt*source(i*dx,j*dy,n*dt);
}
