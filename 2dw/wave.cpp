#include "various.cpp"
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstdlib>
using namespace std;
using namespace boost::numeric::ublas;

double T, Lx, Ly, dt;
//double b;
int Nx,Ny,Nt;
double dx, dy;
double ddx, ddy, ddt;
double(*c)(double x, double y);
double(*source)(double x, double y,double t);
double(*initc)(double x, double y);
double(*V)(double x, double y);
//double deltaDiffB(int i1, int i, int i_1, int j1, int j, int j_1, int n, matrix<double> *uc);
//void diff(matrix<double> *u,matrix<double> *uc,matrix<double> *ul, double a, double b, double c,int n);
double diff_2(int i, int j, int n, matrix<double> *uc, matrix<double> *ul);
double RHS(int i, int j, int n, matrix<double> * uc);
double firstStep(int i, int j, int n, matrix<double> *uc);
int main(int argc, char** argv)
{
    stringstream str;
    int write_delay;

    // Input parameters ----->
    T = 1;
    write_delay = 100;
    Lx = 4;
    Ly = 4;
    b = 0.01;
    dt = 0.001;
    dx = 0.08;
    dy = 0.08;
    // cmd line override
    if (argc==2) T = atof(argv[1]);

    // functions
    initc = &exact_initc;
    c = &exact_coeff;
    source = &exact_source;
    V = &exact_V;
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
    matrix<double> exact(Nx,Ny);
    matrix<double> error(Nx,Ny);

    double E;
    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {
            u(i,j) = 0;
            unext(i,j) = 0;
            ulast(i,j) = initc(dx*i,dy*j);
        }
    }
    // n = 0 step (first time step)
    #pragma omp parallel for
    for (int i = 0; i<Nx ;i++) {
        for (int j=0; j<Ny; j++) {
            u(i,j) = ulast(i,j) + dt*V(i*dx,j*dy)*(1-b*dt/2) + 0.5*RHS(i,j,0,&ulast)*ddt;
            //exact(i,j) = exact_sol(i*dx,j*dy,0);
        }
    }
    // time loop
    for (int n=1; n<Nt; n++) {
    #pragma omp parallel for
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
                unext(i,j) = (1/(1+b*dt/2.))*( 2*u(i,j) - ulast(i,j)*(1-b*dt/2.) + RHS(i,j,n,&u)*ddt);
                //exact(i,j) = exact_sol(i*dx,j*dy,n*dt);
            }
        }
        E = 0;
    #pragma omp parallel for
        for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
                ulast(i,j) = u(i,j);
                u(i,j) = unext(i,j);
                //error(i,j) = u(i,j) - exact(i,j);
                //E += (u(i,j) - exact(i,j))*(u(i,j) - exact(i,j));
            }
        }
        E = dx*dy*sqrt(E);

        if (n % write_delay == 0) {
            str.str(std::string(""));
            str << "test.d" << padnumber(n/write_delay,'0');
            //print_matrix(Nx,Ny,&u,str.str().c_str(),1);
            // print error
            //cout << "E: " << E << endl;
        }
        // cout << "E: " << E << endl;
    }

    return 0;
}
double RHS(int i, int j, int n, matrix<double> * uc)
{
    int i1 = i+1; int i_1 = i-1; int j1 = j+1; int j_1 = j-1;
    if (i_1 < 0) i_1 = i1;
    if (j_1 < 0) j_1 = j1;
    if (i1 > Nx-1) i1 = i_1;
    if (j1 > Ny-1) j1 = j_1;
    return source(i*dx,j*dy,n*dt) + (1/ddx)*(  c((i+0.5)*dx,j*dy)*((*uc)(i1,j) - (*uc)(i,j)) - c((i-0.5)*dx,j*dy)*((*uc)(i,j) - (*uc)(i_1,j)) )
                                  + (1/ddy)*(  c(i*dx,(j+0.5)*dy)*((*uc)(i,j1) - (*uc)(i,j)) - c(i*dx,(j-0.5)*dy)*((*uc)(i,j) - (*uc)(i,j_1)) );
}
