#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

int Nr;

inline int idx(int i) {
	return (i+10*Nr)%Nr;
}

inline double calcC(double x, double y) {
	return 1.0;
}

int main(int args, char* argv[]) {
	Nr = args > 1 ? atoi(argv[1]) : 100;
	double t_max = args > 2 ? atof(argv[2]) : 1.0;

	double r_min = -1.0;
	double r_max = 1.0;

	double dr = (r_max-r_min)/(Nr-1); // Step length in space
	double c_max = 1.1; // Used to determine dt and Nt

	double k = 0.5; // We require k<=0.5 for stability
	double dt = dr/c_max*sqrt(k);
	int Nt = t_max/dt+1;
	double d2 = dt*dt/(dr*dr);
	
	vec x  = zeros<vec>(Nr,1); // Positions to calculate initial conditions
	vec y  = zeros<vec>(Nr,1);

	// um=u^{n-1}, up=u^{n+1}, u=u^{n}
	mat um = zeros<mat>(Nr,Nr);
	mat u = zeros<mat>(Nr,Nr);
	mat up = zeros<mat>(Nr,Nr);

	FILE *file = 0;
	char *filename = new char[50];

	double _x,_y,c;
	
	// Calculate initial conditions
	for(int i=0;i<Nr;i++) {
		x(i) = r_min+i*dr;
		for(int j=0;j<Nr;j++) {
			y(j) = r_min+j*dr;

			_x = x(i);
			_y = y(j);

			um(i,j) = 0.5*exp(-(pow(_x,2)+pow(_y,2))*30);
			u(i,j)  = 0.5*exp(-(pow(_x,2)+pow(_y,2))*30);
		}
	}

	// Timesteps yo, I actually don't print out the first two timesteps because it's ugly code
	for(int n=1;n<Nt-1;n++) {
		sprintf(filename,"data/test.d%.07d",n-1);
		file = fopen(filename,"w");

		for(int i=0;i<Nr;i++) {
			for(int j=0;j<Nr;j++) {
				c = calcC(x(i),y(j));
				up(i,j) = c*d2*( u(idx(i+1),idx(j)) + u(idx(i-1),idx(j)) - 2*u(idx(i),idx(j)) ) 
					+ c*d2*( u(idx(i),idx(j+1)) + u(idx(i),idx(j-1)) -2*u(idx(i),idx(j)) )
					- um(idx(i),idx(j)) + 2*u(idx(i),idx(j));
				fprintf(file,"%f ",up(i,j));
			}
			fprintf(file,"\n");
		}
		fclose(file);

		um = u;
		u = up;
	}

	printf("Simulation finished after %d timesteps\n",Nt);

	return 0;
}