#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

int Nr;
vec x;
vec y;
mat H;
mat world;
mat u_p;
mat u_;
mat u_m;
int k; // Used for fast calculation of bounds

bool height_from_file;
bool world_from_file;

mat readBMP(char* filename)
{
    int i;
    FILE* f = fopen(filename, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int fileSize = *(int*)&info[2];
    
    int width = *(int*)&info[18];
    int height = *(int*)&info[22];
    int pixels = width*height;
    int bytesLeft = fileSize-54; // 54 bytes for the header, we need the rest

    unsigned char* data = new unsigned char[bytesLeft];
    fread(data, sizeof(unsigned char), bytesLeft, f); // read the rest of the data at once
    fclose(f);

    mat img = zeros<mat>(height,width);

    int pixelCount = 0;
    int pixelIndex = 0;

    int row = height-1; // BMP-files start with the lower left pixel, row for row
    int col = 0;

    for(i = 0; i < pixels; i++)
    {
    		int r = data[pixelIndex+2]; // RGB values are interchanged in the file format
    		int g = data[pixelIndex+1];
    		int b = data[pixelIndex];

    		double avg = (r+g+b)/3.0/255.0; // If we have black/white only, the average is 0 (black) to 255 (white)
    		img(row,col) = avg;

    		pixelCount++;
    		pixelIndex+=3; // Each pixel has 3 bytes, RGB
    		
    		if(++col == width) {
    			// BMP-format is stupid since it wants each row to have 4n bytes, so it adds 
    			// the remaining pixels before next row
    			int padding = col % 4;

    			col = 0;
    			pixelIndex+=padding;
    			row--;
    		}
    }

    return img;
}

inline int idx(int i) {
	return (i+10*Nr)%Nr;
}

inline double u(int i,int j, int di=0, int dj=0) {
	if(world(i,j)) {
		return u_(idx(i-di),idx(j-dj));
	} 

	return u_(idx(i+di),idx(j+dj));
}

inline double um(int i,int j, int di=0, int dj=0) {
	if(world(i,j)) {
		return u_m(idx(i-di),idx(j-dj));
	} 

	return u_m(idx(i+di),idx(j+dj));
}

inline double calcC(int i, int j) {
	return 1.0;

	i=(i+10*Nr)%Nr;
	j=(j+10*Nr)%Nr;
	return H(i,j);
}

int main(int args, char* argv[]) {
	Nr = args > 1 ? atoi(argv[1]) : 100;
	double t_max = args > 2 ? atof(argv[2]) : 1.0;
	
	height_from_file = args > 3 ? atoi(argv[3]) : false;
	world_from_file = args > 4;
	
	if(height_from_file) {
		char *heightFilename = argv[3];
		H = readBMP(heightFilename);
	}

	if(world_from_file) {
		char *worldFilename = argv[4];
		world = readBMP(worldFilename);
	}
	else world = zeros<mat>(Nr,Nr);

	double r_min = -1.0;
	double r_max = 1.0;

	double dr = (r_max-r_min)/(Nr-1); // Step length in space
	double c_max = 1.0; // Used to determine dt and Nt

	double k = 0.5; // We require k<=0.5 for stability
	double dt = dr/c_max*sqrt(k);
	int Nt = t_max/dt+1;
	double d2 = dt*dt/(dr*dr);
	
	x  = zeros<vec>(Nr,1); // Positions to calculate initial conditions
	y  = zeros<vec>(Nr,1);

	u_m = zeros<mat>(Nr,Nr);
	u_ = zeros<mat>(Nr,Nr);
	u_p = zeros<mat>(Nr,Nr);

	FILE *file = 0;
	char *filename = new char[50];

	double _x,_y,cx_m, cx_p,cy_m, cy_p, c;
	
	// Calculate initial conditions
	double x0 = -0.5;
	double y0 = 0.0;
	double stddev = 0.001;
	for(int i=0;i<Nr;i++) {
		x(i) = r_min+i*dr;
		for(int j=0;j<Nr;j++) {
			y(j) = r_min+j*dr;

			_x = x(i)-x0;
			_y = y(j)-y0;

			u_m(i,j) = 0.05*exp(-(pow(_x,2)+pow(_y,2))/(2*stddev));
			u_(i,j)  = 0.05*exp(-(pow(_x,2)+pow(_y,2))/(2*stddev));
		}
	}

	// Timesteps yo, I actually don't print out the first two timesteps because it's ugly code
	for(int n=1;n<Nt-1;n++) {
		sprintf(filename,"data/test.d%.07d",n-1);
		file = fopen(filename,"w");

		for(int i=0;i<Nr;i++) {
			for(int j=0;j<Nr;j++) {
				c = calcC(i,j);

				cx_m = 0.5*(c+calcC(i-1,j));
				cx_p = 0.5*(c+calcC(i+1,j));
				cy_m = 0.5*(c+calcC(i,j-1));
				cy_p = 0.5*(c+calcC(i,j+1));

				double ddx = cx_p*( u(i,j,1) - u(i,j) ) - cx_m*( u(i,j) - u(i,j,-1) );
				double ddy = cy_p*( u(i,j,0,1) - u(i,j) ) - cy_m*( u(i,j) - u(i,j,0,-1) );

				u_p(i,j) = world(i,j) ? 0 : d2*(ddx + ddy) - um(i,j) + 2*u(i,j);

				/*
				u_p(i,j) = c*d2*( u(i+1,j) + u(i-1,j) - 2*u(i,j) ) 
				 	+ c*d2*( u(i,j+1) + u(i,j-1) -2*u(i,j) )
				 	- um(i,j) + 2*u(i,j);
				*/
				
				fprintf(file,"%f ",u_p(i,j));
			}
			fprintf(file,"\n");
		}
		fclose(file);

		u_m = u_;
		u_ = u_p;
	}

	printf("Simulation finished after %d timesteps\n",Nt);

	return 0;
}