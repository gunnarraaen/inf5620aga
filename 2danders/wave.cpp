#include <iostream>
#include "math.h"
#include "time.h"
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

int Nr;

imat readBMP(char* filename)
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
    printf("Image is %d bytes\n",fileSize);
    printf("Image contains %d pixels (%d,%d)\n",pixels,height,width);

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

			um(i,j) = sin(5*M_PI*exp(-(pow(_x,2)+pow(_y,2))*0.1))*exp(-(pow(_x,2)+pow(_y,2))*7);
			u(i,j)  = sin(5*M_PI*exp(-(pow(_x,2)+pow(_y,2))*0.1))*exp(-(pow(_x,2)+pow(_y,2))*7);
			
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