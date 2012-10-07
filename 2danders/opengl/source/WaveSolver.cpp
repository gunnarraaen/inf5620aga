#include <WaveSolver.h>
#include <Armadillo.h>
using namespace arma;

// Mods the index on system size
inline int idx(int i) {
	return (i+10*Nr)%Nr;
}

// Calculates the value of u_(i+di,j+dj) but takes care of the boundary conditions and world
inline double u(int i,int j, int di=0, int dj=0) {
	if(world(idx(i+di),idx(j+dj))) {
		return u_(idx(i-di),idx(j-dj));
	} 

	return u_(idx(i+di),idx(j+dj));
}

// Calculates the value of u_prev(i+di,j+dj) but takes care of the boundary conditions and world
inline double uprev(int i,int j, int di=0, int dj=0) {
	if(world(idx(i+di),idx(j+dj))) {
		return u_prev(idx(i-di),idx(j-dj));
	} 

	return u_prev(idx(i+di),idx(j+dj));
}

// Gives the c-value in the point i,j
inline double calcC(int i, int j) {
	i=(i+10*Nr)%Nr;
	j=(j+10*Nr)%Nr;
	return H(i,j);
}

WaveSolver::WaveSolver() {
	Nr = 200;
	double t_max = 50.0;
	
	height_from_file = false;
	world_from_file = false;
	
	H = ones<mat>(Nr,Nr); 			// H contains values between 0 and 1 which are the velocities in each point
	world = zeros<mat>(Nr,Nr); 		// World contains boolean values 0 or 1 where 1 is wall

	r_min = -1.0;
	r_max = 1.0;               	// Square grid

	dr = (r_max-r_min)/(Nr-1); 	// Step length in space
	double c_max = 1.0;           		// Used to determine dt and Nt

	double k = 0.5;               		// We require k<=0.5 for stability
	dt = dr/c_max*sqrt(k); 		// This guarantees (I guess) stability if c_max is correct
	Nt = t_max/dt; 
	dtdt_drdr = dt*dt/(dr*dr); 	// Constant that is used in the calculation
	
	x  = zeros<vec>(Nr,1); 				// Positions to calculate initial conditions
	y  = zeros<vec>(Nr,1);				

	u_next = zeros<mat>(Nr,Nr);			// U_{n-1}
	u_ = zeros<mat>(Nr,Nr);  			// U_{n}
	u_prev = zeros<mat>(Nr,Nr); 		// U_{n+1}

	double _x,_y; // Temp variables for speed'n read
	
	// Calculate initial conditions
	double x0 = 0.5;
	double y0 = 0;
	double stddev = 0.001;

	for(int i=0;i<Nr;i++) {
		x(i) = r_min+i*dr;
		for(int j=0;j<Nr;j++) {
			y(j) = r_min+j*dr; 

			_x = x(i)-x0; // The x- and y-center can have an offset
			_y = y(j)-y0;

			u_prev(i,j) = 0.05*exp(-(pow(_x,2)+pow(_y,2))/(2*stddev));
			u_(i,j)  = 0.05*exp(-(pow(_x,2)+pow(_y,2))/(2*stddev));
		}
	}

	return 0;
}

void WaveSolver::step() {
	double cx_m, cx_p,cy_m, cy_p, c; // Temp variables for speed'n read

	for(int i=0;i<Nr;i++) {
		for(int j=0;j<Nr;j++) {
			c = calcC(i,j);

			cx_m = 0.5*(c+calcC(i-1,j)); // Calculate the 4 c's we need. We need c_{i \pm 1/2,j} and c_{i,j \pm 1/2}
			cx_p = 0.5*(c+calcC(i+1,j));
			cy_m = 0.5*(c+calcC(i,j-1));
			cy_p = 0.5*(c+calcC(i,j+1));

			double ddx = cx_p*( u(i,j,1) - u(i,j) ) - cx_m*( u(i,j) - u(i,j,-1) );
			double ddy = cy_p*( u(i,j,0,1) - u(i,j) ) - cy_m*( u(i,j) - u(i,j,0,-1) );

			u_next(i,j) = dtdt_drdr*(ddx + ddy) - uprev(i,j) + 2*u(i,j);
		}
	}

	u_prev = u_;
	u_ = u_next;
}

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