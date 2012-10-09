#include <WaveSolver.h>
#include <CWave.h>
#include <CIniFile.h>

// Mods the index on system size
inline int WaveSolver::idx(int i) {
	return (i+10*Nr)%Nr;
}

// Calculates the value of u_(i+di,j+dj) but takes care of the boundary conditions and world
inline double WaveSolver::u(int i,int j, int di, int dj) {
	if(walls(idx(i+di),idx(j+dj))) {
		return u_(idx(i-di),idx(j-dj));
	} 

	return u_(idx(i+di),idx(j+dj));
}

// Calculates the value of u_prev(i+di,j+dj) but takes care of the boundary conditions and world
inline double WaveSolver::uprev(int i,int j, int di, int dj) {
	if(walls(idx(i+di),idx(j+dj))) {
		return u_prev(idx(i-di),idx(j-dj));
	} 

	return u_prev(idx(i+di),idx(j+dj));
}

// Gives the c-value in the point i,j
inline double WaveSolver::calcC(int i, int j) {
	i=(i+10*Nr)%Nr;
	j=(j+10*Nr)%Nr;
	return max(1-ground(i,j),1.0);
}



void WaveSolver::copyToGrid(AObject& grid) {
	for (int i=0;i<Nr;i++)
		for (int j=0;j<Nr;j++) {
			grid.getGridPos(i,j)->pos.z = u_(i,j);
		}	
}

void WaveSolver::changeGroundZIncrease(double deltaZ) {
	ground += deltaZ;
	for(int i=0;i<Nr;i++) {
		for(int j=0;j<Nr;j++) {
			walls(i,j) = ground(i,j) >= 1.0;
		}
	}
}

WaveSolver::WaveSolver(CIniFile &ini) {
	Nr = ini.getint("grid_size");
	dampingFactor = ini.getdouble("damping");
	walls = zeros<mat>(Nr,Nr);

	// Read in the ground info from bmp file
	// The values can be adjusted and scaled in the wave.ini
	ground = readBMP((char*)ini.getstring("ground_file").c_str());
	ground += ini.getdouble("ground_z_increase");
	ground *= ini.getdouble("ground_z_scale");;
	for(int i=0;i<Nr;i++) {
		for(int j=0;j<Nr;j++) {
			walls(i,j) = ground(i,j) >= 1.0;
		}
	}
	
	max_value = 0;
	r_min = -1.0;
	r_max = 1.0;               			// Square grid
	time = 0;

	dr = (r_max-r_min)/(Nr-1);		 	// Step length in space
	double c_max = 1.0;       			// Used to determine dt and Nt

	double k = 0.5;               		// We require k<=0.5 for stability
	dt = 0.9*dr/sqrt(2*c_max); 			// This guarantees (I guess) stability if c_max is correct
	
	dtdt_drdr = dt*dt/(dr*dr); 			// Constant that is used in the calculation
	
	x  = zeros<vec>(Nr,1); 				// Positions to calculate initial conditions
	y  = zeros<vec>(Nr,1);				

	u_next = zeros<mat>(Nr,Nr);			// U_{n+1}
	u_ = zeros<mat>(Nr,Nr);  			// U_{n}
	u_prev = zeros<mat>(Nr,Nr); 		// U_{n-1}

	double _x,_y; // Temp variables for speed'n read

	// Calculate initial conditions, gaussian at x0,y0
	double x0 = ini.getdouble("x0");
	double y0 = ini.getdouble("y0");
	double stddev = ini.getdouble("stddev");

	for(int i=0;i<Nr;i++) {
		x(i) = r_min+i*dr;
		for(int j=0;j<Nr;j++) {
			y(j) = r_min+j*dr; 

			_x = x(i)-x0; 					// The x- and y-center can have an offset
			_y = y(j)-y0;

			u_prev(i,j) = exp(-(pow(_x,2)+pow(_y,2))/(2*stddev*stddev));
			u_(i,j)     = exp(-(pow(_x,2)+pow(_y,2))/(2*stddev*stddev));
			max_value = max(max_value,abs(u_(i,j)));
		}
	}

	// Normalize these to the amplitude we have chosen in ini
	double amplitude = ini.getdouble("amplitude");
	u_prev *= amplitude/max_value;
	u_ *= amplitude/max_value;
	max_value = amplitude;
}

void WaveSolver::step() {
	time += dt;
	
	double cx_m, cx_p,cy_m, cy_p, c, ddx, ddy, ddt_rest, source; 		// Temp variables for speed'n read
	double factor = 1.0/(1+0.5*dampingFactor*dt);
	int i,j;

	#pragma omp parallel for private(cx_m, cx_p,cy_m, cy_p, c, ddx, ddy, ddt_rest, source,i,j) num_threads(4)
	for(i=0;i<Nr;i++) {
		for(j=0;j<Nr;j++) {
			c = calcC(i,j);

			cx_m = 0.5*(c+calcC(i-1,j)); 	// Calculate the 4 c's we need. We need c_{i \pm 1/2,j} and c_{i,j \pm 1/2}
			cx_p = 0.5*(c+calcC(i+1,j));
			cy_m = 0.5*(c+calcC(i,j-1));
			cy_p = 0.5*(c+calcC(i,j+1));

			ddx = cx_p*( u(i,j,1,0)   - u(i,j) ) - cx_m*( u(i,j) - u(i,j,-1,0) );
			ddy = cy_p*( u(i,j,0,1)   - u(i,j) ) - cy_m*( u(i,j) - u(i,j,0,-1) );
			ddt_rest = -(1-0.5*dampingFactor*dt)*uprev(i,j) + 2*u(i,j);
			source = 0;

			// Set value to zero if we have a wall.
			// u_next(i,j) = walls(i,j) ? 0 : factor*(dtdt_drdr*(ddx + ddy) + ddt_rest + source);
			u_next(i,j) = walls(i,j) ? 0 : factor*(dtdt_drdr*(ddx + ddy) + ddt_rest + source);
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
    
    int width  = *(int*)&info[18];
    int height = *(int*)&info[22];
    int pixels = width*height;

    int bytesLeft = fileSize-54; // 54 bytes for the header, we need the rest

    unsigned char* data = new unsigned char[bytesLeft];
    fread(data, sizeof(unsigned char), bytesLeft, f); // read the rest of the data at once
    fclose(f);

    mat img = zeros<mat>(height,width);

    int pixelCount = 0;
    int pixelIndex = 0;

    int row = 0; // BMP-files start with the lower left pixel, row for row
    int col = 0;

    for(i = 0; i < pixels; i++)
    {
		int r = data[pixelIndex+2]; // RGB values are interchanged in the file format
		int g = data[pixelIndex+1];
		int b = data[pixelIndex];

		double avg = 1.0*(r+g+b)/3.0/255.0; // If we have black/white only, the average is 0 (black) to 255 (white)
		img(col,row) = avg;

		pixelCount++;
		pixelIndex+=3; // Each pixel has 3 bytes, RGB
		
		if(++col == width) {
			// BMP-format is stupid since it wants each row to have 4n bytes, so it adds 
			// the remaining pixels before next row
			int padding = col % 4;

			col = 0;
			pixelIndex+=padding;
			row++;
		}
    }

    return img;
}