#include <WaveSolver.h>
#include <CWave.h>
#include <CIniFile.h>

// Mods the index on system size
inline int WaveSolver::idx(int i) {
	return (i+10*Nr)%Nr;
}

// Calculates the value of u_(i+di,j+dj) but takes care of the boundary conditions and world
inline double WaveSolver::u(int i,int j, int di, int dj) {
	if(world(idx(i+di),idx(j+dj))) {
		return u_(idx(i-di),idx(j-dj));
	} 

	return u_(idx(i+di),idx(j+dj));
}

// Calculates the value of u_prev(i+di,j+dj) but takes care of the boundary conditions and world
inline double WaveSolver::uprev(int i,int j, int di, int dj) {
	if(world(idx(i+di),idx(j+dj))) {
		return u_prev(idx(i-di),idx(j-dj));
	} 

	return u_prev(idx(i+di),idx(j+dj));
}

// Gives the c-value in the point i,j
inline double WaveSolver::calcC(int i, int j) {
	i=(i+10*Nr)%Nr;
	j=(j+10*Nr)%Nr;
	return H(i,j);
}

WaveSolver::WaveSolver(CIniFile &ini) {
	render_wall = ini.getbool("render_wall");
	Nr = ini.getint("grid_size");
	dampingFactor = ini.getdouble("damping");
	
	H = readBMP((char*)ini.getstring("velocity_file").c_str());	
	world = readBMP((char*)ini.getstring("wall_file").c_str());	

	max_value = 0;
	r_min = -1.0;
	r_max = 1.0;               			// Square grid
	time = 0;

	dr = (r_max-r_min)/(Nr-1);		 	// Step length in space
	double c_max = H.max();           		// Used to determine dt and Nt

	double k = 0.5;               		// We require k<=0.5 for stability
	dt = dr/sqrt(2*c_max); 				// This guarantees (I guess) stability if c_max is correct
	
	dtdt_drdr = dt*dt/(dr*dr); 			// Constant that is used in the calculation
	
	x  = zeros<vec>(Nr,1); 				// Positions to calculate initial conditions
	y  = zeros<vec>(Nr,1);				

	u_next = zeros<mat>(Nr,Nr);			// U_{n-1}
	u_ = zeros<mat>(Nr,Nr);  			// U_{n}
	u_prev = zeros<mat>(Nr,Nr); 		// U_{n+1}

	double _x,_y; // Temp variables for speed'n read

	// Calculate initial conditions
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
	
	max_value = 0;

	double cx_m, cx_p,cy_m, cy_p, c; 		// Temp variables for speed'n read
	double factor = 1.0/(1+0.5*dampingFactor*dt);
	for(int i=0;i<Nr;i++) {
		for(int j=0;j<Nr;j++) {
			c = calcC(i,j);

			cx_m = 0.5*(c+calcC(i-1,j)); 	// Calculate the 4 c's we need. We need c_{i \pm 1/2,j} and c_{i,j \pm 1/2}
			cx_p = 0.5*(c+calcC(i+1,j));
			cy_m = 0.5*(c+calcC(i,j-1));
			cy_p = 0.5*(c+calcC(i,j+1));

			double ddx = cx_p*( u(i,j,1) - u(i,j) ) - cx_m*( u(i,j) - u(i,j,-1) );
			double ddy = cy_p*( u(i,j,0,1) - u(i,j) ) - cy_m*( u(i,j) - u(i,j,0,-1) );
			double ddt_rest = -(1-0.5*dampingFactor*dt)*uprev(i,j) + 2*u(i,j);

			u_next(i,j) = world(i,j) ? 0 : factor*(dtdt_drdr*(ddx + ddy) + ddt_rest);
			max_value = max(max_value,abs(u_next(i,j)));
		}
	}

	u_prev = u_;
	u_ = u_next;
}

void WaveSolver::RenderWall(int i,int j) {
	double z_top = 0.2;
	double z_bot = 0.0;

	glColor4f(0.8, 0.8, 0.8, 1.0);

	glBegin(GL_POLYGON);
	glNormal3f(0,0,1);
	glVertex3f( x(i)     ,  y(j)   , z_top );
	glVertex3f( x(i+1)   ,  y(j)   , z_top );
	glVertex3f( x(i)     ,  y(j+1) , z_top );
	glVertex3f( x(i+1)   ,  y(j+1) , z_top );
	glEnd();

	glBegin(GL_POLYGON);
	glNormal3f(0,0,-1);
	glVertex3f( x(i)     ,  y(j)     , z_bot );
	glVertex3f( x(i+1)   ,  y(j)     , z_bot );
	glVertex3f( x(i)     ,  y(j+1)   , z_bot );
	glVertex3f( x(i+1)   ,  y(j+1)   , z_bot );
	glEnd();

	glBegin(GL_POLYGON);
	glNormal3f(1,0,0);
	glVertex3f( x(i+1)     ,  y(j)     , z_bot );
	glVertex3f( x(i+1)     ,  y(j)     , z_top );
	glVertex3f( x(i+1)     ,  y(j+1)   , z_bot );
	glVertex3f( x(i+1)     ,  y(j+1)   , z_top );
	glEnd();

	glBegin(GL_POLYGON);
	glNormal3f(-1,0,0);
	glVertex3f( x(i)     ,  y(j)     , z_bot );
	glVertex3f( x(i)     ,  y(j)     , z_top );
	glVertex3f( x(i)     ,  y(j+1)   , z_bot );
	glVertex3f( x(i)     ,  y(j+1)   , z_top );
	glEnd();

	glBegin(GL_POLYGON);
	glNormal3f(0,1,0);
	glVertex3f( x(i+1)  ,  y(j+1)   , z_bot );
	glVertex3f( x(i)    ,  y(j+1)   , z_top );
	glVertex3f( x(i+1)  ,  y(j+1)   , z_bot );
	glVertex3f( x(i)    ,  y(j+1)   , z_top );
	glEnd();

	glBegin(GL_POLYGON);
	glNormal3f(0,-1,0);
	glVertex3f( x(i+1)  ,  y(j)   , z_bot );
	glVertex3f( x(i)    ,  y(j)   , z_top );
	glVertex3f( x(i+1)  ,  y(j)   , z_bot );
	glVertex3f( x(i)    ,  y(j)   , z_top );
	glEnd();
}

// Fonts
void* glutFonts[7] = { 
    GLUT_BITMAP_9_BY_15, 
    GLUT_BITMAP_8_BY_13, 
    GLUT_BITMAP_TIMES_ROMAN_10, 
    GLUT_BITMAP_TIMES_ROMAN_24, 
    GLUT_BITMAP_HELVETICA_10, 
    GLUT_BITMAP_HELVETICA_12, 
    GLUT_BITMAP_HELVETICA_18 
}; 

void glutPrint(float x, float y, void* font, char* text, float r, float g, float b, float a) 
{ 
    if(!text || !strlen(text)) return; 
    bool blending = false; 
    if(glIsEnabled(GL_BLEND)) blending = true; 
    glEnable(GL_BLEND); 
    glColor4f(r,g,b,a); 
    glRasterPos2f(x,y); 
    while (*text) { 
        glutBitmapCharacter(font, *text); 
        text++; 
    } 
    if(!blending) glDisable(GL_BLEND); 
}

// Divide each square (i,j), (i+1,j), (i,j+1), (i+1,j+1) into two triangles and draw them
void WaveSolver::Render() {
	glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    double r,g,b;
    
    for(int i=0;i<Nr-1;i++) {
        for(int j=0;j<Nr-1;j++) {
            if(world(i,j)) glColor4f(1, 1, 1, 1.0);
            else {
	        	r = (u_(i,j) > 0) ? u_(i,j)/max_value : 0;
	        	g = 1.0-abs(u_(i,j))/max_value;
	        	b = (u_(i,j) < 0) ? -u_(i,j)/max_value : 0;
	        	
	            glColor4f(r, g, b, 1.0);
	        }
            


            // The first triangle
            // ***
            // **
            // *
            glBegin(GL_TRIANGLES);

            glVertex3f( x(i)   ,  y(j)   , u_(i,j)   );
            glVertex3f( x(i+1) ,  y(j)   , u_(i+1,j) );
            glVertex3f( x(i)   ,  y(j+1) , u_(i,j+1) );

            glEnd();

            // The other triangle
            //   *
            //  **
            // ***
            glBegin(GL_TRIANGLES);

            glVertex3f( x(i+1)   ,  y(j+1)   , u_(i+1,j+1)   );
            glVertex3f( x(i+1)   ,  y(j)     , u_(i+1,j) );
            glVertex3f( x(i)     ,  y(j+1)   , u_(i,j+1) );

            glEnd();
            
            if(render_wall && world(i,j)) RenderWall(i,j);
        }
    }

    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    
    char *str = new char[50];
    sprintf(str,"Maximum amplitude = %.3f",max_value);
    glutPrint(-4.0,5.5,glutFonts[5],str,1,1,1,1);

    sprintf(str,"dt                          = %.3f",dt);
    glutPrint(-3.675,4.8,glutFonts[5],str,1,1,1,1);

    sprintf(str,"Grid size                = %d",Nr);
    glutPrint(-3.370,4.15,glutFonts[5],str,1,1,1,1);

    sprintf(str,"damping                 = %.3f",dampingFactor);
    glutPrint(-3.13,3.65,glutFonts[5],str,1,1,1,1);
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

		double avg = (r+g+b)/3.0/255.0; // If we have black/white only, the average is 0 (black) to 255 (white)
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