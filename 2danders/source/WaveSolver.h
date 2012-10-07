#pragma once

class CWave;

#include <Armadillo>
#include <COpenGL.h>
using namespace arma;

mat readBMP(char* filename);

class WaveSolver {
public:
	int Nr;
	mat H;
	mat world;
	mat u_next;
	mat u_;
	mat u_prev;
	vec x;
	vec y;
	double time;
	double dt;
	double dr;
	double dtdt_drdr;
	double r_min;
	double r_max;
	double max_value;

	WaveSolver();
	void step();
	void Render();
	double u(int i,int j, int di=0, int dj=0);
	double uprev(int i,int j, int di=0, int dj=0);
	double calcC(int i, int j);
	int idx(int i);


};