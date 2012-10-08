#pragma once

class CWave;
// class CIniFile;

#include <Armadillo>
#include <COpenGL.h>
#include <CIniFile.h>
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
	double dampingFactor;
	double time;
	double dt;
	double dr;
	double dtdt_drdr;
	double r_min;
	double r_max;
	double max_value;
	bool   render_wall;
	bool   render_ground;
	bool   render_wave;

	WaveSolver() {}
	WaveSolver(CIniFile &ini);
	void step();
	void Render();
	void RenderWave(int i,int j);
	void RenderGround(int i,int j);
	void RenderWall(int i, int j, int di,int dj);

	double u(int i,int j, int di=0, int dj=0);
	double uprev(int i,int j, int di=0, int dj=0);
	double calcC(int i, int j);
	int idx(int i);


};
