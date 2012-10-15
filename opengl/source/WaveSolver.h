#pragma once

class CWave;
// class CIniFile;

#include <Armadillo>
#include <COpenGL.h>
#include <CIniFile.h>
#include <AObject.h>
#include <Perlin.h>

using namespace arma;

mat readBMP(char* filename);

class WaveSolver {
public:
	int Nr;
	mat ground;
	mat walls;
	mat u_next;
	mat u_;
	mat u_prev;
	mat source;

	vec x;
	vec y;
	double dampingFactor;
	double t;
	double dt;
	double dr;
	double dtdt_drdr;
	double r_min;
	double r_max;
	double max_value;
	double avg_u;
	int mx,my;
	double Lx,Ly,omega;
	
	WaveSolver() {}
	WaveSolver(CIniFile &ini);
	void step();
	void changeGroundZIncrease(double delta_z);
	void calculateWalls();
	void calculateMean();
	void createRandomGauss();
	void calculateSource();
	void copyToGrid(AObject& grid);
	double u(int i,int j, int di=0, int dj=0);
	double uprev(int i,int j, int di=0, int dj=0);
	double calcC(int i, int j);
	int idx(int i);


};
