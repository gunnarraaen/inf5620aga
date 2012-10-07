#pragma once

mat readBMP(char* filename);

class WaveSolver {
public:
	static int Nr;
	static int Nt;
	static mat H;
	static mat world;
	static mat u_next;
	static mat u_;
	static mat u_prev;
	static vec x;
	static vec y;
	static bool height_from_file;
	static bool world_from_file;
	static bool save_to_file;
	static double time;
	static double dt;
	static double dr;
	static double dtdt_drdr;
	static double r_min;
	static double r_max;

	WaveSolver();
	void step();

};
