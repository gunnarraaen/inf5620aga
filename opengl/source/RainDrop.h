#include <CVector.h>
#include <armadillo>
#include <COpenGL.h>
using namespace arma;

class RainDrop {
public:
	int i;
	int j;
	double z;
	RainDrop(int _i, int _j, double _z);
	void render(vec &x, vec &y);
};