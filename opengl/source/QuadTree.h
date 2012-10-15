#include <CVector.h>
#include <armadillo>

using namespace arma;


class QuadTree {
public:
	int x0,y0,size, lod;
	QuadTree *NW;
	QuadTree *SW;
	QuadTree *SE;
	QuadTree *NE;

	public static const lodmax = 10;

	QuadTree();
	Setup(int _x0, int _y0, int _size, int _lod);

	Render(int depth_max);
};