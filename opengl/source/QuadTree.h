#include <CVector.h>
#include <armadillo>

using namespace arma;


class QuadTree {
public:
	int x0,y0,size;
	QuadTree *NW;
	QuadTree *SW;
	QuadTree *SE;
	QuadTree *NE;

	QuadTree();
	Render(int depth_max);
};