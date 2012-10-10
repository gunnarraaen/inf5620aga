#include <QuadTree.h>

QuadTree::QuadTree(int _x0, int _y0, int _size) {
	x0 = _x0;
	y0 = _y0;
	size = _size;

	xCenter = x0 + 0.5*size;
	yCenter = y0 + 0.5*size;

	if(0.5*size > 0) {
		NW = QuadTree(x0,y0,size/2);
		SW = QuadTree(x0,y0+size/2,size/2);
		SE = QuadTree(x0+size/2,y0+size/2,size/2);
		NW = QuadTree(x0+size/2,y0,size/2);
	}


}