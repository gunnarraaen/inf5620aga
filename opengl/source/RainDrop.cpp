#include <RainDrop.h>

RainDrop::RainDrop(int _i, int _j, double _z) {
	i = _i;
	j = _j;
	z = _z;
}

void RainDrop::render(vec &x, vec &y) {
	glVertex3f(x(i), y(j), z);
}