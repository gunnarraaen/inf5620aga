#include <RainDrop.h>
#include <GLUT/glut.h>
RainDrop::RainDrop(int _i, int _j, double _z) {
	i = _i;
	j = _j;
	z = _z;
}

void RainDrop::render(vec &x, vec &y) {
	glPushMatrix();
	glTranslatef(x(i),y(j),z);
	glutSolidSphere(0.01, 4, 4);
	glPopMatrix();
}