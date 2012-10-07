#include <CPoint.h>

CPoint::CPoint(CVector p,CVector* cam) {
    pos = p;
    camera = cam;
}

void CPoint::Update() {

}

void CPoint::Render() {
	glVertex3f(pos.x, pos.y, pos.z);
}
