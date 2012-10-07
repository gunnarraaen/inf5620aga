#pragma once

#include <CApplication.h>

class CPoint {
public:
	CVector pos;
	CVector* camera;

	CPoint(CVector p, CVector* camOA);
	void Render();
	void Update();
}; 
