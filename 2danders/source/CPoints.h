#pragma once

#include <CPoint.h>
#include <vector>
#include <COpenGL.h>
#include <CApplication.h>

using namespace std;

class CPoints {
public:
	vector<CPoint> points;
	static CVector* camera;

	void Initialize(int, COpenGL& gl);
	void Render();
	void Update();
};
