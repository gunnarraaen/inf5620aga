#pragma once
#include <CVector.h>
#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;

class AVertex {
public:
	CVector pos,normal,color;
};

class AFace {
public:
	int face[3];
	CVector normal;

	void calculateNormal(vector<AVertex>& v) {
		CVector& p1 = v[face[0]].pos;
		CVector& p2 = v[face[1]].pos;
		CVector& p3 = v[face[2]].pos;
		normal = (p3-p1).Cross(p2-p1).Normalize();
	}

};


class AObject: public AVertex {
public:
	vector<AVertex> vertices;
	vector<AFace>     faces;

	int       		gridSize;
	bool			hasNormals, hasColors;


	void InitializeGrid(int _gridSize,double _realSize);
	void copyGridFromBMP(mat);
	void RenderTriangles();
	void calculateGridFaceNormals();
	void calculateGridVertexNormals();
	AVertex* getGridPos(const int i, const int j);

};