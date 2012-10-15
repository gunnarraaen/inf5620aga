#pragma once
#include <CVector.h>
#include <vector>
#include <COpenGL.h>
#include <armadillo>
using namespace std;
using namespace arma;

class AFace;

class AVertex {
public:
	CVector pos,normal,color;
	vector<int> faces;

	void calculateNormal(vector<AFace> &faceList);
	void renderVertex() {
		glNormal3f(normal.x, normal.y, normal.z);
		glVertex3f(pos.x, pos.y, pos.z);
	}
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

class AQuadPoint {
	public:
	vector<CVector> points;
	CVector center;
	AQuadPoint() {
		points.resize(4);
	}

};


class AObject: public AVertex {
public:
	vector<AVertex> vertices;
	vector<AFace>     faces;

	int       		gridSize;
	bool			hasNormals, hasColors;
	mat             vertexFaceList;

	void InitializeGrid(int _gridSize,double _realSize);
	void copyGridFromBMP(mat);
	void RenderTriangles();
	void calculateGridFaceNormals();
	void calculateGridVertexNormals();

	void renderAsQuad(double size, int lod, AQuadPoint qp, CVector camera);
	AVertex* getGridPos(const int i, const int j);

};