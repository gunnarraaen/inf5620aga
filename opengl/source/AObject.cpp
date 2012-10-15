#include <AObject.h>
#include <COpenGL.h>

void AObject::InitializeGrid(int _gridSize,double _realSize) {
	gridSize = _gridSize;
	int N = gridSize;

	vertices.resize(N*N);
	faces.resize((N-1)*(N-1)*2);

	for (int i=0;i<N;i++) {
		for (int j=0;j<N;j++) {
			AVertex v;

			double x = ((i-N/2.0)/(double)N)*_realSize;
			double y = ((j-N/2.0)/(double)N)*_realSize;

			v.pos.Set(x,y,0);
			vertices[i+j*N] = v;
		}
	}

	int face = 0;
	for (int i=0;i<(N-1);i++) {
		for (int j=0;j<(N-1);j++) {
			int c = i + j*N;

			faces[face].face[0] = c;
			faces[face].face[1] = c+1;
			faces[face].face[2] = c+N+1;
			
			vertices[c].faces.push_back(face);
			vertices[c+1].faces.push_back(face);
			vertices[c+N+1].faces.push_back(face);
			
			face++;

			faces[face].face[0] = c;
			faces[face].face[1] = c+N+1;
			faces[face].face[2] = c+N;
			vertices[c].faces.push_back(face);
			vertices[c+N+1].faces.push_back(face);
			vertices[c+N].faces.push_back(face);

			face++;


		}		
	}

	calculateGridFaceNormals();
}




void AVertex::calculateNormal(vector<AFace>& faceList) {
	normal.Set(0,0,0);
	
	int size=faces.size();

	for(int i=0;i<size;i++) {
		normal = normal + faceList[faces[i]].normal;
	}

	normal = normal / ((size>0) ? size : 1);
}

void AObject::calculateGridVertexNormals() {
	int N = gridSize;
	for(int i=0;i<N*N;i++) 
		vertices[i].calculateNormal(faces);
}

void AObject::calculateGridFaceNormals() {
	for (int i=0;i<faces.size();i++)
		faces[i].calculateNormal(vertices);
}


void AObject::RenderTriangles() {
	glPushMatrix();
	glTranslatef(pos.x,pos.y,pos.z);

	glBegin(GL_TRIANGLES);
	for (int i=0;i<faces.size();i++) {
		for (int j=0;j<3;j++) {
			AVertex *v = &vertices[ faces [i].face[j]];
			glNormal3f( v->normal.x, v->normal.y, v->normal.z);
			glVertex3f( v->pos.x, v->pos.y, v->pos.z);
		}
	}
	glEnd();
	glPopMatrix();
}

void AObject::copyGridFromBMP(mat bmp) {
	for(int i=0;i<gridSize;i++) {
		for(int j=0;j<gridSize;j++) {
			getGridPos(i,j)->pos.z = bmp(i,j);
		}
	}

	calculateGridFaceNormals();
	calculateGridVertexNormals();
}


AVertex* AObject::getGridPos(const int i, const int j) {
	if (i>gridSize && i<0) return 0;
	if (j>gridSize && j<0) return 0;
	return &vertices[i + j*gridSize];
}

void AObject::renderAsQuad(double size, int lod, AQuadPoint qp, CVector camera) {

	qp.center.Set(0,0,0);

	for (int i=0;i<qp.points.size();i++)
		qp.center = qp.center + qp.points[i];

	qp.center = qp.center/(double)qp.points.size();
//	cout << lod << " " << qp.center;

	double len = (getGridPos(qp.center.x,qp.center.y)->pos - camera).Length();
	len = pow(len, 1.0);
	if (len*lod<size && lod<=8) {
	//if (lod<=8) {
	AQuadPoint q1;
		q1.points[0] = qp.points[0];
		q1.points[1] = CVector(qp.center.x, qp.points[0].y,0);
		q1.points[2] = qp.center;
		q1.points[3] = CVector(qp.points[0].x, qp.center.y,0);
		renderAsQuad(size, lod+1, q1, camera);

		AQuadPoint q2;
		q2.points[0] = CVector(qp.center.x, qp.points[1].y,0);
		q2.points[1] = qp.points[1];
		q2.points[2] = CVector(qp.points[1].x, qp.center.y,0);
		q2.points[3] = qp.center;
		renderAsQuad(size, lod+1, q2, camera);

		AQuadPoint q3;
		q3.points[0] = qp.center;
		q3.points[1] = CVector(qp.points[1].x, qp.center.y,0);
		q3.points[2] = qp.points[2];
		q3.points[3] = CVector(qp.center.x, qp.points[2].y,0);
		renderAsQuad(size, lod+1, q3, camera);

		AQuadPoint q4;
		q4.points[0] = CVector(qp.points[0].x, qp.center.y,0);
		q4.points[1] = qp.center;
		q4.points[2] = CVector(qp.center.x, qp.points[2].y,0);
		q4.points[3] = CVector(qp.points[3].x, qp.points[3].y,0);
		renderAsQuad(size, lod+1, q4, camera);

		return;
	}


	glBegin(GL_TRIANGLE_FAN);

	AVertex* p = getGridPos(qp.center.x,qp.center.y);		

	p->renderVertex();

	for (int k=0;k<qp.points.size()+1;k++) {
		int kk = k % qp.points.size();

		int i = qp.points[kk].x;
		int j = qp.points[kk].y;
		if (i==gridSize) i--;
		if (j==gridSize) j--;
		
		AVertex* p = getGridPos(i,j);
		if (p!=0) {
			p->renderVertex();
		}

	}
	glEnd();
	
}


