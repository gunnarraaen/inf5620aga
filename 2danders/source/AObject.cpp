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

/*			if (hasColors) glColor3f( v->color.x, v->color.y, v->color.z);		
			if (hasNormals) glNormal3f( v->normal.x, v->normal.y, v->normal.z);	*/	
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