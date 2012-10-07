#include <CPoints.h>

bool sortFunc(CPoint a, CPoint b) {
    return ((a.pos-*a.camera).Length()>(b.pos-*b.camera).Length());
}  

void CPoints::Initialize(int N, COpenGL &ogl) {
    double x,y,z;
    for(int i=0;i<N;i++) {
        
        x = -0.5+i/((double)N);
        for(int j=0;j<N;j++) {
            y = -0.5+j/((double)N);

            z = 0.3*exp(-(x*x+y*y)*30);

            CVector pos = CVector(x,y,z);
            CPoint p = CPoint(pos,&ogl.camera);
            points.push_back(p);
        }
    }
}

void CPoints::Update() {
    // for (int i=0;i<points.size();i++)
        // points[i].Update();

    // sort(points.begin(), points.end(), sortFunc);
}

inline int idx(int i,int j,int N) {
    return i+j*N;
}

void CPoints::Render() {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glColor4f(1.0, 1.0,1.0, 1.0);
    CVector v1;
    CVector v2;
    CVector v3;
    
    int len = sqrt(points.size());
    glColor4f(1.0, 1,1, 1.0);

    for(int i=0;i<len-1;i++) {
        for(int j=0;j<len-1;j++) {
            glBegin(GL_TRIANGLES);
            v1 = points[idx(i,j,len)].pos;
            v2 = points[idx(i+1,j,len)].pos;
            v3 = points[idx(i,j+1,len)].pos;
            
            glColor4f(0.2+v1.z, 0,0, 1.0);

            glVertex3f(v1.x,v1.y,v1.z);
            glVertex3f(v2.x,v2.y,v2.z);
            glVertex3f(v3.x,v3.y,v3.z);

            glEnd();

            glBegin(GL_TRIANGLES);
            v1 = points[idx(i+1,j+1,len)].pos;
            v2 = points[idx(i+1,j,len)].pos;
            v3 = points[idx(i,j+1,len)].pos;

            glVertex3f(v1.x,v1.y,v1.z);
            glVertex3f(v2.x,v2.y,v2.z);
            glVertex3f(v3.x,v3.y,v3.z);

            glEnd();
        }
    }

    // glEnd();
    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    
}


