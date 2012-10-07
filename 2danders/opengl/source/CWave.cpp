#include <CWave.h>
#include <vector>

Variables CWave::var;

void CWave::Update(void) {
    Initialize_();
    glutPostRedisplay();
    InternalUpdate();
    Events();

    var.time+=0.002;

    var.points.Update();
}  

void CWave::Initialize_() {
    if (Initialized)
    return;

    var.time = 0;

    var.points.Initialize(200, ogl);
    Initialized = true;
}  

void CWave::RenderAxis() {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glBegin(GL_LINES);
    glColor4f(1.0, 0.0,0.0, 1.0);
    glVertex3f(0,0,0);
    glVertex3f(1,0,0);
    glEnd();

    glBegin(GL_LINES);
    glColor4f(0.0, 1.0,0.0, 1.0);
    glVertex3f(0,0,0);
    glVertex3f(0,1,0);
    glEnd();

    glBegin(GL_LINES);
    glColor4f(0.0, 0.0,1.0, 1.0);
    glVertex3f(0,0,0);
    glVertex3f(0,0,1);
    glEnd();

    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
}

void CWave::Display (void) {
    if (!Initialized)
    return;

    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |GL_STENCIL_BUFFER_BIT); 

    double A = 10;

    ogl.setperspective(75);
    double x = cos(var.time);
    double y = sin(var.time);

    ogl.camera = CVector(x,y,1);
    ogl.ypr = CVector(0,0,90); // Rotate camera

    ogl.target = CVector(0,0,0);
    ogl.setup_camera();

    var.points.Render();
    RenderAxis();

    glutSwapBuffers(); 
    Events();
}


void CWave::Events ()  {
    bool go = true;
    if (key==27)
        exit(1);
}


