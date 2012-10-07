#include <CWave.h>
#include <vector>

Variables CWave::var;

void CWave::Update(void) {
    Initialize_();
    glutPostRedisplay();
    InternalUpdate();
    Events();

    var.solver.step();
    var.solver.step();
    var.solver.step();
}  

void CWave::Initialize_() {
    if (Initialized)
    return;
    
    var.solver = WaveSolver();

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

    ogl.setperspective(60);
    double x = cos(var.solver.time);
    double y = sin(var.solver.time);

    ogl.camera = CVector(0.0,-2.0,1.3);
    // ogl.camera = CVector(x,y,1.3);
    ogl.ypr = CVector(0,0,90); // Rotate camera

    ogl.target = CVector(0,0,0);
    ogl.setup_camera();

    // var.points.Render();
    var.solver.Render();
    // RenderAxis();

    glutSwapBuffers(); 
    Events();
}


void CWave::Events ()  {
    bool go = true;
    if (key==27)
        exit(1);
}


