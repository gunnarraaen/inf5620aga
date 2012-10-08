#include <CWave.h>
#include <vector>
#include <GLUT/glut.h>

Variables CWave::var;

void CWave::Update(void) {
    if(!Initialized) Initialize_();
    glutPostRedisplay();
    InternalUpdate();
    Events();

    for(int i=0;i<var.speed;i++) 
        var.solver.step();
}  

void CWave::Initialize_() {
    var.solver = WaveSolver();
    var.time = 0;
    var.speed = 3;

    Initialized = true;
}  

void CWave::Display (void) {
    if (!Initialized)
    return;

    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT); 

    ogl.setperspective(60);
    
    ogl.camera = CVector(0.0,-2.0,1.3);
    
    ogl.ypr = CVector(0,0,90); // Rotate camera

    ogl.target = CVector(0,0,0);
    ogl.setup_camera();

    var.solver.Render();

    glutSwapBuffers(); 
    Events();
}

void CWave::Events ()  {
    bool go = true;
    if (key==27)
        exit(1);
    if (key=='1')
        var.speed = 1;
    if (key=='2')
        var.speed = 2;
    if (key=='3')
        var.speed = 3;
    if (key=='4')
        var.speed = 4;
    if (key=='5')
        var.speed = 5;
    if (key=='6')
        var.speed = 6;
    if (key=='7')
        var.speed = 7;
    if (key=='8')
        var.speed = 8;
    if (key=='9')
        var.speed = 9;
}
