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
    var.solver = WaveSolver(ini);
    var.time = 0;
    var.speed = 3;  
    var.waveGrid.InitializeGrid(ini.getint("grid_size"),2.0);
    var.groundGrid.InitializeGrid(ini.getint("grid_size"),2.0);
    var.groundGrid.copyGridFromBMP(var.solver.H);
    var.groundGrid.pos.z = -0.5;


    var.waveShader.Initialize("waveShader");
    var.groundShader.Initialize("groundShader");

    Initialized = true;
}  


void CWave::renderWave() {

    var.solver.copyToGrid(var.waveGrid);
    var.waveGrid.calculateGridFaceNormals();

    float A = 5;
    float B = 2;

    var.waveShader.lightpos.Set(A*cos(var.solver.time*B),A*sin(var.solver.time*B),1);
    var.waveShader.Start();

    glColor4f(0,1,0,1);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    var.waveGrid.RenderTriangles();

    var.waveShader.End();
}

void CWave::renderGround() {

    var.groundShader.lightpos  = var.waveShader.lightpos;
    var.groundShader.Start();

    glColor4f(0,1,0,1);
    glDisable(GL_BLEND);
    var.groundGrid.RenderTriangles();

    var.groundShader.End();
}

void CWave::Display (void) {
    if (!Initialized)
    return;

    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT); 

    ogl.setperspective(60);
    
    double x = 2*cos(var.solver.time);
    double y = 2*sin(var.solver.time);

     ogl.camera = CVector(0.0,-2.0,1.3);
    //ogl.camera = CVector(x,y,1.3);
    
    ogl.ypr = CVector(0,0,90); // Rotate camera

    ogl.target = CVector(0,0,0);
    ogl.setup_camera();

    // var.solver.Render();
    renderGround();
    renderWave();    
    
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
    if (key=='g')
        var.solver.render_ground = !var.solver.render_ground;
    if (key=='w')
        var.solver.render_wave = !var.solver.render_wave;
    if (key=='a')
        var.solver.render_wall = !var.solver.render_wall;

    key = '0';

}
