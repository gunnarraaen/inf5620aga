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
    var.speed = 2;
    var.theta = 3*M_PI/2;
    var.cam_radius  = 2.0;
    var.waveGrid.InitializeGrid(ini.getint("grid_size"),2.0);
    var.groundGrid.InitializeGrid(ini.getint("grid_size"),2.0);
    var.groundGrid.copyGridFromBMP(var.solver.ground);
    var.groundGrid.pos.z = -1;
    var.render_ground = true;
    var.render_wall = true;
    var.render_wave = true;

    var.waveShader.Initialize("waveShader");
    var.groundShader.Initialize("groundShader");

    Initialized = true;
}  

void CWave::RenderWall(int i,int j, int di, int dj) {
    vec &x = var.solver.x;
    vec &y = var.solver.x;
    int Nr = var.solver.Nr;

    double z_top = 0.05;
    double z_bot = -1.0;

    glColor4f(0.05, 0.05, 0.05, 1.0);

    glBegin(GL_POLYGON);
    
    glVertex3f( x(0)      ,  y(0)    , z_bot );
    glVertex3f( x(Nr-1)   ,  y(0)    , z_bot );
    glVertex3f( x(Nr-1)   ,  y(0)    , z_top );
    glVertex3f( x(0)      ,  y(0)    , z_top );
    glEnd();

    glBegin(GL_POLYGON);
    
    glVertex3f( x(0)      ,  y(0)       , z_bot );
    glVertex3f( x(0)      ,  y(Nr-1)    , z_bot );
    glVertex3f( x(0)      ,  y(Nr-1)    , z_top );
    glVertex3f( x(0)      ,  y(0)       , z_top );
    glEnd();

    glBegin(GL_POLYGON);
    
    glVertex3f( x(0)      ,  y(Nr-1)    , z_bot );
    glVertex3f( x(Nr-1)   ,  y(Nr-1)    , z_bot );
    glVertex3f( x(Nr-1)   ,  y(Nr-1)    , z_top );
    glVertex3f( x(0)      ,  y(Nr-1)    , z_top );
    glEnd();

    glBegin(GL_POLYGON);
    
    glVertex3f( x(Nr-1)      ,  y(Nr-1)    , z_bot );
    glVertex3f( x(Nr-1)      ,  y(0)       , z_bot );
    glVertex3f( x(Nr-1)      ,  y(0)       , z_top );
    glVertex3f( x(Nr-1)      ,  y(Nr-1)    , z_top );
    
    glEnd();
}

void CWave::renderWalls() {
    glDisable(GL_BLEND);
    
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    int Nr = var.solver.Nr;
    RenderWall(0,0,Nr-1,1);
    RenderWall(0,0,1,Nr-1);
    RenderWall(Nr-1,0,-1,Nr-1);
    RenderWall(0,Nr-1,Nr-1,-1);
}

void CWave::renderWave() {

    var.solver.copyToGrid(var.waveGrid);
    var.waveGrid.calculateGridFaceNormals();
    var.waveGrid.calculateGridVertexNormals();

    double w = 0.5;
    var.waveShader.lightpos.Set(2*cos(var.solver.time*w),2*sin(var.solver.time*w),1);
    
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
    
    double x = var.cam_radius*cos(var.theta);
    double y = var.cam_radius*sin(var.theta);
    
    // ogl.camera = CVector(0.0,-2.0,1.3);
    ogl.camera = CVector(x,y,1.3);
    
    ogl.ypr = CVector(0,0,90); // Rotate camera

    ogl.target = CVector(0,0,0);
    ogl.setup_camera();

    // var.solver.Render();
    if(var.render_wall)   renderWalls();
    if(var.render_ground) renderGround();
    if(var.render_wave)   renderWave();
    
    
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
        var.speed = 6;
    if (key=='6')
        var.speed = 8;
    if (key=='7')
        var.speed = 14;
    if (key=='8')
        var.speed = 20;
    if (key=='9')
        var.speed = 30;
    if (key=='q')
        var.render_ground = !var.render_ground;
    if (key=='w')
        var.render_wave = !var.render_wave;
    if (key=='e')
        var.render_wall = !var.render_wall;
    if (key=='a')
        var.theta -= 0.05;
    if (key=='d')
        var.theta += 0.05;
    if (key=='+')
        var.cam_radius -= 0.05;
    if (key=='-')
        var.cam_radius += 0.05;

    key = '0';

}
