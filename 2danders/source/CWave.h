#pragma once

#include <WaveSolver.h>
#include <CApplication.h>
#include <CDots.h>
#include <CPoints.h>
#include <AObject.h>
#include <CWaveShader.h>

using namespace std;

class Variables  {
public:
    WaveSolver solver;
    AObject    waveGrid;
    AObject    groundGrid;
    int speed;
    double theta;
    CWaveShader waveShader;
    CGroundShader groundShader;
    bool render_ground;
    bool render_wall;
    bool render_wave;

};

class CWave : public CApplication {
public:
    static Variables var;

    static void Initialize_();

    // Rasterizing
    static void Display(void);

    // Update physics / OpenGL
    static void Update(void);

    // Handle events (keys, mouse)
    static void Events();

    // internal stuff
private:    
    static void renderWave();
    static void renderGround();
    static void renderWalls();
    static void RenderWall(int i,int j, int di, int dj);

};
