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
    double time;
    int speed;
    CWaveShader waveShader;
    CGroundShader groundShader;

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

};
