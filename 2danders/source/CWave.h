#pragma once

#include <CApplication.h>
#include <CDots.h>
#include <CPoints.h>
#include <WaveSolver.h>

using namespace std;

class Variables  {
public:
    WaveSolver solver;
    double time;
    int speed;
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
};
