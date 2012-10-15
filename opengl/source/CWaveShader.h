#pragma once

#include <CShaders.h>
#include <armadillo>

using namespace arma;

class CWaveShader : public CShaderParent {
 public:
  CWaveShader();  

  CVector lightpos; 
  CVector targetdir; 

  void Initialize(string);
  void Start();
  void End();

};

class CGroundShader : public CShaderParent {
 public:
  CGroundShader();  

  CVector lightpos; 
  CVector targetdir; 

  void Initialize(string);
  void Start();
  void End();

};
