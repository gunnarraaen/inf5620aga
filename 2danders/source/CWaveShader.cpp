#include <CWaveShader.h>

void CWaveShader::Initialize(string s) {

  ((CShaderParent*)this)->Initialize(s);

}


void CWaveShader::Start() {
  Shader->begin();

  Shader->sendUniform3f((char*)"lightpos",lightpos.x, lightpos.y, lightpos.z);

}

void CWaveShader::End() {
  Shader->disable_multitextures();
  Shader->end();
}

CWaveShader::CWaveShader() {
  vert = string(
  	"uniform vec3 lightpos; \n"
	"varying vec3 normal; \n"
	"varying vec3 myPos; \n"
	"void main(void) \n"
	"{ \n"
	"	normal = gl_Normal; \n"
	"	myPos = gl_Vertex.xyz; \n"
	"   gl_TexCoord[0] = gl_MultiTexCoord0;\n"
	"   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
	"}\n");

 
  frag = string(
  	"uniform vec3 lightpos; \n"
	"varying vec3 normal; \n"
	"varying vec3 myPos; \n"
	"void main(void)\n"
	"{\n "
	"  vec4 val = vec4(0.2,myPos.z*10.0,1,1);"
	"  float light = clamp(dot(normalize(lightpos), normal), 0.0, 1.0);"
	"  gl_FragColor = val*light; \n"
	"  gl_FragColor.w = 0.7;"
	"}\n");




}



void CGroundShader::Initialize(string s) {

  ((CShaderParent*)this)->Initialize(s);

}


void CGroundShader::Start() {
  Shader->begin();

  Shader->sendUniform3f((char*)"lightpos",lightpos.x, lightpos.y, lightpos.z);

}

void CGroundShader::End() {
  Shader->disable_multitextures();
  Shader->end();
}

CGroundShader::CGroundShader() {
  vert = string(
  	"uniform vec3 lightpos; \n"
	"varying vec3 normal; \n"
	"varying vec3 myPos; \n"
	"void main(void) \n"
	"{ \n"
	"	normal = gl_Normal; \n"
	"	myPos = gl_Vertex.xyz; \n"
	"   gl_TexCoord[0] = gl_MultiTexCoord0;\n"
	"   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
	"}\n");

 
  frag = string(
  	"uniform vec3 lightpos; \n"
	"varying vec3 normal; \n"
	"varying vec3 myPos; \n"
	"void main(void)\n"
	"{\n "
	"  vec4 val = vec4(0.7,0.5,0.3,1);"
	"  float light = clamp(dot(normalize(lightpos), normal), 0.0, 1.0);"
	"  gl_FragColor = val*light; \n"
	"  gl_FragColor.w = 1.0;"
	"}\n");




}