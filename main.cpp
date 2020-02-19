#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include "STF.h"
using namespace std;

STF stf;

// rendering projection parameters
const static int WINDOW_WIDTH = 800;
const static int WINDOW_HEIGHT = 600;
const static double VIEW_WIDTH = 1.5*800.f;
const static double VIEW_HEIGHT = 1.5*600.f;

const static float H = 16.f; // kernel radius
const static float DT = 0.4f; // integration timestep

void InitSTF(void) {

  // create solid
  stf.objects.solids.clear();
  stf.objects.solids.resize(1);

  // create nodes
  int Nx = 4;
  int Ny = 12;
  int dx = 40;
  int dy = 40;
  stf.objects.solids[0].nodes.resize(Nx*Ny);
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      stf.objects.solids[0].nodes[Nx*i+j].position[0] = 0.5*VIEW_WIDTH+dx*j-0*i;
      stf.objects.solids[0].nodes[Nx*i+j].position[1] = 0.5*VIEW_HEIGHT+dy*i+0*j;
    }
  }

  // create elements
  stf.objects.solids[0].elements.resize((Nx-1)*(Ny-1));
  for (int i = 0; i < (Ny-1); i++) {
    for (int j = 0; j < (Nx-1); j++) {
      stf.objects.solids[0].elements[(Nx-1)*i+j] = new Quadrilateral();
      stf.objects.solids[0].elements[(Nx-1)*i+j]->nodeIDs = { Nx*i+j, Nx*i+j+1, Nx*(i+1)+j+1, Nx*(i+1)+j };
    }
  }

  // create material
  stf.objects.solids[0].material.rho = 0.001;
  stf.objects.solids[0].material.E = 50.0;
  stf.objects.solids[0].material.nu = 0.2;
  stf.objects.solids[0].material.dVdT =-0.5;
  stf.objects.solids[0].material.alpha = 0.001;
  stf.objects.solids[0].material.fiber[0] = 0.0;
  stf.objects.solids[0].material.fiber[1] = 1.0;

  // create body forces
  stf.environment.bodyForces.resize(2);

  // create gravity force
  Gravity* gravityForce = new Gravity();
  gravityForce->gravity = 1.0;
  stf.environment.bodyForces[0] = gravityForce;

  // create bouyant force
  BuoyantForce* bouyantForce = new BuoyantForce();
  bouyantForce->gravity = 1.0;
  bouyantForce->rho_fluid = 0.002f;
  bouyantForce->y_surface = 0.0*VIEW_HEIGHT;
  stf.environment.bodyForces[1] = bouyantForce;

  // create damping force
  stf.environment.dampingForces.resize(1);
  MassDamping* damping = new MassDamping();
  damping->viscosity = 0.001f;
  damping->y_surface = 0.0*VIEW_HEIGHT;
  stf.environment.dampingForces[0] = damping;
  
  // create boundary
  stf.environment.boundaries.resize(3);
  stf.environment.boundaries[0] = new Boundary();
  stf.environment.boundaries[0]->basePoint[0] = 0;
  stf.environment.boundaries[0]->basePoint[1] = 0;
  stf.environment.boundaries[0]->normal[0] = 0;
  stf.environment.boundaries[0]->normal[1] = 1;
  stf.environment.boundaries[1] = new Boundary();
  stf.environment.boundaries[1]->basePoint[0] = 0;
  stf.environment.boundaries[1]->basePoint[1] = 0;
  stf.environment.boundaries[1]->normal[0] = 1;
  stf.environment.boundaries[1]->normal[1] = 0;
  stf.environment.boundaries[2] = new Boundary();
  stf.environment.boundaries[2]->basePoint[0] = VIEW_WIDTH;
  stf.environment.boundaries[2]->basePoint[1] = 0;
  stf.environment.boundaries[2]->normal[0] =-1;
  stf.environment.boundaries[2]->normal[1] = 0;

  // initialize
  stf.initialize();

}

void Update(void) {
  int Nsubincrements = 100;
  float ddt = DT / Nsubincrements;
  for (int i = 0; i < Nsubincrements; i++) {
    stf.timeIntegrate(ddt);
  }

  glutPostRedisplay();
}

void InitGL(void) {
  glClearColor(0.9f,0.9f,0.9f,1);
  glEnable(GL_POINT_SMOOTH);
  glPointSize(H/2.f);
  glMatrixMode(GL_PROJECTION);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void Render(void) {
  glClear(GL_COLOR_BUFFER_BIT);
	
  glLoadIdentity();
  glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1);

  glBegin(GL_QUADS);
  for (auto &e : stf.objects.solids[0].elements) {
    for (auto &i : e->nodeIDs) {
      if (e->active) {
	float fresh = stf.objects.solids[0].nodes[i].temperature;
	float charred = 1.0 - stf.objects.solids[0].nodes[i].temperature;
	glColor4f(0.7f*charred+0.2*fresh, 0.1f*charred+0.2*fresh, 0.1f*charred+0.2*fresh, 1);
	glVertex2f(stf.objects.solids[0].nodes[i].position[0],
		   stf.objects.solids[0].nodes[i].position[1]);
      }
    }
  }
  glEnd();

  glBegin(GL_QUADS);
  glColor4f(0.5f, 0.4f, 0.0f, 0.5f);
  glVertex2f(0,0);
  glVertex2f(VIEW_WIDTH,0);
  glVertex2f(VIEW_WIDTH,0.0*VIEW_HEIGHT);
  glVertex2f(0,0.0*VIEW_HEIGHT);
  glEnd();

  /*
  glColor4f(0.25f, 0.1f, 0.1f, 1);
  glBegin(GL_POINTS);
  for (auto &p : stf.objects.solids[0].nodes) {
    glVertex2f(p.position[0], p.position[1]);
  }
  glEnd();
  */

  glutSwapBuffers();
}

void Keyboard(unsigned char c, __attribute__((unused)) int x, __attribute__((unused)) int y) {   
  switch(c) {
  case 'r': 
  case 'R':  
    InitSTF(); 
    break;
  }
}

void Arrows(int key, __attribute__((unused)) int x, __attribute__((unused)) int y) {   
  switch(key) {
  case GLUT_KEY_LEFT:  
    stf.objects.solids[0].material.Et = -0.5; 
    break;
  case GLUT_KEY_RIGHT: 
    stf.objects.solids[0].material.Et = -0.5; 
    break;
  case GLUT_KEY_UP:    
    stf.objects.solids[0].material.Em = -0.5; 
    break;
  case GLUT_KEY_DOWN: 
    stf.objects.solids[0].material.Em = -0.5; 
    break;
  }
}

int main(int argc, char** argv) {
  glutInitWindowSize(WINDOW_WIDTH,WINDOW_HEIGHT);
  glutInit(&argc, argv);
  glutCreateWindow("STF Library Demo");
  glutDisplayFunc(Render);
  glutIdleFunc(Update);
  glutKeyboardFunc(Keyboard);
  glutSpecialFunc(Arrows);

  InitGL();
  InitSTF();

  glutMainLoop();
  return 0;
}
