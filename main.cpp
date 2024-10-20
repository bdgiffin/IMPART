#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include "STF.h"
#include "ColorMap.h"
using namespace std;

// declare global STF object instance
STF stf;

// declare global color maps:
ColorMap rawColor;
ColorMap thermalColor;
int colorPlot = 0;

// rendering projection parameters
const static int WINDOW_WIDTH = 800;
const static int WINDOW_HEIGHT = 600;
const static float VIEW_WIDTH = 1.5*800.f;
const static float VIEW_HEIGHT = 1.5*600.f;

const static float H = 32.f; // kernel radius
const static float DT = 0.5f; // integration timestep
float vis_time = 0.0; // total simuation time

// Define heating settings
const static float min_heat = 0.0;
const static float max_heat = 100.0;
const static float dheat = 20.0;

void InitSTF(void) {

  // define color map(s)
  rawColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.0,{125.0/256.0,10.0/256.0,10.0/256.0}));
  rawColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.2,{175.0/256.0,64.0/256.0,64.0/256.0}));
  rawColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.3,{215.0/256.0,115.0/256.0,115.0/256.0}));
  rawColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.5,{240.0/256.0,240.0/256.0,192.0/256.0}));
  rawColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.6,{156.0/256.0,123.0/256.0,117.0/256.0}));
  rawColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.9,{84.0/256.0,50.0/256.0,43.0/256.0}));
  rawColor.colorTable.push_back(std::pair<float,std::array<float,3> >(1.0,{33.0/256.0,7.0/256.0,1.0/256.0}));

  thermalColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.0,{3.0/256.0,8.0/256.0,71.0/256.0}));
  thermalColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.2,{68.0/256.0,9.0/256.0,130.0/256.0}));
  thermalColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.4,{158.0/256.0,22.0/256.0,79.0/256.0}));
  thermalColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.6,{222.0/256.0,114.0/256.0,51.0/256.0}));
  thermalColor.colorTable.push_back(std::pair<float,std::array<float,3> >(0.8,{252.0/256.0,241.0/256.0,116.0/256.0}));
  thermalColor.colorTable.push_back(std::pair<float,std::array<float,3> >(1.0,{256.0/256.0,256.0/256.0,256.0/256.0}));
  
  // create solid
  stf.objects.solids.clear();
  stf.objects.solids.resize(1);

  // create nodes
  int Nx = 3*2;
  int Ny = 9*2;
  int dx = 60/2;
  int dy = 60/2;
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
  stf.objects.solids[0].material.E = 60.0;
  stf.objects.solids[0].material.nu = 0.1;
  stf.objects.solids[0].material.dVdT =-0.5;
  stf.objects.solids[0].material.alpha = 0.001;
  stf.objects.solids[0].material.fiber[0] = 0.0;
  stf.objects.solids[0].material.fiber[1] = 1.0;
  
  // create rigid body
  stf.objects.rigid_bodies.clear();
  stf.objects.rigid_bodies.resize(1);

  // create rigid nodes
  Nx = 24*2;
  Ny = 1*2;
  dx = 40/2;
  dy = 60/2;
  float R = 6*dx;
  float x0 = 0.1*VIEW_WIDTH;
  float y0 = 0.2*VIEW_HEIGHT;
  stf.objects.rigid_bodies[0].nodes.resize(Nx*Ny);
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      float x = dx*j;
      float y = dy*i;
      if (x < R) {
	float xc = R;
	float yc = (Ny-1)*dy + R;
	float theta = asin(1.0-x/R);
	float d = R + (Ny-1)*dy - y;
	x = xc - d*sin(theta);
	y = yc - d*cos(theta);
      } else if (x > (dx*(Nx-1)-R)) {
	float xc = dx*(Nx-1)-R;
	float yc = (Ny-1)*dy + R;
	float theta = asin((x - (dx*(Nx-1)-R))/R);
	float d = R + (Ny-1)*dy - y;
	x = xc + d*sin(theta);
	y = yc - d*cos(theta);
      }
      stf.objects.rigid_bodies[0].nodes[Nx*i+j].position[0] = x0+x;
      stf.objects.rigid_bodies[0].nodes[Nx*i+j].position[1] = y0+y;
    }
  }

  // create rigid elements
  stf.objects.rigid_bodies[0].elements.resize((Nx-1)*(Ny-1));
  for (int i = 0; i < (Ny-1); i++) {
    for (int j = 0; j < (Nx-1); j++) {
      stf.objects.rigid_bodies[0].elements[(Nx-1)*i+j] = new Quadrilateral();
      stf.objects.rigid_bodies[0].elements[(Nx-1)*i+j]->nodeIDs = { Nx*i+j, Nx*i+j+1, Nx*(i+1)+j+1, Nx*(i+1)+j };
    }
  }

  // create rigid material
  stf.objects.rigid_bodies[0].material.rho = 0.01;
  stf.objects.rigid_bodies[0].material.alpha = 0.01;

  // create contact surfaces
  stf.objects.contact_surfaces.clear();
  stf.objects.contact_surfaces.resize(1);
  stf.objects.contact_surfaces[0].mu = 0.5;
  stf.objects.contact_surfaces[0].penalty_stiffness = 500.0;
  stf.objects.contact_surfaces[0].damping_viscosity = 50.0;
  stf.objects.contact_surfaces[0].conductivity = 0.01;
  stf.objects.contact_surfaces[0].surface.clear();
  for (int j = 0; j < (Nx-1); j++) {
    int i = 0;
    stf.objects.contact_surfaces[0].surface.push_back(&stf.objects.rigid_bodies[0].nodes[Nx*i+j]);
  }
  for (int i = 0; i < (Ny-1); i++) {
    int j = Nx-1;
    stf.objects.contact_surfaces[0].surface.push_back(&stf.objects.rigid_bodies[0].nodes[Nx*i+j]);
  }
  for (int j = (Nx-1); j > 0; j--) {
    int i = Ny-1;
    stf.objects.contact_surfaces[0].surface.push_back(&stf.objects.rigid_bodies[0].nodes[Nx*i+j]);
  }
  for (int i = (Ny-1); i > 0; i--) {
    int j = 0;
    stf.objects.contact_surfaces[0].surface.push_back(&stf.objects.rigid_bodies[0].nodes[Nx*i+j]);
  }

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

  // create radiation source
  stf.environment.radiationSources.resize(9);
  for (int i = 0; i < 9; i++) {
    stf.environment.radiationSources[i] = new Radiation();
    stf.environment.radiationSources[i]->position[0] = (0.4 + 0.025*i)*VIEW_WIDTH;
    stf.environment.radiationSources[i]->position[1] = 0.12*VIEW_HEIGHT;
  }

  // set nodal constraint stiffness and viscosity
  stf.environment.nodalConstraint.stiffness = 50.0;
  stf.environment.nodalConstraint.viscosity = 50.0;
  
  // create boundary
  stf.environment.boundaries.resize(5);
  stf.environment.boundaries[0] = new Boundary();
  stf.environment.boundaries[0]->basePoint[0] = 0;
  stf.environment.boundaries[0]->basePoint[1] = 0.1*VIEW_HEIGHT;
  stf.environment.boundaries[0]->normal[0] = 0;
  stf.environment.boundaries[0]->normal[1] = 1;
  stf.environment.boundaries[0]->width = VIEW_WIDTH;
  stf.environment.boundaries[0]->temperature = 0.0;
  stf.environment.boundaries[0]->penalty_stiffness = 100.0;
  stf.environment.boundaries[0]->damping_viscosity = 10.0;
  stf.environment.boundaries[0]->conductivity = 0.01;
  stf.environment.boundaries[1] = new Boundary();
  stf.environment.boundaries[1]->basePoint[0] = 0;
  stf.environment.boundaries[1]->basePoint[1] = VIEW_HEIGHT;
  stf.environment.boundaries[1]->normal[0] = 1;
  stf.environment.boundaries[1]->normal[1] = 0;
  stf.environment.boundaries[1]->width = VIEW_HEIGHT;
  stf.environment.boundaries[1]->temperature = 0.0;
  stf.environment.boundaries[1]->penalty_stiffness = 100.0;
  stf.environment.boundaries[1]->damping_viscosity = 10.0;
  stf.environment.boundaries[1]->conductivity = 0.0;
  stf.environment.boundaries[2] = new Boundary();
  stf.environment.boundaries[2]->basePoint[0] = VIEW_WIDTH;
  stf.environment.boundaries[2]->basePoint[1] = 0;
  stf.environment.boundaries[2]->normal[0] =-1;
  stf.environment.boundaries[2]->normal[1] = 0;
  stf.environment.boundaries[2]->width = VIEW_WIDTH;
  stf.environment.boundaries[2]->temperature = 0.0;
  stf.environment.boundaries[2]->penalty_stiffness = 100.0;
  stf.environment.boundaries[2]->damping_viscosity = 10.0;
  stf.environment.boundaries[2]->conductivity = 0.0;
  stf.environment.boundaries[3] = new Boundary();
  stf.environment.boundaries[3]->basePoint[0] = VIEW_WIDTH;
  stf.environment.boundaries[3]->basePoint[1] = VIEW_HEIGHT;
  stf.environment.boundaries[3]->normal[0] = 0;
  stf.environment.boundaries[3]->normal[1] =-1;
  stf.environment.boundaries[3]->width = VIEW_HEIGHT;
  stf.environment.boundaries[3]->temperature = 0.0;
  stf.environment.boundaries[3]->penalty_stiffness = 100.0;
  stf.environment.boundaries[3]->damping_viscosity = 10.0;
  stf.environment.boundaries[3]->conductivity = 0.0;
  stf.environment.boundaries[4] = new Boundary();
  stf.environment.boundaries[4]->basePoint[0] = 0.3*VIEW_WIDTH;
  stf.environment.boundaries[4]->basePoint[1] = 0.2*VIEW_HEIGHT;
  stf.environment.boundaries[4]->normal[0] = 0;
  stf.environment.boundaries[4]->normal[1] = 1;
  stf.environment.boundaries[4]->width = 0.4*VIEW_WIDTH;
  stf.environment.boundaries[4]->temperature = 0.0;
  stf.environment.boundaries[4]->penalty_stiffness = 100.0;
  stf.environment.boundaries[4]->damping_viscosity = 10.0;
  stf.environment.boundaries[4]->conductivity = 0.0;

  // create fluid
  //stf.objects.fluid.Nx = 8*2;
  //stf.objects.fluid.Ny = 6*2;
  //stf.objects.fluid.dx = VIEW_WIDTH/stf.objects.fluid.Nx;
  //stf.objects.fluid.dy = VIEW_HEIGHT/stf.objects.fluid.Ny;
  //stf.objects.fluid.density0 = 1.0;
  //stf.objects.fluid.bulk_modulus = 10.0;
  //stf.objects.fluid.viscosity = 1.0;

  // initialize
  stf.initialize();

}

void Update(void) {
  int Nsubincrements = 50;
  float ddt = DT / Nsubincrements;
  // use a fixed number of time steps to prevent stalling
  for (int i = 0; i < Nsubincrements; i++) {
    stf.timeIntegrate(ddt);
  }
  vis_time += DT;

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

  // draw fluid
  //glBegin(GL_QUADS);
  //auto& fluid = stf.objects.fluid;
  //for (int j = 0; j < fluid.Ny; j++) {
  //  for (int i = 0; i < fluid.Nx; i++) {
  //    int id = fluid.Nx*j+i;
  //    float c = 1.0-1.0*(fluid.density[id]/fluid.density0);
  //    glColor4f(c,c,1.0,0.5);
  //    glVertex2f(fluid.dx*i,fluid.dy*j);
  //    glVertex2f(fluid.dx*(i+1),fluid.dy*j);
  //    glVertex2f(fluid.dx*(i+1),fluid.dy*(j+1));
  //    glVertex2f(fluid.dx*i,fluid.dy*(j+1));
  //  }
  //}
  //glEnd();

  // draw stove and heating elements
  glBegin(GL_QUADS);
  {
    float rgb[] = { 0.5f,0.5f,0.5f };
    float bbox[] = { 0.0f, 0.0f, VIEW_WIDTH, 0.1f*VIEW_HEIGHT };
    glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[1]);
    glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[1]);
    glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[3]);
    glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[3]);
  }
  {
    float rgb[] = { 0.3f,0.3f,0.3f };
    {
      float bbox[] = { 0.3f*VIEW_WIDTH, 0.18f*VIEW_HEIGHT, 0.7f*VIEW_WIDTH, 0.2f*VIEW_HEIGHT };
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[1]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[1]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[3]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[3]);
    }
    {
      float bbox[] = { 0.3f*VIEW_WIDTH, 0.1f*VIEW_HEIGHT, 0.32f*VIEW_WIDTH, 0.18f*VIEW_HEIGHT };
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[1]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[1]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[3]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[3]);
    }
    {
      float bbox[] = { 0.68f*VIEW_WIDTH, 0.1f*VIEW_HEIGHT, 0.7f*VIEW_WIDTH, 0.18f*VIEW_HEIGHT };
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[1]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[1]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[3]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[3]);
    }
    {
      float bbox[] = { 0.45f*VIEW_WIDTH, 0.1f*VIEW_HEIGHT, 0.55f*VIEW_WIDTH, 0.12f*VIEW_HEIGHT };
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[1]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[1]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[3]);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[3]);
    }
  }
  {
    float rgb[] = { 0.4f,0.4f,0.4f };
    float bbox[] = { 0.4f*VIEW_WIDTH, 0.12f*VIEW_HEIGHT, 0.6f*VIEW_WIDTH, 0.14f*VIEW_HEIGHT };
    glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[1]);
    glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[1]);
    glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[2],bbox[3]);
    glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(bbox[0],bbox[3]);
  }

  // draw solids
  //glBegin(GL_QUADS);
  for (auto &e : stf.objects.solids[0].elements) {
    if (e->active) {
      float nValues[4];
      if (colorPlot == 0) {
	float qValues[4];
	for (int j = 0; j < 4; j++) {
	  qValues[j] = e->integrationPoints[j].maxTemperature;
	}
	const float invW1 = 2521.0/1351.0;
	const float invW2 = -0.5;
	const float invW3 = 390.0/2911.0;
	nValues[0] = invW1*qValues[0] + invW2*qValues[1] + invW3*qValues[2] + invW2*qValues[3];
	nValues[1] = invW2*qValues[0] + invW1*qValues[1] + invW2*qValues[2] + invW3*qValues[3];
	nValues[2] = invW3*qValues[0] + invW2*qValues[1] + invW1*qValues[2] + invW2*qValues[3];
	nValues[3] = invW2*qValues[0] + invW3*qValues[1] + invW2*qValues[2] + invW1*qValues[3];
      } else {
	for (int j = 0; j < 4; j++) {
	  nValues[j] = stf.objects.solids[0].nodes[e->nodeIDs[j]].temperature;
	}
      }
      for (int j = 0; j < 4; j++) {
	int i = e->nodeIDs[j];
	//float fresh = e->integrationPoints[j].maxTemperature;
	//float charred = 1.0 - fresh;
	//glColor4f(0.7f*charred+0.2*fresh, 0.1f*charred+0.2*fresh, 0.1f*charred+0.2*fresh, 1);
	float rgb[3];
	if (colorPlot == 0) {
	  rawColor.getColor(nValues[j],rgb);
	} else {
	  thermalColor.getColor(nValues[j],rgb);
	}
	glColor4f(rgb[0],rgb[1],rgb[2],1);
	glVertex2f(stf.objects.solids[0].nodes[i].position[0],
		   stf.objects.solids[0].nodes[i].position[1]);
      }
    }
  }
  for (auto &e : stf.objects.rigid_bodies[0].elements) {
    if (e->active) {
      for (int j = 0; j < 4; j++) {
	int i = e->nodeIDs[j];
	float rgb[3];
	if (colorPlot == 0) {
	  rgb[0] = 0.1;
	  rgb[1] = 0.1;
	  rgb[2] = 0.1;
	} else {
	  thermalColor.getColor(stf.objects.rigid_bodies[0].nodes[i].temperature,rgb);
	}
	glColor4f(rgb[0],rgb[1],rgb[2],1);
	glVertex2f(stf.objects.rigid_bodies[0].nodes[i].position[0],
		   stf.objects.rigid_bodies[0].nodes[i].position[1]);
      }
    }
  }
  for (auto &r : stf.environment.radiationSources) {
    float x = r->position[0];
    float y = r->position[1];
    for (int i = 0; i < 3; i++) {
      float factor = 10*(3.0-i);
      float size = factor*r->heat/max_heat;
      float rgb[3];
      if (colorPlot == 0) {
  	rgb[0] = 1.0;
  	rgb[1] = 1.0-0.02*size;
  	rgb[2] = 0.0;
      } else {
  	thermalColor.getColor(0.1*size,rgb);
      }
      float freq = 1.0-0.2*i;
      float theta = sin(freq*vis_time/DT);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(x+0.1*0.3*theta*size,y-0.3*size);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(x+0.5*size,y+0.1*0.5*theta*size);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(x-0.1*1.0*theta*size,y+1.0*size);
      glColor4f(rgb[0],rgb[1],rgb[2],1); glVertex2f(x-0.5*size,y-0.1*0.5*theta*size);
    }
  }
  glEnd();

  // not sure what this is? maybe water?
  //glBegin(GL_QUADS);
  //glColor4f(0.5f, 0.4f, 0.0f, 0.5f);
  //glVertex2f(0,0);
  //glVertex2f(VIEW_WIDTH,0);
  //glVertex2f(VIEW_WIDTH,0.0*VIEW_HEIGHT);
  //glVertex2f(0,0.0*VIEW_HEIGHT);
  //glEnd();

  /*
  glColor4f(0.25f, 0.1f, 0.1f, 1);
  glBegin(GL_POINTS);
  for (auto &p : stf.objects.solids[0].nodes) {
    glVertex2f(p.position[0], p.position[1]);
  }
  glEnd();
  */

  // Draw cursor position
  if (stf.environment.nodalConstraint.active) {
    glColor4f(0.0f, 0.3f, 0.9f, 0.3f);
    glBegin(GL_POINTS);
    glVertex2f(stf.environment.nodalConstraint.cursor[0],stf.environment.nodalConstraint.cursor[1]);
    glEnd();
  }

  glutSwapBuffers();
}

void CheckDoneness(void) {
  float total    = 0.0;
  float raw      = 0.0; // [0.0, 0.2]
  float rare     = 0.0; // [0.2, 0.3]
  float medium   = 0.0; // [0.3, 0.5]
  float well     = 0.0; // [0.5, 0.6]
  float browned  = 0.0; // [0.6, 0.9]
  float charred  = 0.0; // [0.9, inf]
  for (auto &e : stf.objects.solids[0].elements) {
    if (e->active) {
      for (int j = 0; j < 4; j++) {
	float maxT = e->integrationPoints[j].maxTemperature;
	total += 1.0;
	if        (maxT < 0.2) {
	  raw     += 1.0;
	} else if (maxT < 0.3) {
	  rare    += 1.0;
	} else if (maxT < 0.5) {
	  medium  += 1.0;
	} else if (maxT < 0.6) {
	  well    += 1.0;
	} else if (maxT < 0.9) {
	  browned += 1.0;
	} else if (maxT >= 0.9) {
	  charred += 1.0;
	}
      }
    }
  }
  raw     /= total;
  rare    /= total;
  medium  /= total;
  well    /= total;
  browned /= total;
  charred /= total;
  std::cout << "Doneness: "  << std::endl;
  std::cout << "  Raw:     " << int(100*raw)     << "%" << std::endl;
  std::cout << "  Rare:    " << int(100*rare)    << "%" << std::endl;
  std::cout << "  Medium:  " << int(100*medium)  << "%" << std::endl;
  std::cout << "  Well:    " << int(100*well)    << "%" << std::endl;
  std::cout << "  Browned: " << int(100*browned) << "%" << std::endl;
  std::cout << "  Charred: " << int(100*charred) << "%" << std::endl;
  std::cout << std::endl;
} // CheckDoneness

void MoveMouse(int x, int y) {
  stf.environment.nodalConstraint.cursor[0] = 1.5*x;
  stf.environment.nodalConstraint.cursor[1] = VIEW_HEIGHT-1.5*y;
} // MoveMouse

void Mouse(int button, int state, int x, int y) {
  if (button == 0) {
    if (state == GLUT_DOWN) {
      stf.environment.nodalConstraint.active = true;
      MoveMouse(x,y);
    } else if (state == GLUT_UP) {
      stf.environment.nodalConstraint.active = false;
    }
  }
  glutMotionFunc(MoveMouse);
} // Mouse()

void Keyboard(unsigned char c, __attribute__((unused)) int x, __attribute__((unused)) int y) {
  switch(c) {
  case 'r': 
  case 'R':  
    InitSTF(); 
    break;
  case 'c':
  case 'C':
    colorPlot = (colorPlot+1)%2;
    if (colorPlot == 0) {
      glClearColor(0.9f,0.9f,0.9f,1);
    } else {
      glClearColor(0.0f,0.0f,0.0f,1);
    }
    break;
  case 'd':
  case 'D':
    CheckDoneness();
    break;
  case '+':
    for (auto r : stf.environment.radiationSources) r->heat = std::min(max_heat,r->heat+dheat);
    break;
  case '-':
    for (auto r : stf.environment.radiationSources) r->heat = std::max(min_heat,r->heat-dheat);
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
  glutMouseFunc(Mouse);
  glutKeyboardFunc(Keyboard);
  glutSpecialFunc(Arrows);

  InitGL();
  InitSTF();

  glutMainLoop();
  return 0;
}
