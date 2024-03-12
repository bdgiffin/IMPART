#ifndef FLUID_H
#define FLUID_H

#include<vector>

class Fluid {
public:

  void initialize() {

    // initialize the grid fields
    density = std::vector<float>(Nx*Ny,density0);
    vx = std::vector<float>(Nx*Ny,0.0);
    vy = std::vector<float>(Nx*Ny,0.0);
    fx = std::vector<float>(Nx*Ny,0.0);
    fy = std::vector<float>(Nx*Ny,0.0);

    // initialize edge fields
    x_mass_flux = std::vector<float>((Nx-1)*Ny,0.0);
    x_momx_flux = std::vector<float>((Nx-1)*Ny,0.0);
    x_momy_flux = std::vector<float>((Nx-1)*Ny,0.0);
    y_mass_flux = std::vector<float>((Ny-1)*Nx,0.0);
    y_momx_flux = std::vector<float>((Ny-1)*Nx,0.0);
    y_momy_flux = std::vector<float>((Ny-1)*Nx,0.0);
    
  } // initialize()

  void computeInternalForces(float dt) {

    // zero out all forces
    float gravity = -100.0;
    for (int i = 0; i < Nx*Ny; i++) {
      fx[i] = 0.0;
      fy[i] = gravity * density[i];
    } // for i = ...

    // loop over all x-direction edges
    for (int i = 0; i < Ny; i++) {
      for (int e = 0; e < (Nx-1); e++) {
	// get the indicies of the grid cells
	int id1 = Nx*i+e;
	int id2 = Nx*i+e+1;
	
	// compute the density jump (gradient) at the current edge
	float drho = (density[id2] - density[id1])/dx;

	// compute the pressure
	float pressure = -bulk_modulus*drho;
	
	// compute the velocity jump (gradient) at the current edge
	float dvx = (vx[id2] - vx[id1])/dx;
	float dvy = (vy[id2] - vy[id1])/dx;

	// compute the viscous stresses
	float sxx = viscosity*dvx;
	float sxy = viscosity*dvy;

	// distribute forces to each grid cell
	fx[id1] -= pressure*dy - sxx*dy;
	fy[id1] += sxy*dy;
	fx[id2] += pressure*dy - sxx*dy;
	fy[id2] -= sxy*dy;
      }
    }

    // loop over all y-direction edges
    for (int e = 0; e < (Ny-1); e++) {
      for (int i = 0; i < Nx; i++) {
	// get the indicies of the grid cells
	int id1 = Nx*e+i;
	int id2 = Nx*(e+1)+i;
	
	// interpolate the density at the current edge
	float rho = 0.5*(density[id1] + density[id2]);

	// compute the pressure
	float pressure = bulk_modulus*(rho-density0);
	
	// compute the velocity gradient at the current edge
	float dvx = (vx[id2] - vx[id1])/dy;
	float dvy = (vy[id2] - vy[id1])/dy;

	// compute the viscous stresses
	float sxy = viscosity*dvx;
	float syy = viscosity*dvy;

	// distribute forces to each grid cell
	fx[id1] += sxy*dx;
	fy[id1] -= pressure*dx - syy*dx;
	fx[id2] -= sxy*dx;
	fy[id2] += pressure*dx - syy*dx;
      }
    }
    
  } // computeInternalForces()

  void timeIntegrate(float dt) {

    // loop over all grid cells and update velocity
    for (int i = 0; i < Nx*Ny; i++) {
      float mass = density[i]*dx*dy;
      float ax = fx[i]/mass;
      float ay = fy[i]/mass;
      vx[i] += ax*dt;
      vy[i] += ay*dt;
    } // for i = ...

    // loop over all x-direction edges and compute fluxes
    //for (int i = 0; i < Ny; i++) {
    //  for (int e = 0; e < (Nx-1); e++) {
    //	// get the indicies of the grid cells
    //	int id1 = Nx*i+e;
    //	int id2 = Nx*i+e+1;
    //	
    //	// interpolate the normal velocity at the current edge
    //	float vn = 0.5*(vx[id1] + vx[id2]);
    //
    //	// compute upwind edge flux
    //	int ie = (Nx-1)*i+e;
    //	if (vn > 0.0) {
    //	  x_mass_flux[ie] = density[id1]*dy*vn*dt;
    //	  x_momx_flux[ie] = density[id1]*vx[id1]*dy*vn*dt;
    //	  x_momy_flux[ie] = density[id1]*vy[id1]*dy*vn*dt;
    //	} else {
    //	  x_mass_flux[ie] = density[id2]*dy*vn*dt;
    //	  x_momx_flux[ie] = density[id2]*vx[id2]*dy*vn*dt;
    //	  x_momy_flux[ie] = density[id2]*vy[id2]*dy*vn*dt;
    //	}
    //  }
    //}

    // loop over all y-direction edges
    for (int e = 0; e < (Ny-1); e++) {
      for (int i = 0; i < Nx; i++) {
	// get the indicies of the grid cells
	int id1 = Nx*e+i;
	int id2 = Nx*(e+1)+i;
	
	// interpolate the normal velocity at the current edge
	float vn = 0.5*(vy[id1] + vy[id2]);

	// compute upwind edge flux
	int ie = Nx*e+i;
	if (vn > 0.0) {
	  y_mass_flux[ie] = density[id1]*dx*vn*dt;
	  y_momx_flux[ie] = density[id1]*vx[id1]*dx*vn*dt;
	  y_momy_flux[ie] = density[id1]*vy[id1]*dx*vn*dt;
	} else {
	  y_mass_flux[ie] = density[id2]*dx*vn*dt;
	  y_momx_flux[ie] = density[id2]*vx[id2]*dx*vn*dt;
	  y_momy_flux[ie] = density[id2]*vy[id2]*dx*vn*dt;
	}
      }
    }

    // loop over all grid cells and perform advection step
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
	// get current grid cell index
	int id = Nx*j+i;

	// get current grid cell's (extensive) conserved quantities
	float mass = density[id]*dx*dy;
	float momx = density[id]*vx[id]*dx*dy;
	float momy = density[id]*vy[id]*dx*dy;

	// update grid cell's conserved quantities via edge fluxes
	//if (i > 0) { // left edge flux
	//  int ie = (Nx-1)*j+(i-1);
	//  mass += x_mass_flux[ie];
	//  momx += x_momx_flux[ie];
	//  momy += x_momy_flux[ie];
	//}
	//if (i < (Nx-1)) { // right edge flux
	//  int ie = (Nx-1)*j+i;
	//  mass -= x_mass_flux[ie];
	//  momx -= x_momx_flux[ie];
	//  momy -= x_momy_flux[ie];
	//}
	if (j > 0) { // bottom edge flux
	  int ie = Nx*(j-1)+i;
	  mass += y_mass_flux[ie];
	  momx += y_momx_flux[ie];
	  momy += y_momy_flux[ie];
	}
	if (j < (Nx-1)) { // top edge flux
	  int ie = Nx*j+i;
	  mass -= y_mass_flux[ie];
	  momx -= y_momx_flux[ie];
	  momy -= y_momy_flux[ie];
	}

	// compute updated (intensive) grid cell quantities
	density[id] = mass/(dx*dy);
	if (density[id] < 0.1*density0) density[id] = 0.1*density0;
	mass        = density[id]*dx*dy;
	vx[id]      = momx/mass;
	vy[id]      = momy/mass;
      }
    }

  } // timeIntegrate()

  float Nx;           // number of grid cells in the x-direction
  float Ny;           // number of grid cells in the y-direction
  float dx;           // x-dimension of each grid cell
  float dy;           // y-dimension of each grid cell
  float density0;     // initial density of the fluid
  float bulk_modulus; // bulk modulus of the fluid
  float viscosity;    // kinematic viscosity of the fluid
  
  std::vector<float> density; // density in each grid cell
  std::vector<float> vx;      // x-velocity in each grid cell
  std::vector<float> vy;      // y-velocity in each grid cell
  std::vector<float> fx;      // x-force in each grid cell
  std::vector<float> fy;      // y-force in each grid cell
  
  std::vector<float> x_mass_flux; // x-edge fluxes
  std::vector<float> x_momx_flux; // x-edge fluxes
  std::vector<float> x_momy_flux; // x-edge fluxes
  std::vector<float> y_mass_flux; // y-edge fluxes
  std::vector<float> y_momx_flux; // y-edge fluxes
  std::vector<float> y_momy_flux; // y-edge fluxes
}; // Fluid

#endif // FLUID_H
