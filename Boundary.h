#ifndef BOUNDARY_H
#define BOUNDARY_H

#include<math.h>
#include "Node.h"

class Boundary {
public:

  void initialize() {
    mu = 0.5;
    temperature = 1.0;
  } // initialize()
  
  void timeIntegrate(float dt) {
    
  } // timeIntegrate()

  void enforce(Node& node, float dt) {
    float gap = (node.position[0] - basePoint[0]) * normal[0]
              + (node.position[1] - basePoint[1]) * normal[1];
    if (gap < 0.0) {
      // enforce contact
      node.position[0] -= gap * normal[0];
      node.position[1] -= gap * normal[1];
      float vgap = node.velocity[0] * normal[0] 
                 + node.velocity[1] * normal[1];
      node.velocity[0] -= 1.1 * vgap * normal[0];
      node.velocity[1] -= 1.1 * vgap * normal[1];
      float N = - mu * gap;
      float hgap = node.velocity[0] * normal[1] 
                 - node.velocity[1] * normal[0];
      if (abs(hgap) * dt < N) {
	// stick
	node.velocity[0] -= hgap * normal[1];
	node.velocity[1] += hgap * normal[0];
      } else {
        // slip
	node.velocity[0] *= 0.99999;
	node.velocity[1] *= 0.99999;
      }
      
      // conduct heat
      float tgap = node.temperature - temperature;
      node.heating += - conductivity * tgap;
    } // if (gap < 0.0)
  } // enforce

  void enforce_penalty(Node& node, float dt) {
    float gap = (node.position[0] - basePoint[0]) * normal[0]
              + (node.position[1] - basePoint[1]) * normal[1];
    if (gap < 0.0) {
      float vgap = node.velocity[0] * normal[0] 
                 + node.velocity[1] * normal[1];
      // enforce contact by applying a penalty force
      float penalty_force = - penalty_stiffness * gap - damping_viscosity * vgap;
      node.force[0] += penalty_force * normal[0];
      node.force[1] += penalty_force * normal[1];
      float hgap = node.velocity[0] * normal[1] 
                 - node.velocity[1] * normal[0];
      float slip_force = penalty_stiffness * hgap;
      if (abs(slip_force) > mu * penalty_force) {
	slip_force = (mu * penalty_force) * slip_force / abs(slip_force);
      }
      node.force[0] -= slip_force * normal[1];
      node.force[1] += slip_force * normal[0];
      
      // conduct heat
      float tgap = node.temperature - temperature;
      node.heating += - conductivity * tgap;
    } // if (gap < 0.0)
  } // enforce_penalty

  float basePoint[2];
  float normal[2];
  float temperature;
  float mu; // coefficient of friction
  float penalty_stiffness; // contact penalty stiffness
  float damping_viscosity; // contact damping viscosity
  float conductivity; // thermal conductivity
}; // Boundary

#endif // BOUNDARY_H
