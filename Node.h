#ifndef NODE_H
#define NODE_H

#include<vector>

class Node {
public:

  void initialize() {
    mass = 0.0;
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    force[0] = 0.0;
    force[1] = 0.0;
    temperature = 0.0;
  } // initialize()

  // sum mass contribution
  void sumMass(Node node) {
    mass += node.mass;
  } // sumMass()

  // zero out the current nodal forces
  void zeroForces() {
    force[0] = 0.0;
    force[1] = 0.0;
    heating  = 0.0;
  } // zeroForces()

  // sum force contributions
  void sumForces(Node node) {
    force[0] += node.force[0];
    force[1] += node.force[1];
    heating  += node.heating;
  } // zeroForces()

  // integrate velocity and position in time
  void timeIntegrate(float dt) {
    float dt_imass = dt / mass;
    velocity[0] += dt_imass * force[0];
    velocity[1] += dt_imass * force[1];
    //velocity[0] *= 0.99999;
    //velocity[1] *= 0.99999;
    position[0] += dt * velocity[0];
    position[1] += dt * velocity[1];
    temperature += dt * heating;
  } // timeIntegrate()

  float mass;
  float position[2];
  float velocity[2];
  float force[2];
  float temperature;
  float heating;
}; // Node

#endif // NODE_H
