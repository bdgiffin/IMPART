#ifndef OBJECTS_H
#define OBJECTS_H

#include<vector>
#include "BodyForce.h"
#include "Damping.h"
#include "Particle.h"
#include "Solid.h"

class Objects {
public:

  void initialize() {

    // initialize all particles
    for (uint p = 0; p < particles.size(); p++) {
      particles[p].initialize();
    } // for p = ...

    // initialize all solids
    for (uint s = 0; s < solids.size(); s++) {
      solids[s].initialize();
    } // for s = ...

  } // initialize()

  void timeIntegrate(float dt) {

    // timeIntegrate all particles
    for (uint p = 0; p < particles.size(); p++) {
      particles[p].timeIntegrate(dt);
    } // for p = ...

    // timeIntegrate all solids
    for (uint s = 0; s < solids.size(); s++) {
      solids[s].timeIntegrate(dt);
    } // for s = ...

  } // timeIntegrate()

  void computeInternalForces(float dt) {

    // compute internal forces all solids
    for (uint s = 0; s < solids.size(); s++) {
      solids[s].computeInternalForces(dt);
    } // for s = ...

  } // computeInternalForces()

  void applyBodyForce(BodyForce* bodyForce) {

    // apply body force to all solids
    for (uint s = 0; s < solids.size(); s++) {
      solids[s].applyBodyForce(bodyForce);
    } // for s = ...

  } // applyBodyForce()

  void applyDamping(Damping* damping, float dt) {

    // apply damping to all solids
    for (uint s = 0; s < solids.size(); s++) {
      solids[s].applyDamping(damping,dt);
    } // for s = ...

  } // applyDamping()

  void applyBoundaryCondition(Boundary* boundary, float dt) {

    // apply boundary condition to all solids
    for (uint s = 0; s < solids.size(); s++) {
      solids[s].applyBoundaryCondition(boundary,dt);
    } // for s = ...

  } // applyBoundaryCondition()

  // fluid particles
  float radius;
  float mass;
  std::vector<Node> particles;

  // solid objects
  std::vector<Solid> solids;
}; // Objects

#endif // OBJECTS_H
