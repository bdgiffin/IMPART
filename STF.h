#ifndef STF_H
#define STF_H

#include "Objects.h"
#include "Environment.h"

class STF {
public:

  void initialize() {
    environment.initialize();
    objects.initialize();
  } // initialize()

  void timeIntegrate(float dt) {
    environment.timeIntegrate(dt);
    
    // compute internal forces in all objects
    objects.computeInternalForces(dt);

    // apply external body forces to all objects
    for (uint i = 0; i < environment.bodyForces.size(); i++) {
      objects.applyBodyForce(environment.bodyForces[i]);
    } // i = ...

    // update unconstrained positions of all objects
    objects.timeIntegrate(dt);

    // apply damping to all objects
    for (uint i = 0; i < environment.dampingForces.size(); i++) {
      objects.applyDamping(environment.dampingForces[i],dt);
    } // i = ...

    // apply boundary conditions and enforce contact
    // apply external body forces to all objects
    for (uint i = 0; i < environment.boundaries.size(); i++) {
      objects.applyBoundaryCondition(environment.boundaries[i],dt);
    } // i = ...
  } // timeIntegrate()
  
  Objects objects;
  Environment environment;
}; // STF

#endif // STF_H
