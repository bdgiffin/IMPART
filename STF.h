#ifndef STF_H
#define STF_H

#include "Objects.h"
#include "Environment.h"

class STF {
public:

  void initialize() {
    dtold = 0.0;
    environment.initialize();
    objects.initialize();
  } // initialize()

  void timeIntegrate(float dtmax) {
    
    // compute internal forces in all objects and compute stable time step
    float dtcrit = objects.computeInternalForces(dtold);

    // disallow growth of the time step by more than 1% each step
    if ((dtold > 0.0) && (dtmax > 1.01*dtold)) dtmax = 1.01*dtold;

    // selectively limit the time step
    float dt;
    if (dtcrit < dtmax) {
      dt = dtcrit;
    } else {
      dt = dtmax;
    }

    // update the environment
    environment.timeIntegrate(dt);

    // apply external body forces to all objects
    for (int i = 0; i < environment.bodyForces.size(); i++) {
      objects.applyBodyForce(environment.bodyForces[i]);
    } // i = ...

    // apply penalty-enforced boundary conditions
    for (int i = 0; i < environment.boundaries.size(); i++) {
      objects.applyPenaltyBoundaryCondition(environment.boundaries[i],dt);
    } // i = ...

    // apply contact forces between all objects
    objects.applyContactForces(dt);

    // apply nodal constraint
    objects.applyNodalConstraint(environment.nodalConstraint,dt);

    // update unconstrained positions of all objects
    objects.timeIntegrate(dt);

    // apply damping to all objects
    for (int i = 0; i < environment.dampingForces.size(); i++) {
      objects.applyDamping(environment.dampingForces[i],dt);
    } // i = ...

    // apply boundary conditions and enforce contact
    // apply external body forces to all objects
    for (int i = 0; i < environment.boundaries.size(); i++) {
      objects.applyBoundaryCondition(environment.boundaries[i],dt);
    } // i = ...

    // save the old time step size
    dtold = dt;

  } // timeIntegrate()

  float dtold;
  Objects objects;
  Environment environment;
}; // STF

#endif // STF_H
