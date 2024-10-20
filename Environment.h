#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include<vector>
#include "BodyForce.h"
#include "Damping.h"
#include "Boundary.h"
#include "Radiation.h"

class Environment {
public:
  
  void initialize() {

    // initialize all body forces
    for (int b = 0; b < bodyForces.size(); b++) {
      bodyForces[b]->initialize();
    } // for b = ...

    // initialize all damping forces
    for (int d = 0; d < dampingForces.size(); d++) {
      dampingForces[d]->initialize();
    } // for d = ...

    // initialize all boundaries
    for (int b = 0; b < boundaries.size(); b++) {
      boundaries[b]->initialize();
    } // for b = ...

    // initialize all radiation sources
    for (int r = 0; r < radiationSources.size(); r++) {
      radiationSources[r]->initialize();
    } // for r = ...

    // initialize nodal constraint
    nodalConstraint.initialize();

  } // initialize()

  void timeIntegrate(float dt) {

    // timeIntegrate all bodyForces
    for (int b = 0; b < bodyForces.size(); b++) {
      bodyForces[b]->timeIntegrate(dt);
    } // for b = ...

    // timeIntegrate all dampingForces
    for (int d = 0; d < dampingForces.size(); d++) {
      dampingForces[d]->timeIntegrate(dt);
    } // for d = ...

    // timeIntegrate all boundaries
    for (int b = 0; b < boundaries.size(); b++) {
      boundaries[b]->timeIntegrate(dt);
    } // for b = ...

    // timeIntegrate all radiation sources
    for (int r = 0; r < radiationSources.size(); r++) {
      radiationSources[r]->timeIntegrate(dt);
    } // for r = ...

  } // timeIntegrate()

  NodalConstraint         nodalConstraint;
  std::vector<BodyForce*> bodyForces;
  std::vector<Damping*>   dampingForces;
  std::vector<Boundary*>  boundaries;
  std::vector<Radiation*> radiationSources;
}; // Environment

#endif // ENVIRONMENT_H
