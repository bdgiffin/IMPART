#ifndef NODAL_CONSTRAINT_H
#define NODAL_CONSTRAINT_H

#include "Node.h"

class NodalConstraint {
public:

  void initialize() {
    active = false;
    node = NULL;
  } // initialize()

  bool  active;
  float cursor[2];
  Node* node;
  float theta0; // reference rotation for rigid body
  float stiffness;
  float viscosity;
}; // NodalConstraint

#endif // NODAL_CONSTRAINT_H
