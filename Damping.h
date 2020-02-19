#ifndef DAMPING_H
#define DAMPING_H

#include<vector>
#include<math.h>
#include "Node.h"
#include "MaterialPoint.h"

class Damping {
public:

  void initialize() {

  } // initialize()

  void timeIntegrate(float dt) {

  } // timeIntegrate()

  virtual void applyDamping(Node& node, float dt) = 0;

}; // Damping

class MassDamping : public Damping {
public:

  void applyDamping(Node& node, float dt) {

    if (node.position[1] < y_surface) {
      // apply damping to all nodes below surface
      float drag = 0.5 * node.mass * viscosity * sqrt(node.velocity[0] * node.velocity[0] + node.velocity[1] * node.velocity[1]);
      float fx = - drag * node.velocity[0];
      float fy = - drag * node.velocity[1];
      node.velocity[0] += fx * dt;
      node.velocity[1] += fy * dt;
      node.position[0] += fx * dt * dt;
      node.position[1] += fy * dt * dt;
    }

  } // applyDamping()

  float y_surface; // y-coordinate of fluid surface
  float viscosity; // viscous damping
}; // MassDamping

#endif // BODY_FORCE_H
