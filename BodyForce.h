#ifndef BODY_FORCE_H
#define BODY_FORCE_H

#include<vector>
#include "Node.h"
#include "MaterialPoint.h"

class BodyForce {
public:

  void initialize() {

  } // initialize()

  void timeIntegrate(float dt) {

  } // timeIntegrate()

  virtual void applyForce(MaterialPoint& integrationPoint, std::vector<Node>& elementNodes) = 0;

}; // BodyForce

class Gravity : public BodyForce {
public:

  void applyForce(MaterialPoint& integrationPoint, std::vector<Node>& elementNodes) {

    // sum the nodal force contributions
    for (int a = 0; a < elementNodes.size(); a++) {
      elementNodes[a].force[1] -= gravity * integrationPoint.mass;
    } // a = ...

  } // applyForce()

  float gravity; // gravitational acceleration
}; // Gravity

class BuoyantForce : public BodyForce {
public:

  void applyForce(MaterialPoint& integrationPoint, std::vector<Node>& elementNodes) {
    
    // compute the interpolated current coordinate
    float x[2] = { 0.0, 0.0 };
    for (int a = 0; a < elementNodes.size(); a++) {
      for (int i = 0; i < 2; i++) {
	x[i] += elementNodes[a].position[i] * integrationPoint.shapeFunctionValues[a];
      } // i = ...
    } // a = ...

    // compute the bouyant force at the specified coordinate
    float b = 0.0;
    if (x[1] < y_surface) {
      b = gravity * rho_fluid;
    }

    // sum the nodal force contributions
    for (int a = 0; a < elementNodes.size(); a++) {
      elementNodes[a].force[1] += b * integrationPoint.integrationWeight;
    } // a = ...

  } // applyForce()

  float gravity;   // gravitational acceleration
  float rho_fluid; // fluid density
  float y_surface; // y-coordinate of fluid surface
}; // BuoyantForce

#endif // BODY_FORCE_H
