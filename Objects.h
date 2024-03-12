#ifndef OBJECTS_H
#define OBJECTS_H

#include<vector>
#include<limits>
#include "BodyForce.h"
#include "Damping.h"
#include "Particle.h"
#include "Solid.h"
#include "RigidBody.h"
#include "Contact.h"
#include "Fluid.h"

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

    // initialize all rigid bodies
    for (uint r = 0; r < rigid_bodies.size(); r++) {
      rigid_bodies[r].initialize();
    } // for r = ...

    // initialize the fluid
    fluid.initialize();

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

    // timeIntegrate all rigid bodies
    for (uint r = 0; r < rigid_bodies.size(); r++) {
      rigid_bodies[r].timeIntegrate(dt);
    } // for r = ...

    // timeIntegrate the fluid
    fluid.timeIntegrate(dt);
    
  } // timeIntegrate()

  float computeInternalForces(float dt) {

    // set the critical stable time step
    float dtcrit = std::numeric_limits<float>::max();

    // compute internal (thermal-mechanical) forces all solids
    for (uint s = 0; s < solids.size(); s++) {
      float dts = solids[s].computeInternalForces(dt);
      if (dts < dtcrit) dtcrit = dts;
    } // for s = ...

    // compute internal (thermal) forces all rigid bodies
    for (uint r = 0; r < rigid_bodies.size(); r++) {
      rigid_bodies[r].computeInternalForces(dt);
    } // for r = ...

    // compute internal forces for the fluid
    fluid.computeInternalForces(dt);

    // return the limiting stable time step
    return dtcrit;

  } // computeInternalForces()

  void applyNodalConstraint(NodalConstraint& nodalConstraint, float dt) {

    // apply to all rigid bodies
    for (uint r = 0; r < rigid_bodies.size(); r++) {
      rigid_bodies[r].applyNodalConstraint(nodalConstraint,dt);
    } // for r = ...

  } // applyNodalConstraint()

  void applyBodyForce(BodyForce* bodyForce) {

    // apply body force to all solids
    for (uint s = 0; s < solids.size(); s++) {
      solids[s].applyBodyForce(bodyForce);
    } // for s = ...

    // apply body force to all rigid bodies
    for (uint r = 0; r < rigid_bodies.size(); r++) {
      rigid_bodies[r].applyBodyForce(bodyForce);
    } // for r = ...

  } // applyBodyForce()

  void applyDamping(Damping* damping, float dt) {

    // apply damping to all solids
    for (uint s = 0; s < solids.size(); s++) {
      solids[s].applyDamping(damping,dt);
    } // for s = ...

    // apply damping to all rigid bodies
    for (uint r = 0; r < rigid_bodies.size(); r++) {
      rigid_bodies[r].applyDamping(damping,dt);
    } // for r = ...

  } // applyDamping()

  void applyBoundaryCondition(Boundary* boundary, float dt) {

    // apply boundary condition to all solids
    for (uint s = 0; s < solids.size(); s++) {
      solids[s].applyBoundaryCondition(boundary,dt);
    } // for s = ...

  } // applyBoundaryCondition()

  void applyPenaltyBoundaryCondition(Boundary* boundary, float dt) {

    // apply boundary condition to all rigid bodies
    for (uint r = 0; r < rigid_bodies.size(); r++) {
      rigid_bodies[r].applyPenaltyBoundaryCondition(boundary,dt);
    } // for r = ...

  } // applyPenaltyBoundaryCondition()

  void applyContactForces(float dt) {

    // enforce solid node-to-surface contact
    for (uint c = 0; c < contact_surfaces.size(); c++) {
      for (uint s = 0; s < solids.size(); s++) {
	contact_surfaces[c].enforce_penalty(solids[s].nodes,dt);
      } // for s = ...
    } // for c = ...

  } // applyContact()

  // fluid particles
  float radius;
  float mass;
  std::vector<Node> particles;

  // solid objects
  std::vector<Solid> solids;

  // rigid bodies
  std::vector<RigidBody> rigid_bodies;

  // contact surfaces
  std::vector<Contact> contact_surfaces;

  // fluid
  Fluid fluid;
}; // Objects

#endif // OBJECTS_H
