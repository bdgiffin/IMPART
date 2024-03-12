#ifndef RIGID_BODY_H
#define RIGID_BODY_H

#include<vector>
#include<limits>
#include "BodyForce.h"
#include "Damping.h"
#include "Boundary.h"
#include "RigidMaterialModel.h"
#include "Node.h"
#include "Element.h"
#include "NodalConstraint.h"

class RigidBody {
public:

  void initialize() {

    // initialize material
    material.initialize();

    // initialize all nodes
    for (uint a = 0; a < nodes.size(); a++) {
      nodes[a].initialize();
    } // for a = ...
    
    // initialize all elements
    for (uint e = 0; e < elements.size(); e++) {

      // gather local node data for the current element
      std::vector<Node> elementNodes;
      for (uint a = 0; a < elements[e]->nodeIDs.size(); a++) {
	elementNodes.push_back(nodes[elements[e]->nodeIDs[a]]);
      } // for a = ...

      // initialize each element in its undeformed configuration
      elements[e]->initialize(elementNodes);

      // loop over all integration points, initialize the material state,
      // and simultaneously accumulate masses at each local element node
      for (uint i = 0; i < elements[e]->integrationPoints.size(); i++) {
	material.initializeState(elements[e]->integrationPoints[i],elementNodes);
      } // for i = ...

      // scatter local mass contributions from the current element to the nodes 
      for (uint a = 0; a < elements[e]->nodeIDs.size(); a++) {
	nodes[elements[e]->nodeIDs[a]].sumMass(elementNodes[a]);
      } // for a = ...

    } // for e = ...

    // compute the total mass and initial location of the center of mass
    mass        = 0.0;
    position[0] = 0.0;
    position[1] = 0.0;
    for (auto& node : nodes) {
      mass        += node.mass;
      position[0] += node.mass*node.position[0];
      position[1] += node.mass*node.position[1];
    }
    position[0] /= mass;
    position[1] /= mass;

    // initialize the total rotation angle
    theta = 0.0;

    // compute the mass-moment of inertia about the center of mass
    mass_moment = 0.0;
    for (auto& node : nodes) {
      float rx = node.position[0]-position[0];
      float ry = node.position[1]-position[1];
      mass_moment += node.mass*(rx*rx + ry*ry);
    }

    // initialize the rigid translational and angular velocity
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    omega       = 0.0;
    
  } // initialize()

  void computeInternalForces(float dt) {

    // zero out force and torque applied directly to the rigid body
    force[0] = 0.0;
    force[1] = 0.0;
    torque   = 0.0;

    // zero out all nodal forces
    for (uint a = 0; a < nodes.size(); a++) {
      nodes[a].zeroForces();
    } // for a = ...

    // No internal mechanical forces are genereated by rigid elements,
    // but heat conduction is still evolved on the finite element partition
    
    // loop over all elements
    for (uint e = 0; e < elements.size(); e++) {

      // gather local node data for the current element
      std::vector<Node> elementNodes;
      for (uint a = 0; a < elements[e]->nodeIDs.size(); a++) {
	elementNodes.push_back(nodes[elements[e]->nodeIDs[a]]);
	elementNodes[a].zeroForces();
      } // for a = ...

      // loop over all integration points, update the material state,
      // and simultaneously accumulate forces at each local element node
      for (uint i = 0; i < elements[e]->integrationPoints.size(); i++) {
	material.updateState(elements[e]->integrationPoints[i],elementNodes,dt);
      } // for i = ...

      // scatter local force contributions from the current element to the nodes
      if (elements[e]->active) {
	for (uint a = 0; a < elements[e]->nodeIDs.size(); a++) {
	  nodes[elements[e]->nodeIDs[a]].sumForces(elementNodes[a]);
	} // for a = ...
      }

    } // for e = ...

  } // computeInternalForces()

  int find_nearest(float (&x)[2]) {
    // proximity search to find the index of the node closest to (x,y)
    int nearest = 0;
    float min_distance = std::numeric_limits<float>::max();
    for (uint i = 0; i < nodes.size(); i++) {
      float dx = x[0] - nodes[i].position[0];
      float dy = x[1] - nodes[i].position[1];
      float distance = sqrt(dx*dx + dy*dy);
      if (distance < min_distance) {
	min_distance = distance;
	nearest = i;
      }
    }
    return nearest;
  } // find_nearest()

  void applyNodalConstraint(NodalConstraint& nodalConstraint, float dt) {

    // only apply nodal constraint if it is active
    if (nodalConstraint.active) {
      // check to see if a node has been selected
      if (nodalConstraint.node == NULL) {
	int nearest = find_nearest(nodalConstraint.cursor);
	float dx = nodalConstraint.cursor[0] - nodes[nearest].position[0];
	float dy = nodalConstraint.cursor[1] - nodes[nearest].position[1];
	float distance = sqrt(dx*dx + dy*dy);
	if (distance < 20.0) {
	  nodalConstraint.node = &nodes[nearest];
	  nodalConstraint.theta0 = theta;
	}
      }
      // only apply nodal constraint if a node was found
      if (nodalConstraint.node != NULL) {
	float dx = nodalConstraint.cursor[0] - nodalConstraint.node->position[0];
	float dy = nodalConstraint.cursor[1] - nodalConstraint.node->position[1];
	nodalConstraint.node->force[0] += nodalConstraint.stiffness * dx - nodalConstraint.viscosity * nodalConstraint.node->velocity[0];
	nodalConstraint.node->force[1] += nodalConstraint.stiffness * dy - nodalConstraint.viscosity * nodalConstraint.node->velocity[1];
	torque += 1000000.0 * (- nodalConstraint.stiffness * (theta - nodalConstraint.theta0) - nodalConstraint.viscosity * omega);
      }
    } else {
      nodalConstraint.node = NULL;
    }

  } // applyNodalConstraint()

  void applyDamping(Damping* damping, float dt) {
    
    // apply damping to all nodes
    for (uint a = 0; a < nodes.size(); a++) {
      damping->applyDamping(nodes[a], dt);
    } // for a = ...

  } // applyDamping()

  void applyBodyForce(BodyForce* bodyForce) {
    
    // loop over all elements
    for (uint e = 0; e < elements.size(); e++) {

      // gather local node data for the current element
      std::vector<Node> elementNodes;
      for (uint a = 0; a < elements[e]->nodeIDs.size(); a++) {
	elementNodes.push_back(nodes[elements[e]->nodeIDs[a]]);
	elementNodes[a].zeroForces();
      } // for a = ...

      // loop over all integration points, evaluate the local body force,
      // and accumulate forces at each local element node
      for (uint i = 0; i < elements[e]->integrationPoints.size(); i++) {
	bodyForce->applyForce(elements[e]->integrationPoints[i],elementNodes);
      } // for i = ...

      // scatter local force contributions from the current element to the nodes 
      for (uint a = 0; a < elements[e]->nodeIDs.size(); a++) {
	nodes[elements[e]->nodeIDs[a]].sumForces(elementNodes[a]);
      } // for a = ...

    } // for e = ...

  } // applyBodyForce()

  void timeIntegrate(float dt) {

    // sum contributions to the net external force and torque
    for (auto& node : nodes) {
      float rx = node.position[0]-position[0];
      float ry = node.position[1]-position[1];
      force[0] += node.force[0];
      force[1] += node.force[1];
      torque   += rx*node.force[1]-ry*node.force[0];
    }

    // update the rigid translational and angular velocity
    float dt_imass = dt / mass;
    velocity[0] += dt_imass * force[0];
    velocity[1] += dt_imass * force[1];
    omega       += dt * torque / mass_moment;

    // compute and apply the incremental rigid rotation and translation to each node
    // using the Rodriguez rotation formula
    float s = sin(omega*dt);
    float c = cos(omega*dt);
    float cm1 = c - 1.0;
    for (auto& node : nodes) {
      float rx = node.position[0]-position[0];
      float ry = node.position[1]-position[1];
      node.position[0] += dt*velocity[0] + cm1*rx - s*ry;
      node.position[1] += dt*velocity[1] + s*rx   + cm1*ry;
      // zero the individual nodal forces and velocities,
      // but still update the nodal temperatures individually
      node.velocity[0] = 0.0;
      node.velocity[1] = 0.0;
      node.force[0] = 0.0;
      node.force[1] = 0.0;
      node.timeIntegrate(dt);
      // set the nodal velocity consistent with the rigid body velocity
      node.velocity[0] = velocity[0] - velocity[2]*ry;
      node.velocity[1] = velocity[1] + velocity[2]*rx;
    }

    // update the center of mass location and the total rotation angle
    position[0] += dt*velocity[0];
    position[1] += dt*velocity[1];
    theta       += dt*omega;

  } // timeIntegrate()

  void applyPenaltyBoundaryCondition(Boundary* boundary, float dt) {
    
    // update the nodal positions to satisfy the imposed boundary condition
    for (uint a = 0; a < nodes.size(); a++) {
      boundary->enforce_penalty(nodes[a], dt);
    } // for a = ...

  } // applyPenaltyBoundaryCondition()

  float mass;        // total mass
  float mass_moment; // mass-moment of inertia about the center of mass
  float position[2]; // center of mass
  float theta;       // total rotation
  float velocity[2]; // velocity at the center of mass
  float omega;       // angular velocity
  float force[2];    // net force at the center of mass
  float torque;      // net externally applied torque
  RigidMaterialModel material;
  std::vector<Node> nodes;
  std::vector<Element*> elements;
}; // RigidBody

#endif // RIGID_BODY_H
