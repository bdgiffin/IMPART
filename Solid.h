#ifndef SOLID_H
#define SOLID_H

#include<vector>
#include<limits>
#include "BodyForce.h"
#include "Damping.h"
#include "Boundary.h"
#include "MaterialModel.h"
#include "Node.h"
#include "Element.h"

class Solid {
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
    
  } // initialize()

  float computeInternalForces(float dt) {

    // zero out all nodal forces
    for (uint a = 0; a < nodes.size(); a++) {
      nodes[a].zeroForces();
    } // for a = ...

    // set the shortest dimension
    float Lcrit = std::numeric_limits<float>::max();
    
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
	// check for element inversion, and optionally delete inverted element
	if (elements[e]->integrationPoints[i].relativeVolume < 1.0e-4) {
	  elements[e]->active = false;
	} // if (J < tol)
      } // for i = ...

      // scatter local force contributions from the current element to the nodes
      if (elements[e]->active) {
	for (uint a = 0; a < elements[e]->nodeIDs.size(); a++) {
	  nodes[elements[e]->nodeIDs[a]].sumForces(elementNodes[a]);
	} // for a = ...
	float Le = elements[e]->shortest_dimension(elementNodes);
	if (Le < Lcrit) Lcrit = Le;
      }

    } // for e = ...

    // return the critical stable time step
    return 1.0 * Lcrit / material.c;

  } // computeInternalForces()

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

    // use a semi-implicit time-integrator to update the nodal positions
    for (uint a = 0; a < nodes.size(); a++) {
      nodes[a].timeIntegrate(dt);
    } // for a = ...

  } // timeIntegrate()

  void applyBoundaryCondition(Boundary* boundary, float dt) {
    
    // update the nodal positions to satisfy the imposed boundary condition
    for (uint a = 0; a < nodes.size(); a++) {
      boundary->enforce(nodes[a], dt);
    } // for a = ...

  } // applyBoundaryCondition()

  MaterialModel material;
  std::vector<Node> nodes;
  std::vector<Element*> elements;
}; // Solid

#endif // SOLID_H
