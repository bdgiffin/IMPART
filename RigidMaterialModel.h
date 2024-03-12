#ifndef RIGID_MATERIAL_MODEL_H
#define RIGID_MATERIAL_MODEL_H

#include<vector>
#include<math.h>
#include "Node.h"
#include "MaterialPoint.h"

class RigidMaterialModel {
public:

  void initialize() {
    
    // do nothing
    
  } // initialize()

  void initializeState(MaterialPoint& integrationPoint, std::vector<Node>& elementNodes) {
    
    // initialize stress state
    integrationPoint.stress[0] = 0.0;
    integrationPoint.stress[1] = 0.0;
    integrationPoint.stress[2] = 0.0;

    // initialize relative volume
    integrationPoint.relativeVolume = 1.0;

    // accumulate nodal masses
    integrationPoint.mass = rho * integrationPoint.integrationWeight;
    for (uint a = 0; a < elementNodes.size(); a++) {
      elementNodes[a].mass += integrationPoint.mass * integrationPoint.shapeFunctionValues[a];
    } // for a = ...

    // initialize max temperature
    integrationPoint.maxTemperature = 0.0;

  } // initializeState()

  void updateState(MaterialPoint& integrationPoint, std::vector<Node>& elementNodes, float dt) {

    // interpolate the temperature
    float T = 0.0;
    for (uint a = 0; a < elementNodes.size(); a++) {
      T += elementNodes[a].temperature * integrationPoint.shapeFunctionValues[a];
    } // a = ...

    // ============================================================================ //

    // compute the heat flux vector
    integrationPoint.flux[0] = 0.0;
    integrationPoint.flux[1] = 0.0;
    for (uint a = 0; a < elementNodes.size(); a++) {
      for (int i = 0; i < 2; i++) {
	integrationPoint.flux[i] -= alpha * elementNodes[a].temperature * integrationPoint.shapeFunctionGradients[a][i];
      } // i = ...
    } // a = ...

    // ============================================================================ //

    // sum the heat flux divergence contributions
    for (uint a = 0; a < elementNodes.size(); a++) {
      for (int i = 0; i < 2; i++) {
	elementNodes[a].heating += integrationPoint.flux[i] * integrationPoint.shapeFunctionGradients[a][i] * integrationPoint.integrationWeight;
      } // i = ...
    } // a = ...

    // ============================================================================ //

  } // updateState()

  float rho;      // material density
  float alpha;    // thermal diffusivity = k / (c * rho)
}; // RigidMaterialModel

#endif // RIGID_MATERIAL_MODEL_H
