#ifndef ELEMENT_H
#define ELEMENT_H

#include<vector>
#include<math.h>
#include "Node.h"
#include "MaterialPoint.h"

class Element {
public:

  virtual void initialize(std::vector<Node>& elementNodes) = 0;

  virtual float shortest_dimension(std::vector<Node>& elementNodes) = 0;

  bool active = true;
  std::vector<int> nodeIDs;
  std::vector<MaterialPoint> integrationPoints;
}; // Element
  
class Quadrilateral : public Element {
public:

  virtual void initialize(std::vector<Node>& elementNodes) {
    // declare const data for all quadrilateral elements
    constexpr float gaussPoint = 1.0 / sqrt(3.0);
    const float XA[4] = { -1.0, +1.0, +1.0, -1.0 };
    const float YA[4] = { -1.0, -1.0, +1.0, +1.0 };
    const float XQ[4] = { -gaussPoint, +gaussPoint, +gaussPoint, -gaussPoint };
    const float YQ[4] = { -gaussPoint, -gaussPoint, +gaussPoint, +gaussPoint };
    
    // create and initialize all integration point data
    integrationPoints.resize(4);

    // loop over all Gauss points
    for (int q = 0; q < 4; q++) {

      // initialize integration point shape function data
      integrationPoints[q].shapeFunctionValues.resize(4);
      integrationPoints[q].shapeFunctionGradients.resize(4);

      // loop over all nodes
      float jacobian[2][2] = { { 0.0, 0.0 }, { 0.0, 0.0 } };
      float dNa_dxi[4][2];
      for (int a = 0; a < 4; a++) {

	// compute shape function values at the current Gauss point
	integrationPoints[q].shapeFunctionValues[a] = 0.25 * (1.0 + XA[a] * XQ[q]) * (1.0 + YA[a] * YQ[q]);

        // compute shape function gradients w.r.t. the element's parent coordinates at the current Gauss point
        dNa_dxi[a][0] = 0.25 * XA[a] * (1.0 + YA[a] * YQ[q]);
	dNa_dxi[a][1] = 0.25 * YA[a] * (1.0 + XA[a] * XQ[q]);

	// sum contributions to the element Jacobian at the current Gauss point 
	for (int i = 0; i < 2; i++) {
	  for (int j = 0; j < 2; j++) {
	    jacobian[i][j] += elementNodes[a].position[i] * dNa_dxi[a][j];
	  } // j = ...
	} // i = ...
	
      } // a = ...

      // compute the determinant of the Jacobian at the current Gauss point
      integrationPoints[q].integrationWeight = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
      float jinv = 1.0 / integrationPoints[q].integrationWeight;

      // store the shape function gradients w.r.t. the element's initial coordinates at the current Gauss point
      for (int a = 0; a < 4; a++) {
	integrationPoints[q].shapeFunctionGradients[a][0] = jinv * (dNa_dxi[a][0] * jacobian[1][1] - dNa_dxi[a][1] * jacobian[1][0]);
	integrationPoints[q].shapeFunctionGradients[a][1] = jinv * (dNa_dxi[a][1] * jacobian[0][0] - dNa_dxi[a][0] * jacobian[0][1]);
      } // a = ...

    } // q = ...

  } // initialize()

  virtual float shortest_dimension(std::vector<Node>& elementNodes) {
    // compute the volume of the element
    float diag13x = elementNodes[2].position[0] - elementNodes[0].position[0];
    float diag13y = elementNodes[2].position[1] - elementNodes[0].position[1];
    float diag24x = elementNodes[3].position[0] - elementNodes[1].position[0];
    float diag24y = elementNodes[3].position[1] - elementNodes[1].position[1];
    float volume = diag13x*diag24y - diag13y*diag24x;

    // compute the length of each diagonal
    float diag13r = sqrt(diag13x*diag13x + diag13y*diag13y);
    float diag24r = sqrt(diag24x*diag24x + diag24y*diag24y);

    // return the shortest dimension
    if (diag13r > diag24r) {
      return volume / diag13r;
    } else {
      return volume / diag24r;
    }
  } // shortest_dimension

}; // Quadrilateral

#endif // ELEMENT_H
