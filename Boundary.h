#ifndef BOUNDARY_H
#define BOUNDARY_H

#include<math.h>
#include "Node.h"

class Boundary {
public:

  void initialize() {
    mu = 0.5;
    temperature = 1.0;
  } // initialize()
  
  void timeIntegrate(float dt) {
    
  } // timeIntegrate()

  void enforce(Node& node, float dt) {
    float gap = (node.position[0] - basePoint[0]) * normal[0]
              + (node.position[1] - basePoint[1]) * normal[1];
    if (gap < 0.0) {
      // enforce contact
      node.position[0] -= gap * normal[0];
      node.position[1] -= gap * normal[1];
      float vgap = node.velocity[0] * normal[0] 
                 + node.velocity[1] * normal[1];
      node.velocity[0] -= 1.1 * vgap * normal[0];
      node.velocity[1] -= 1.1 * vgap * normal[1];
      float N = - mu * gap;
      float hgap = node.velocity[0] * normal[1] 
                 - node.velocity[1] * normal[0];
      if (abs(hgap) * dt < N) {
	// stick
	node.velocity[0] -= hgap * normal[1];
	node.velocity[1] += hgap * normal[0];
      } else {
        // slip
	node.velocity[0] *= 0.99999;
	node.velocity[1] *= 0.99999;
      }
      
      // conduct heat
      float tgap = temperature - node.temperature;
      node.temperature += tgap * 0.00001;
    } // if (gap < 0.0)
  } // enforce

  float basePoint[2];
  float normal[2];
  float temperature;
  float mu; // coefficient of friction
}; // Boundary

#endif // BOUNDARY_H
