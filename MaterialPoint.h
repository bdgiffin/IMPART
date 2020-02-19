#ifndef MATERIAL_POINT_H
#define MATERIAL_POINT_H

#include<vector>
#include<array>
#include "MaterialPoint.h"

class MaterialPoint {
public:
  float stress[3];
  float relativeVolume;
  float flux[2];
  float mass;
  float integrationWeight;
  std::vector<float>    shapeFunctionValues;
  std::vector<std::array<float,2> > shapeFunctionGradients;
}; // MaterialPoint

#endif // MATERIAL_POINT_H
