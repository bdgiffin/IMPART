#ifndef MATERIAL_MODEL_H
#define MATERIAL_MODEL_H

#include<vector>
#include<math.h>
#include "Node.h"
#include "MaterialPoint.h"

class MaterialModel {
public:

  void initialize() {
    
    // compute Lame parameters
    G = 0.5 * E / (1.0 + nu);
    kappa = E / ( 3.0 * (1.0 - 2.0 * nu) );
    Em = 0.0;
    Et = 0.0;

    // compute the sound speed
    c = sqrt(kappa/rho);

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

    // compute the mechanical deformation gradient Fe = dx/dX
    float Fe[2][2] = { { 0.0, 0.0 }, { 0.0, 0.0 } };
    for (uint a = 0; a < elementNodes.size(); a++) {
      for (int i = 0; i < 2; i++) {
	for (int j = 0; j < 2; j++) {
	  Fe[i][j] += elementNodes[a].position[i] * integrationPoint.shapeFunctionGradients[a][j];
	} // j = ...
      } // i = ...
    } // a = ...

    // interpolate the temperature
    float T = 0.0;
    for (uint a = 0; a < elementNodes.size(); a++) {
      T += elementNodes[a].temperature * integrationPoint.shapeFunctionValues[a];
    } // a = ...

    // ============================================================================ //

    // update the maximum temperature
    if (T > integrationPoint.maxTemperature) integrationPoint.maxTemperature = T;

    // define the fiber stretch tensor: Fm = Fm0 * fiber (x) fiber + Fm1 * (I - fiber) (x) (I - fiber)
    float Fm0 = exp(Em); // normal stretch
    float Fm1 = exp(Et); // transverse stretch
    float Fm[2][2] = { { Fm0 * fiber[0] * fiber[0] + Fm1 * fiber[1] * fiber[1], Fm0 * fiber[0] * fiber[1] - Fm1 * fiber[1] * fiber[0] }, 
                       { Fm0 * fiber[1] * fiber[0] - Fm1 * fiber[0] * fiber[1], Fm0 * fiber[1] * fiber[1] + Fm1 * fiber[0] * fiber[0] } };

    // apply muscle stretch: F = Fe * Fm
    float F[2][2] = { { Fe[0][0] * Fm[0][0] + Fe[0][1] * Fm[1][0], Fe[0][0] * Fm[0][1] + Fe[0][1] * Fm[1][1] }, 
                      { Fe[1][0] * Fm[0][0] + Fe[1][1] * Fm[1][0], Fe[1][0] * Fm[0][1] + Fe[1][1] * Fm[1][1] } };
    
    // decay muscle stretches
    Em *= 0.999999;
    Et *= 0.999999;

    // compute determinant, and its inverse
    float detF = (F[0][0] * F[1][1] - F[0][1] * F[1][0]);
    float idetF = 1.0 / detF;

    // compute volumetric strain rate
    float idt = 0.0;
    if (dt > 0.0) idt = 1.0 / dt;
    float dvol = log(detF / integrationPoint.relativeVolume) * idt;

    // compute bulk viscosity pressure
    float viscosity = 10.0;
    float q = - viscosity * dvol;

    // compute the B tensor
    float B11 = F[0][0] * F[0][0] + F[0][1] * F[0][1]; // F_1i * F_1i = F_11 * F_11 + F_12 * F_12
    float B22 = F[1][0] * F[1][0] + F[1][1] * F[1][1]; // F_2i * F_2i = F_21 * F_21 + F_22 * F_22
    float B12 = F[0][0] * F[1][0] + F[0][1] * F[1][1]; // F_1i * F_2i = F_11 * F_21 + F_12 * F_22

    // compute the trace of B, and subsequently dev(B)
    float trB = (B11 + B22 + 1.0) / 3.0;
    B11 -= trB;
    B22 -= trB;

    // compute the deviatoric stress
    float GdivJ = G * pow(idetF, 2.0 / 3.0);
    B11 *= GdivJ;
    B22 *= GdivJ;
    B12 *= GdivJ;
    
    // compute the pressure and the Kirchhoff stress, including thermal strain and bulk viscosity
    float p = - kappa * (detF * (detF - 1.0) - dVdT * integrationPoint.maxTemperature);
    B11 -= p + q;
    B22 -= p + q;

    // ============================================================================ //

    // compute and store the Cauchy stress
    integrationPoint.stress[0] = idetF * B11;
    integrationPoint.stress[1] = idetF * B22;
    integrationPoint.stress[2] = idetF * B12;

    // store the relative volume
    integrationPoint.relativeVolume = detF;

    // compute the 1st Piola-Kirchhoff stress (P = tau * F^{-T})
    float P[2][2] = { { (integrationPoint.stress[0] * F[0][0] - integrationPoint.stress[2] * F[0][1]), (-integrationPoint.stress[0] * F[1][0] + integrationPoint.stress[2] * F[1][1]) }, 
                      { (integrationPoint.stress[2] * F[0][0] - integrationPoint.stress[1] * F[0][1]), (-integrationPoint.stress[2] * F[1][0] + integrationPoint.stress[1] * F[1][1]) } };

    // compute the heat flux vector
    integrationPoint.flux[0] = 0.0;
    integrationPoint.flux[1] = 0.0;
    for (uint a = 0; a < elementNodes.size(); a++) {
      for (int i = 0; i < 2; i++) {
	integrationPoint.flux[i] -= alpha * elementNodes[a].temperature * integrationPoint.shapeFunctionGradients[a][i];
      } // i = ...
    } // a = ...

    // ============================================================================ //

    // sum the nodal force contributions
    for (uint a = 0; a < elementNodes.size(); a++) {
      for (int i = 0; i < 2; i++) {
	for (int j = 0; j < 2; j++) {
	  elementNodes[a].force[i] -= P[i][j] * integrationPoint.shapeFunctionGradients[a][j] * integrationPoint.integrationWeight;
	} // j = ...
      } // i = ...
    } // a = ...

    // sum the heat flux divergence contributions
    for (uint a = 0; a < elementNodes.size(); a++) {
      for (int i = 0; i < 2; i++) {
	elementNodes[a].heating += integrationPoint.flux[i] * integrationPoint.shapeFunctionGradients[a][i] * integrationPoint.integrationWeight;
      } // i = ...
    } // a = ...

    // ============================================================================ //

  } // updateState()

  float rho;      // material density
  float E;        // elastic (Young's) modulus
  float nu;       // Poisson's ratio
  float G;        // shear modulus
  float kappa;    // bulk modulus
  float c;        // sound speed
  float dVdT;     // coefficient of thermal expansion
  float alpha;    // thermal diffusivity = k / (c * rho)
  float Em;       // axial muscle strain
  float Et;       // transverse muscle strain
  float fiber[2]; // fiber direction
}; // MaterialModel

#endif // MATERIAL_MODEL_H
