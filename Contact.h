#ifndef CONTACT_H
#define CONTACT_H

#include<math.h>
#include<vector>
#include<limits>
#include "Node.h"

class Contact {
public:

  void initialize() {

  } // initialize()

  int find_nearest(Node& node) {
    // proximity search to find the index of the surface node closest to the input node
    int nearest = 0;
    float min_distance = std::numeric_limits<float>::max();
    for (int i = 0; i < surface.size(); i++) {
      float dx = node.position[0] - surface[i]->position[0];
      float dy = node.position[1] - surface[i]->position[1];
      float distance = sqrt(dx*dx + dy*dy);
      if (distance < min_distance) {
	min_distance = distance;
	nearest = i;
      }
    }
    return nearest;
  } // find_nearest()

  void enforce_penalty(std::vector<Node>& nodes, float dt) {
    for (auto& node : nodes) enforce_penalty(node,dt);
  } // enforce_penalty()

  void enforce_penalty(Node& node, float dt) {
    // find the nearest node on the surface, and get the two adjoining segments
    int nearest = find_nearest(node);
    Node* segments[3];
    segments[1] = surface[nearest];
    if (nearest == 0) {
      segments[0] = surface[surface.size()-1];
      segments[2] = surface[nearest+1];
    } else if (nearest == int(surface.size())-1) {
      segments[0] = surface[nearest-1];
      segments[2] = surface[0];
    } else {
      segments[0] = surface[nearest-1];
      segments[2] = surface[nearest+1];
    }

    // compute the length and normal of each segment
    float L[2];
    float Nx[2];
    float Ny[2];
    float xi[2];
    float gap[2];
    bool  on[2];
    for (int i = 0; i < 2; i++) {
      Ny[i] = -(segments[i+1]->position[0] - segments[i]->position[0]);
      Nx[i] = +(segments[i+1]->position[1] - segments[i]->position[1]);
      L[i]  = sqrt(Nx[i]*Nx[i] + Ny[i]*Ny[i]);
      Nx[i] /= L[i];
      Ny[i] /= L[i];
      // get the relative xi coordinate on the current segment
      float rax = segments[i]->position[0] - node.position[0];
      float ray = segments[i]->position[1] - node.position[1];
      xi[i] = (rax*Ny[i]-ray*Nx[i]) / L[i];
      on[i] = (xi[i] >= 0.0) && (xi[i] <= 1.0);
      // compute the gap (a positive gap implies penetration)
      gap[i] = rax*Nx[i]+ray*Ny[i];
    }

    // apply contact force and heating for a given node-to-segment interaction
    auto apply_contact_force = [&](int i, float x, float dx, float dy, float d, float scale) {
      // get the velocity relative to the projected point on the surface
      float vx = node.velocity[0] - ((1.0-x) * segments[i+0]->velocity[0] + x * segments[i+1]->velocity[0]);
      float vy = node.velocity[1] - ((1.0-x) * segments[i+0]->velocity[1] + x * segments[i+1]->velocity[1]);
      float vgap = vx*dx + vy*dy;
      float hgap = vx*dy - vy*dx;
      // the temperature gap relative to the projected point on the surface
      float tgap = node.temperature - ((1.0-x) * segments[i+0]->temperature + x * segments[i+1]->temperature);
      // enforce contact by applying a penalty force
      float penalty_force = scale * (penalty_stiffness * d - damping_viscosity * vgap);
      float slip_force = scale * penalty_stiffness * hgap;
      if (abs(slip_force) > mu * penalty_force) {
        slip_force = (mu * penalty_force) * slip_force / abs(slip_force);
      }
      float fx = penalty_force * dx - slip_force * dy;
      float fy = penalty_force * dy + slip_force * dx;
      float h  = - scale * conductivity * tgap;
      node.force[0]           += fx;
      node.force[1]           += fy;
      node.heating            += h;
      segments[i+0]->force[0] -= (1.0-x) * fx;
      segments[i+0]->force[1] -= (1.0-x) * fy;
      segments[i+0]->heating  -= (1.0-x) * h;
      segments[i+1]->force[0] -= x * fx;
      segments[i+1]->force[1] -= x * fy;
      segments[i+1]->heating  -= x * h;
    };

    // handle different cases based on the penetrating conditions
    if (on[0] && on[1]) {
      if ((gap[0] > 0.0) && (gap[1] > 0.0)) {
    	float segment_weight[2];
    	segment_weight[0] = gap[1]/(gap[0]+gap[1]);
    	segment_weight[1] = 1.0 - segment_weight[0];
    	apply_contact_force(0,xi[0],Nx[0],Ny[0],gap[0],segment_weight[0]);
    	apply_contact_force(1,xi[1],Nx[1],Ny[1],gap[1],segment_weight[1]);
      }
    } else if (on[0]) {
      if (gap[0] > 0.0) apply_contact_force(0,xi[0],Nx[0],Ny[0],gap[0],1.0);
    } else if (on[1]) {
      if (gap[1] > 0.0) apply_contact_force(1,xi[1],Nx[1],Ny[1],gap[1],1.0);
    } else {
      if ((gap[0] > 0.0) && (gap[1] > 0.0)) {
      	// enforce contact directly with the nearest node
      	float gx = segments[1]->position[0] - node.position[0];
      	float gy = segments[1]->position[1] - node.position[1];
      	float gn = sqrt(gx*gx + gy*gy);
      	gx /= gn;
      	gy /= gn;
	apply_contact_force(0,1.0,gx,gy,gn,1.0);
      }
    }
    
  } // enforce_penalty()

  std::vector<Node*> surface; // closed counter-clockwise cycle of nodes defining the contiguous contact surface
  float mu; // coefficient of friction
  float penalty_stiffness; // contact penalty stiffness
  float damping_viscosity; // contact damping viscosity
  float conductivity; // thermal conductivity
}; // Contact

#endif // CONTACT_H
