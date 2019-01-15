// Include file for the photon map render code

// #ifndef PHOTON_H
// #define PHOTON_H

struct Photon
{
  R3Vector normal;
  R3Point source;
  R3Point position;
  R3Vector direction;
  R3Vector out_direction;
  RNRgb power;
  int bounces;
}; 

bool traceRayDiffuse(R3Scene *scene, RNScalar *prev_ior, R3Ray ray, R3Point *point, R3SceneElement **element, R3Vector *normal, RNScalar termination_rate_ray_trace, RNRgb *power_multiplier);

RNRgb EstimateFlux(R3Kdtree<Photon *> *photon_map,  R3Point point, int num_photons, RNScalar max_distance, RNRgb diffuseBrdf);

void getR3CircleAxes(R3Vector normal, R3Vector *axis1, R3Vector *axis2);
