// Include file for the photon map render code

// #ifndef RENDER_H
// #define RENDER_H

#include "photonmap.h"



R2Image *RenderImage(R3Scene *scene, R3Kdtree<Photon *> *photon_map, R3Kdtree<Photon *> *caustic_map, int width, int height, int print_verbose, int num_samples, RNScalar general_search_range, RNScalar caustic_search_range, RNScalar tone_map_const, int num_photon_estimate);

