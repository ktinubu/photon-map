// Source file for photonmap renderer
// This implementation is a simple raycaster.
// Replace it with your own code.



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Graphics/R3Graphics.h"
#include <iostream>
#include "photonmap.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <string> 
#include <chrono>




// program variables
static RNScalar termination_rate = 0.001;
static RNScalar camera_index_of_refraction = 1.0;


////////////////////////////////////////////////////////////////////////
// Function to render image with photon mapping
////////////////////////////////////////////////////////////////////////

// rgb to grayscale consts
R3Matrix gray_conv_matrix = R3Matrix(1.0, 0.956, 0.621, 1.0, -0.272, -0.647, 1.0, -1.106, 1.703).Inverse();
R3Vector gray_conv_coeffs =  R3Vector(gray_conv_matrix[0][0], gray_conv_matrix[0][1], gray_conv_matrix[0][2]);

R2Image *
RenderImage(R3Scene *scene,
  R3Kdtree<Photon *> *photon_map,
  R3Kdtree<Photon *> *caustic_map,
  int width,
  int height,
  int print_verbose,
  int num_samples,
  RNScalar max_estimate_dist_proportion_global,
  RNScalar max_estimate_dist_proportion_caustic,
  RNScalar reinhard_tone_map_a,
  int num_photon_estimate)

{
  assert(photon_map);
  assert(caustic_map);

  // Start statistics
  RNTime start_time;
  start_time.Read();
  int ray_count = 0;
  // Allocate image
  R2Image *image = new R2Image(width, height);
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    return NULL;
  }
  R3SceneElement *element;
  R3Point point;
  R3Vector normal;
  RNScalar roulette_multiplier = RNScalar(1)/ 1 - termination_rate; // becuase of russian roulette
  int num_rendered_pixels = 0;
  int total_pixels = width * height;
  std::vector<RNRgb> pixels;

  // precompute axes for area lights
  R3Vector axes1[scene->NLights()];
  R3Vector axes2[scene->NLights()];
  for (int k = 0; k < scene->NLights(); k++) {
    R3Light *light = scene->Light(k);
    if (light->ClassID() == R3AreaLight::CLASS_ID()) {
      R3AreaLight *area_light = (R3AreaLight *) light;
      getR3CircleAxes(area_light->Direction(), &axes1[k], &axes2[k]);
    }
  }
  // Draw intersection point and normal for some rays
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      RNRgb color = RNRgb(0,0,0);

      for (int s = 0; s < num_samples; s ++) {
        R3Ray ray = scene->Viewer().WorldRay(i, j);
        // std::cout<<ray.Point(0)[0]<< ", " << ray.Point(0)[1] << ", " << ray.Point(0)[2] <<std::endl;

        RNScalar prev_ior = camera_index_of_refraction;
        RNRgb power_multiplier =  RNRgb(1,1,1);
        if (!traceRayDiffuse(scene, &prev_ior, ray, &point, &element, &normal, termination_rate, &power_multiplier)) {
          continue;
        }
        normal.Normalize();
        // add indirect lighting contribution with photon map

        const R3Material *material = (element) ? element->Material() : &R3default_material;
        const R3Brdf *brdf = (material) ? material->Brdf() : &R3default_brdf;


        //power_multiplier =  brdf->Diffuse();
        //std::cout<<power_multiplier[0]<< ", " << power_multiplier[1] << ", " << power_multiplier[2] <<std::endl;

        const RNRgb& diff_brdf = power_multiplier / RN_PI;
        color += roulette_multiplier * EstimateFlux(photon_map, point, num_photon_estimate, max_estimate_dist_proportion_global * scene->BBox().DiagonalRadius(), diff_brdf);
        // add caustic contribution
        color += roulette_multiplier * EstimateFlux(caustic_map, point, num_photon_estimate, max_estimate_dist_proportion_caustic * scene->BBox().DiagonalRadius(), diff_brdf);
        // add emitted light
        color += brdf->Emission();

        // add direct light contribution with path tracing
        for (int k = 0; k < scene->NLights(); k++) {
        R3Light *light = scene->Light(k);
          if (light->ClassID() == R3PointLight::CLASS_ID()) {
            R3PointLight *point_light = (R3PointLight *) light;
            R3SceneElement *intersecting_element = element;
            if (scene->Intersects(R3Ray(point_light->Position(), point), NULL, &element, NULL, NULL, NULL, NULL) && element == intersecting_element) {              
              const RNRgb& Ic = point_light->Color() / (R3SquaredDistance(point_light->Position(), point));
              R3Vector L = point_light->DirectionFromPoint(point);
              RNScalar NL = normal.Dot(L);
              if (RNIsNegativeOrZero(NL)) {
                continue;
              }
              // color += roulette_multiplier * NL * diff_brdf * Ic;
            }
          } else if (light->ClassID() == R3AreaLight::CLASS_ID()) {
            // monte carlo estimate to evaluate direct illumination
              // as desdcribed in https://www.cs.utah.edu/~shirley/papers/rw91.pdf
            R3AreaLight *area_light = (R3AreaLight *) light;
            R3Point source_pos;
            RNScalar r1;
            RNScalar r2;
            do {
              r1 = (RNRandomScalar() * 2) - 1;
              r2 = (RNRandomScalar() * 2) - 1;
            } while(r1 * r1 + r2 * r2 > 1);
            source_pos = area_light->Position();
            source_pos += (r1 * axes1[k] * area_light->Radius()) + (r2 * axes2[k] * area_light->Radius());
            source_pos += area_light->Direction() * RN_EPSILON;
            R3SceneElement *intersecting_element = element;
            R3Vector light_to_point = point - source_pos;
            light_to_point.Normalize();
            RNScalar cos_light = light_to_point.Dot(area_light->Direction());
            RNScalar pdf = RNScalar(1)/ (RN_PI * area_light->Radius() * area_light->Radius());
            if (RNIsNegativeOrZero(cos_light)) {
              continue;
            }
            if (scene->Intersects(R3Ray(source_pos, point), NULL, &element, NULL, NULL, NULL, NULL) && element == intersecting_element) {              
              const RNRgb& Ic = area_light->Color() / R3SquaredDistance(source_pos, point);
              RNScalar cos_point = normal.Dot(-light_to_point);
              color += roulette_multiplier * diff_brdf * Ic * cos_point * cos_light / pdf;
            }
          } else if (light->ClassID() == R3DirectionalLight::CLASS_ID()) {
            R3DirectionalLight *dir_light = (R3DirectionalLight *) light;
            R3Vector dir_light_dir = dir_light->Direction();
            dir_light_dir.Normalize();
            R3SceneElement *intersecting_element = element;
            if (scene->Intersects(R3Ray(point - (dir_light_dir * 2 *scene->BBox().DiagonalRadius()), point), NULL, &element, NULL, NULL, NULL, NULL) && element == intersecting_element) {
              const RNRgb& Ic = dir_light->Color();
              RNScalar NL = normal.Dot(-dir_light->Direction());
              if (RNIsNegativeOrZero(NL)) {
                continue;
              }
              color += roulette_multiplier * NL * diff_brdf * Ic;
            }
          } else if (light->ClassID() == R3SpotLight::CLASS_ID()) {
            R3SpotLight *spot_light = (R3SpotLight *) light;
            R3SceneElement *intersecting_element = element;
            R3Vector central_direction = spot_light->Direction();
            central_direction.Normalize();
            if (normal.Dot(central_direction) < cos(spot_light->CutOffAngle()) && scene->Intersects(R3Ray(spot_light->Position(), point), NULL, &element, NULL, NULL, NULL, NULL) && element == intersecting_element) {              
              const RNRgb& Ic = spot_light->Color() / (R3SquaredDistance(spot_light->Position(), point));
              R3Vector L = spot_light->DirectionFromPoint(point);
              RNScalar NL = normal.Dot(L);
              if (RNIsNegativeOrZero(NL)) {
                continue;
              }
              color += roulette_multiplier * NL * diff_brdf * Ic;
            }
          } else {
            std::cout<<"unrecognized light"<<std::endl;
            assert(false);
          }
        }
      }
      color /= num_samples;   

      pixels.push_back(color);
      num_rendered_pixels++;
    }
    std::cout<<double(num_rendered_pixels * 100) / total_pixels<<"% of pixels rendered"<<std::endl;
  }

  // Print statistics
  if (print_verbose) {
    printf("Rendered image ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Rays = %d\n", ray_count);
    fflush(stdout);
  }

  std::string width_str = std::to_string(width);
  std::string height_str = std::to_string(height);
  std::string now = std::to_string(std::chrono::time_point_cast<std::chrono::milliseconds>(
    std::chrono::system_clock::now()).time_since_epoch().count());
  std::string grey_conv = std::to_string(double(gray_conv_coeffs[0])) +  "=" + std::to_string(double(gray_conv_coeffs[1])) + "=" + std::to_string(double(gray_conv_coeffs[2]));
  std::string tone_map = std::to_string(reinhard_tone_map_a);
  std::ofstream f;
  f.open(width_str + ":" + height_str + ":"  + grey_conv + ":" + tone_map + ":"  + now + ".csv");
  f << "pixel number,r,g,b\n";
  int pix_count = 0;
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      f << std::to_string(double(pixels[pix_count][0])) + "," + std::to_string(pixels[pix_count][1]) + "," + std::to_string(pixels[pix_count][2]) + "\n";
      pix_count++;
    }
  } 
  std::cout<<"applying tone mapping..."<<std::endl;
  // apply tone mapping from Reinhard '02

  pix_count = 0;
  // get avg luminance
  RNScalar avg_lum = 0.0;
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      RNRgb color = pixels[pix_count];
      RNScalar curr_lum = gray_conv_coeffs.Dot(R3Vector(color[0], color[1], color[2]));
      avg_lum += log(RN_EPSILON + curr_lum);
      pix_count++;
    }
  }

  avg_lum /= (height*width);
  avg_lum = exp(avg_lum);

  // get maximum scaled luminance
  RNScalar max_t_pix = 0.0;
  pix_count = 0;
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      RNRgb color = pixels[pix_count];
      RNScalar curr_lum = gray_conv_coeffs.Dot(R3Vector(color[0], color[1], color[2]));
      RNScalar t_pix = reinhard_tone_map_a * curr_lum / avg_lum;
      if (t_pix > max_t_pix) {
        max_t_pix = t_pix;
      }
      pix_count++;
    }
  }

  // apply operator
  RNScalar global_max_color = 0.0;
  pix_count = 0;
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      RNRgb color = pixels[pix_count];
      RNScalar curr_lum = gray_conv_coeffs.Dot(R3Vector(color[0], color[1], color[2]));
      RNScalar scaling = 1;
      RNScalar tone_map = 0;
      if (RNIsPositive(curr_lum)) {
        RNScalar t_pix = reinhard_tone_map_a * curr_lum / avg_lum;
        tone_map = t_pix * (1 + (t_pix / pow(max_t_pix, 2))) / (1 + t_pix);
        scaling = tone_map / curr_lum; 
      }
      color *= scaling;
      // get_max_color
      RNScalar max_color =  std::max({color[0],color[1], color[2]});      
      if (max_color > global_max_color) {
        global_max_color = max_color;
      }
      pixels[pix_count] = color;
      pix_count++;
    }
  }

  // normalize pixel values
  pix_count = 0;
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      RNRgb color = pixels[pix_count];
      color /= global_max_color;
      pixels[pix_count] = color;
      std::cout<<color[0]<< ", " << color[1] << ", " << color[2] <<std::endl;

      image->SetPixelRGB(i, j, color);
      pix_count++;
    }
  }
  return image;
}

  



