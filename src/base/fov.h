//
// Created by elizaveta on 24.07.2020.
//

#ifndef COLMAP_FOV_H
#define COLMAP_FOV_H

#include "util/logging.h"
#include "base/camera.h"

namespace colmap {

struct FOVOptions {
  // Whether to check FOV for restrictions
  bool check_fov = false;

  int min_degrees = -1;
  int max_degrees = -1;

  FOVOptions(){};
  FOVOptions(int min_degrees, int max_degrees)
      : check_fov(true), min_degrees(min_degrees), max_degrees(max_degrees){};

  void Check() const {
    if (check_fov) {
      CHECK_GT(max_degrees, min_degrees);
      CHECK_GE(min_degrees, 0);
      CHECK_LE(max_degrees, 360);
    }
  }

  // Returns true if check_fov = false, else checks
  // min_degrees < FOV < max_degrees
  bool CorrectFOV(double fov) const {
    return !check_fov || (min_degrees < fov && fov < max_degrees);
  }
};

double FindFOV(double half_dimention, double focal_distance);
double FOV2FocalLength(double half_dimention, double fov);
double HalfDiag(const Camera& camera);

}  // namespace colmap

#endif  // COLMAP_FOV_H
