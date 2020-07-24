#include "base/fov.h"
#include <cmath>

namespace colmap {

double FindFOV(double half_dimention, double focal_distance) {
  const double pi = 3.14159265;
  return 2 * atan(half_dimention / focal_distance) * 180 / pi;
}

double FOV2FocalLength(double half_dimention, double fov_degrees) {
  const double pi = 3.14159265;
  double fov_radian =  fov_degrees * pi / 180;
  return  half_dimention / tan(fov_radian / 2);
}

double HalfDiag(const Camera& camera) {
  return sqrt(pow(camera.Width(), 2) + pow(camera.Height(), 2)) / 2;
}

}  // namespace colmap