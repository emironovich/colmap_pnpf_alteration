//
// Created by elizaveta on 30.11.2019.
//

#ifndef COLMAP_ABSOLUTE_POSE_PNPF_H
#define COLMAP_ABSOLUTE_POSE_PNPF_H

#include <array>
#include <vector>

#include <Eigen/Core>

#include "util/alignment.h"
#include "util/types.h"

#include <utility>

namespace colmap {

// Solver for the P3.5Pf (Perspective-3.5-Point with unknown focal
// length) problem.
//
// The algorithm is based on the following paper:
//
//    Changchang Wu. P3.5P: Pose Estimation with Unknown Focal Length.
class P35PfEstimator {
 public:
  // The 2D image feature observations.
  typedef Eigen::Vector2d X_t;
  // The observed 3D features in the world frame.
  typedef Eigen::Vector3d Y_t;
  // Estimated camera parameters: focal length, rotation and translation.
  struct M_t {double f; Eigen::Matrix3d R; Eigen::Vector3d T;} ;

  // The minimum number of samples needed to estimate a model.
  static const int kMinNumSamples = 4;

  // Estimate solutions of the P3.5Pf problem from a set of four 2D-3D point
  // correspondences.
  //
  // The algorithm is applicable to pinhole-like cameras with the calibration
  // matrix of the form K = diag(f, f, 1), i.e. camera with no distortion,
  // principal point at (0, 0) and focal length equalling to f.
  //
  // @param points2D   Normalized 2D image points as 4x2 matrix.
  //  // @param points3D   3D world points as 4x3 matrix.
  //
  // @return           Vector of solutions for estimated parameters: {f, R, T}.
  //                   Number of solutions is not greater than 10.
  static std::vector<M_t> Estimate(const std::vector<X_t>& points2D,
                                   const std::vector<Y_t>& points3D);

  // Calculate the squared reprojection error given a set of 2D-3D point
  // correspondences and a parameters set.
  //
  // @param points2D     2D image points as Nx2 matrix.
  // @param points3D     3D world points as Nx3 matrix.
  // @param params       {f, R, T}  focal distance, rotation and translation
  // @param residuals    Output vector of residuals.
  static void Residuals(const std::vector<X_t>& points2D,
                        const std::vector<Y_t>& points3D,
                        const M_t& params, std::vector<double>* residuals);
};
}

#endif  // COLMAP_ABSOLUTE_POSE_PNPF_H
