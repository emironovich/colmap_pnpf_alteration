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

//todo fix comments
//Super solver

class P35PfEstimator {
 public:
  // The 2D image feature observations.
  typedef Eigen::Vector2d X_t;
  // The observed 3D features in the world frame.
  typedef Eigen::Vector3d Y_t;
  // The transformation from the world to the camera frame.
  // todo fix commment
  struct M_t {double f; Eigen::Matrix3d R; Eigen::Vector3d C;} ;

  // The minimum number of samples needed to estimate a model.
  static const int kMinNumSamples = 4;

  // Estimate the most probable solution of the P3P problem from a set of
  // three 2D-3D point correspondences.
  //
  // @param points2D   Normalized 2D image points as 3x2 matrix.
  // @param points3D   3D world points as 3x3 matrix.
  //
  // @return           Most probable pose as length-1 vector of a 3x4 matrix.
  static std::vector<M_t> Estimate(const std::vector<X_t>& points2D,
                                   const std::vector<Y_t>& points3D);

  // Calculate the squared reprojection error given a set of 2D-3D point
  // correspondences and a projection matrix.
  //
  // @param points2D     Normalized 2D image points as Nx2 matrix.
  // @param points3D     3D world points as Nx3 matrix.
  // @param proj_matrix  3x4 projection matrix.
  // @param residuals    Output vector of residuals.
  static void Residuals(const std::vector<X_t>& points2D,
                        const std::vector<Y_t>& points3D,
                        const M_t& proj_matrix, std::vector<double>* residuals);
};
}

#endif  // COLMAP_ABSOLUTE_POSE_PNPF_H
