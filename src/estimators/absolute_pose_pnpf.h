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

#include "absolute_pose_pnpf.h"
#include "util/logging.h"
//#include "util/types.h"
#include "estimators/utils.h"
#include "pnpf/solvers.hpp"

#include <utility>

namespace colmap {

template <class Estimator>
class PnPfEstimator {
 public:
  // The 2D image feature observations.
  typedef Eigen::Vector2d X_t;
  // The observed 3D features in the world frame.
  typedef Eigen::Vector3d Y_t;
  // Estimated camera parameters: focal length, rotation and translation.
  struct M_t {
    double f;
    Eigen::Matrix3d R;
    Eigen::Vector3d T;
  };

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
                                   const std::vector<Y_t>& points3D) {
    CHECK_EQ(points2D.size(), 4);
    CHECK_EQ(points3D.size(), 4);

    Eigen::Matrix3x4d points3D_world;
    points3D_world.col(0) = points3D[0];
    points3D_world.col(1) = points3D[1];
    points3D_world.col(2) = points3D[2];
    points3D_world.col(3) = points3D[3];

    Eigen::Matrix<double, 2, 4> points2D_world;
    points2D_world.col(0) = points2D[0];
    points2D_world.col(1) = points2D[1];
    points2D_world.col(2) = points2D[2];
    points2D_world.col(3) = points2D[3];

    int n;
    int max_solutions = Estimator::getMaxSolutions();  // todo ????
    double fs[max_solutions];
    Eigen::Matrix3d Rs[max_solutions];
    Eigen::Vector3d Cs[max_solutions];
    Estimator::solve(points3D_world, points2D_world, &n, fs, Rs, Cs);

    std::vector<M_t> results(n);
    for (int i = 0; i < n; ++i) {
      results[i].f = fs[i];
      results[i].R = Rs[i];
      results[i].T = -Rs[i] * Cs[i];
    }

    return results;
  }

  // Calculate the squared reprojection error given a set of 2D-3D point
  // correspondences and a parameters set.
  //
  // @param points2D     2D image points as Nx2 matrix.
  // @param points3D     3D world points as Nx3 matrix.
  // @param params       {f, R, T}  focal distance, rotation and translation
  // @param residuals    Output vector of residuals.
  static void Residuals(const std::vector<X_t>& points2D,
                        const std::vector<Y_t>& points3D, const M_t& params,
                        std::vector<double>* residuals) {
    // calibration matrix
    Eigen::Matrix3d K;
    K << params.f, 0, 0, 0, params.f, 0, 0, 0, 1;

    Eigen::Matrix3x4d proj_matrix;
    proj_matrix << params.R, params.T;
    proj_matrix = K * proj_matrix;
    // todo check
    ComputeSquaredReprojectionError(points2D, points3D, proj_matrix, residuals);
  }
};

// Solver for the P3.5Pf (Perspective-3.5-Point with unknown focal
// length) problem.
//
// The algorithm is based on the following paper:
//
//    Changchang Wu. P3.5P: Pose Estimation with Unknown Focal Length.
class P35PfEstimator : public PnPfEstimator<P35PfEstimator> {
 public:
  static void solve(const Eigen::Matrix3x4d& points3D_world,
                    const Eigen::Matrix<double, 2, 4>& points2D_world, int* n,
                    double* fs, Eigen::Matrix3d* Rs, Eigen::Vector3d* Cs) {
    pnpf::P35PSolver<double> solver;
    solver.solve(points3D_world, points2D_world, n, fs, Rs, Cs);
  }
  static int getMaxSolutions() {
    return pnpf::SolverTraits<pnpf::P35PSolver<double>>::MaxSolutions;
  }
};

class P4PfEstimator : public PnPfEstimator<P4PfEstimator> {
 public:
  static void solve(const Eigen::Matrix3x4d& points3D_world,
                    const Eigen::Matrix<double, 2, 4>& points2D_world, int* n,
                    double* fs, Eigen::Matrix3d* Rs, Eigen::Vector3d* Cs) {
    pnpf::P4PSolver<double> solver;
    solver.solve(points3D_world, points2D_world, n, fs, Rs, Cs);
  }
  static int getMaxSolutions() {
    return pnpf::SolverTraits<pnpf::P4PSolver<double>>::MaxSolutions;
  }
};
// todo: add description

}  // namespace colmap

#endif  // COLMAP_ABSOLUTE_POSE_PNPF_H
