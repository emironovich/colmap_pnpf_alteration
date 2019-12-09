//
// Created by elizaveta on 30.11.2019.
//

#include "absolute_pose_pnpf.h"
#include "util/logging.h"
//#include "util/types.h"
#include "estimators/utils.h"
#include "pnpf/solvers.hpp"
namespace colmap {
std::vector<P35PfEstimator::M_t> P35PfEstimator::Estimate(
    const std::vector<X_t>& points2D, const std::vector<Y_t>& points3D) {
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

  pnpf::P35PSolver<double> solver;
  int n;
  int max_solutions = pnpf::SolverTraits<pnpf::P35PSolver<double>>::MaxSolutions; //todo ????
  double fs[max_solutions];
  Eigen::Matrix3d Rs[max_solutions];
  Eigen::Vector3d Cs[max_solutions];
  solver.solve(points3D_world, points2D_world, &n, fs, Rs, Cs);

  std::vector<P35PfEstimator::M_t> results(n);
  for (int i = 0; i < n; ++i) {
    results[i].f = fs[i];
    results[i].R = Rs[i];
    results[i].T = -Rs[i]*Cs[i];
  }

  return results;
}

void P35PfEstimator::Residuals(const std::vector<X_t>& points2D,
                               const std::vector<Y_t>& points3D,
                               const M_t& params,
                               std::vector<double>* residuals) {
  // calibration matrix
  Eigen::Matrix3d K;
  K << params.f, 0, 0,
       0, params.f, 0,
       0,        0, 1;

  Eigen::Matrix3x4d proj_matrix;
  proj_matrix << params.R, params.T;
  proj_matrix= K * proj_matrix;
  //todo check
  ComputeSquaredReprojectionError(points2D, points3D, proj_matrix, residuals);
}

}  // namespace colmap
