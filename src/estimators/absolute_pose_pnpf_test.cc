// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#define TEST_NAME "base/absolute_pose_pnpf"
#include "util/testing.h"

#include <Eigen/Core>

#include "base/pose.h"
#include "optim/ransac.h"
#include "util/random.h"
#include "estimators/absolute_pose_pnpf.h"
#include <iostream>
using namespace colmap;

BOOST_AUTO_TEST_CASE(TestP35Pf) {
  const size_t num_samples = 10;
  // Create some arbitrary transformation.
  std::random_device dev;
  std::mt19937_64 generator(dev()); // Mersenne Twister generator
  std::uniform_real_distribution<double> uniformDistribution(-1., 1.);
  auto uniform = [&]() { return uniformDistribution(generator); };

  // focal distance
  double f = 200 + 1800 * (uniformDistribution(generator) + 1) / 2;

  // rotation
  Eigen::Vector4d rQuat =  Eigen::Vector4d::NullaryExpr(4, 1, uniform);
  Eigen::Matrix3d R =  QuaternionToRotationMatrix(rQuat);

  // camera center
  Eigen::Vector3d C = (Eigen::Vector3d::NullaryExpr(3, 1, uniform)) / 2;

  // calibration matrix
  Eigen::Matrix3d K;
  K << f, 0, 0,
      0, f, 0,
      0, 0, 1;

  // projection matrix
  Eigen::Matrix3x4d P;
  P.setIdentity();
  P.col(3) = -C;
  P = K * R * P;

  // Generate exact data.
  std::vector<Eigen::Vector3d> points3D;
  std::vector<Eigen::Vector2d> points2D;
  for (size_t i = 0; i < num_samples; ++i) {
    Eigen::Vector3d pnt;
    pnt << 0, 0, 6;
    pnt += 2*(Eigen::Vector3d::NullaryExpr(3, 1, uniform));
    points3D.emplace_back(R.transpose()*pnt + C);
    points2D.emplace_back((P*(points3D.back()).homogeneous()).hnormalized());
  }

  // Robustly estimate transformation using RANSAC.
  RANSACOptions options;
  options.max_error = 10;
  RANSAC<P35PfEstimator> ransac(options);
  const auto report = ransac.Estimate(points2D, points3D);

  BOOST_CHECK_EQUAL(report.success, true);
  BOOST_CHECK_GT(report.num_trials, 0);

  // Make sure original parameters were estimated correctly.

  //Test if correct focal length has been found.
  const double f_diff = abs(f - (report.model).f) / f;
  BOOST_CHECK(f_diff < 1e-2);

  // Test if correct rotation has been determined.
  const double R_diff = (R - (report.model).R).norm()/ 3;
  BOOST_CHECK(R_diff < 1e-2);

  // Test if correct translation has been determined.
  const double T_diff = (-R*C - (report.model).T).norm();
  BOOST_CHECK(T_diff < 1e-1);

  auto points3D_faulty = points3D;
  for (size_t i = 0; i < points3D.size(); ++i) {
    points3D_faulty[i](0) = 20;
  }

  // Test residuals of exact points.
  std::vector<double> residuals;
  P35PfEstimator::Residuals(points2D, points3D, report.model, &residuals);
  for (size_t i = 0; i < residuals.size(); ++i) {
    BOOST_CHECK(residuals[i] < 1e-3);
  }

  // Test residuals of faulty points.
  P35PfEstimator::Residuals(points2D, points3D_faulty, report.model,
                            &residuals);
  for (size_t i = 0; i < residuals.size(); ++i) {
    BOOST_CHECK(residuals[i] > 0.1);
  }
}