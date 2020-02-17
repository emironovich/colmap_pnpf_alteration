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

#define TEST_NAME "optim/ransac" // todo what is that?????
#include "util/testing.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "base/pose.h"
#include "optim/loransac.h"
#include "optim/loransac_pnpf.h"
#include "util/random.h"
#include "estimators/absolute_pose_pnpf.h"
#include "estimators/absolute_pose.h"
#include <random>
#include <unsupported/Eigen/MatrixFunctions>

using namespace colmap;

BOOST_AUTO_TEST_CASE(TestP35PLORANSAC) {
  SetPRNGSeed(0);

  const size_t num_samples = 1000;
  const size_t num_outliers = 400;

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

  // Add some faulty data.
  for (size_t i = 0; i < num_outliers; ++i) {
    points2D[i] = Eigen::Vector2d(RandomReal(-3000.0, -2000.0),
                                   RandomReal(-4000.0, -3000.0));
  }

  // Robustly estimate transformation using RANSAC.
  RANSACOptions options;
  options.max_error = 10;
  LORANSAC<P35PfEstimator, EPNPEstimator> ransac(options);
  const auto report = ransac.Estimate(points2D, points3D);

  BOOST_CHECK_EQUAL(report.success, true);
  BOOST_CHECK_GT(report.num_trials, 0);

  // Make sure outliers were detected correctly.
  BOOST_CHECK_EQUAL(report.support.num_inliers, num_samples - num_outliers);
  for (size_t i = 0; i < num_samples; ++i) {
    if (i < num_outliers) {
      BOOST_CHECK(!report.inlier_mask[i]);
    } else {
      BOOST_CHECK(report.inlier_mask[i]);
    }
  }

  // Make sure original parameters were estimated correctly.

  //Test if correct focal length has been found.
  const double f_diff = abs(f - (report.model).f) / f;
  BOOST_CHECK(f_diff < 1e-6);

  // Test if correct rotation has been determined.
  const double R_diff = (R - (report.model).R).norm()/ 3;
  BOOST_CHECK(R_diff < 1e-4);

  // Test if correct translation has been determined.
  const double T_diff = (-R*C - (report.model).T).norm();
  BOOST_CHECK(T_diff < 1e-3);
}
