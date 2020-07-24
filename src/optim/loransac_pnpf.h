//
// Created by elizaveta on 01.12.2019.
//

#ifndef COLMAP_LORANSAC_PNPF_H
#define COLMAP_LORANSAC_PNPF_H

#include <Eigen/Dense>
#include <cfloat>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "estimators/absolute_pose_pnpf.h"
#include "optim/random_sampler.h"
#include "optim/ransac.h"
#include "optim/support_measurement.h"
#include "util/alignment.h"
#include "util/logging.h"

namespace colmap {

// Implementation of LO-RANSAC (Locally Optimized RANSAC)
// that uses PnPf solvers with EPNP Solver
//
// "Locally Optimized RANSAC" Ondrej Chum, Jiri Matas, Josef Kittler, DAGM 2003.
template <typename Estimator, typename LocalEstimator,
          typename SupportMeasurer = InlierSupportMeasurer,
          typename Sampler = RandomSampler>
class LORANSAC_PnPf : public RANSAC<Estimator, SupportMeasurer, Sampler> {
 public:
  using typename RANSAC<Estimator, SupportMeasurer, Sampler>::Report;

  explicit LORANSAC_PnPf(const RANSACOptions& options);

  // Robustly estimate model with RANSAC (RANdom SAmple Consensus).
  //
  // @param X              Independent variables.
  // @param Y              Dependent variables.
  //
  // @return               The report with the results of the estimation.
  Report Estimate(const std::vector<typename Estimator::X_t>& X,
                  const std::vector<typename Estimator::Y_t>& Y,
                  double half_diag);

  // Objects used in RANSAC procedure.
  using RANSAC<Estimator, SupportMeasurer, Sampler>::estimator;
  LocalEstimator local_estimator;
  using RANSAC<Estimator, SupportMeasurer, Sampler>::sampler;
  using RANSAC<Estimator, SupportMeasurer, Sampler>::support_measurer;

 private:
  using RANSAC<Estimator, SupportMeasurer, Sampler>::options_;
};

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

template <typename Estimator, typename LocalEstimator, typename SupportMeasurer,
          typename Sampler>
LORANSAC_PnPf<Estimator, LocalEstimator, SupportMeasurer,
              Sampler>::LORANSAC_PnPf(const RANSACOptions& options)
    : RANSAC<Estimator, SupportMeasurer, Sampler>(options) {}

template <typename Estimator, typename LocalEstimator, typename SupportMeasurer,
          typename Sampler>
typename LORANSAC_PnPf<Estimator, LocalEstimator, SupportMeasurer,
                       Sampler>::Report
LORANSAC_PnPf<Estimator, LocalEstimator, SupportMeasurer, Sampler>::Estimate(
    const std::vector<typename Estimator::X_t>& X,
    const std::vector<typename Estimator::Y_t>& Y, double half_diag) {
  CHECK_EQ(X.size(), Y.size());
  const size_t num_samples = X.size();

  typename RANSAC<Estimator, SupportMeasurer, Sampler>::Report report;
  report.success = false;
  report.num_trials = 0;

  if (num_samples < Estimator::kMinNumSamples) {
    return report;
  }
  //
  typename SupportMeasurer::Support best_support;
  typename Estimator::M_t best_model;
  bool best_model_is_local = false;

  bool abort = false;

  const double max_residual = options_.max_error * options_.max_error;

  std::vector<double> residuals(num_samples);

  std::vector<typename LocalEstimator::X_t> X_inlier;
  std::vector<typename LocalEstimator::Y_t> Y_inlier;

  std::vector<typename Estimator::X_t> X_rand(Estimator::kMinNumSamples);
  std::vector<typename Estimator::Y_t> Y_rand(Estimator::kMinNumSamples);

  std::vector<typename LocalEstimator::X_t> X_normalized(num_samples);

  sampler.Initialize(num_samples);

  size_t max_num_trials = options_.max_num_trials;
  max_num_trials = std::min<size_t>(max_num_trials, sampler.MaxNumSamples());
  size_t dyn_max_num_trials = max_num_trials;

  for (report.num_trials = 0; report.num_trials < max_num_trials;
       ++report.num_trials) {
    if (abort) {
      report.num_trials += 1;
      break;
    }

    sampler.SampleXY(X, Y, &X_rand, &Y_rand);

    // Estimate model for current subset.
    const std::vector<typename Estimator::M_t> sample_models =
        estimator.Estimate(X_rand, Y_rand, half_diag);

    // Iterate through all estimated models
    for (const auto& sample_model : sample_models) {
      estimator.Residuals(X, Y, sample_model, &residuals);
      CHECK_EQ(residuals.size(), X.size());

      const auto support = support_measurer.Evaluate(residuals, max_residual);

      // Do local optimization if better than all previous subsets.
      if (support_measurer.Compare(support, best_support)) {
        best_support = support;
        best_model = sample_model;
        best_model_is_local = false;

        // Estimate locally optimized model from inliers.
        if (support.num_inliers > Estimator::kMinNumSamples &&
            support.num_inliers >= LocalEstimator::kMinNumSamples) {
          X_inlier.clear();
          Y_inlier.clear();
          X_inlier.reserve(support.num_inliers);
          Y_inlier.reserve(support.num_inliers);
          for (size_t i = 0; i < residuals.size(); ++i) {
            X_normalized[i] = X[i] / best_model.f;
            if (residuals[i] <= max_residual) {
              X_inlier.push_back(X_normalized[i]);
              Y_inlier.push_back(Y[i]);
            }
          }

          const std::vector<typename LocalEstimator::M_t> local_models =
              local_estimator.Estimate(X_inlier, Y_inlier);

          for (const auto& local_model : local_models) {
            typename Estimator::M_t local_model_rt_format;
            local_model_rt_format.f = best_model.f;
            local_model_rt_format.R = local_model.leftCols(3);
            local_model_rt_format.T = local_model.rightCols(1);

            estimator.Residuals(X, Y, local_model_rt_format, &residuals);
            CHECK_EQ(residuals.size(), X.size());

            const auto local_support =
                support_measurer.Evaluate(residuals, max_residual);

            // Check if non-locally optimized model is better.
            if (support_measurer.Compare(local_support, best_support)) {
              best_support = local_support;
              best_model.R = local_model.leftCols(3);
              best_model.T = local_model.rightCols(1);
              best_model_is_local = true;
            }
          }
        }

        dyn_max_num_trials =
            RANSAC<Estimator, SupportMeasurer, Sampler>::ComputeNumTrials(
                best_support.num_inliers, num_samples, options_.confidence,
                options_.dyn_num_trials_multiplier);
      }

      if (report.num_trials >= dyn_max_num_trials &&
          report.num_trials >= options_.min_num_trials) {
        abort = true;
        break;
      }
    }
  }

  report.support = best_support;
  report.model = best_model;

  // No valid model was found
  if (report.support.num_inliers < estimator.kMinNumSamples) {
    return report;
  }

  report.success = true;

  // Determine inlier mask. Note that this calculates the residuals for the
  // best model twice, but saves to copy and fill the inlier mask for each
  // evaluated model. Some benchmarking revealed that this approach is faster.
  if (best_model_is_local) {
    typename LocalEstimator::M_t best_local_model;
    best_local_model << (report.model).R, (report.model).T;

    for (size_t i = 0; i < X.size(); ++i) {
      X_normalized[i] = X[i] / (report.model).f;
    }

    local_estimator.Residuals(X_normalized, Y, best_local_model, &residuals);
  } else {
    estimator.Residuals(X, Y, report.model, &residuals);
  }

  CHECK_EQ(residuals.size(), X.size());

  report.inlier_mask.resize(num_samples);
  for (size_t i = 0; i < residuals.size(); ++i) {
    if (residuals[i] <= max_residual) {
      report.inlier_mask[i] = true;
    } else {
      report.inlier_mask[i] = false;
    }
  }

  return report;
}

}  // namespace colmap

#endif  // COLMAP_LORANSAC_PNPF_H
