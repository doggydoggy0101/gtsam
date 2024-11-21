/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    QGMParams.h
 * @brief   The QGMOptimizer class
 * @author  Bang-Shien Chen
 * 
 * modified by
 * @file    GncParams.h
 * @brief   The GncParams class
 * @author  Jingnan Shi
 * @author  Luca Carlone
 * @author  Frank Dellaert
 */

#pragma once

#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>

namespace gtsam {

/* ************************************************************************* */

template<class BaseOptimizerParameters>
class QGMParams {
 public:
  /// For each parameter, specify the corresponding optimizer: e.g., GaussNewtonParams -> GaussNewtonOptimizer.
  typedef typename BaseOptimizerParameters::OptimizerType OptimizerType;

  /// Verbosity levels
  enum Verbosity {
    SILENT = 0,
    SUMMARY,
    WEIGHTS,
    VALUES
  };

  /// Constructor.
  QGMParams(const BaseOptimizerParameters& baseOptimizerParams)
      : baseOptimizerParams(baseOptimizerParams) {
  }

  /// Default constructor.
  QGMParams()
      : baseOptimizerParams() {
  }

  /// QGM parameters.
  BaseOptimizerParameters baseOptimizerParams;  ///< Optimization parameters used to solve the weighted least squares problem at each QGM iteration
  /// any other specific QGM parameters:
  size_t maxIterations = 100;  ///<  Maximum number of iterations
  double relativeCostTol = 1e-5;  ///< If relative cost change is below this threshold, stop iterating
  Verbosity verbosity = SILENT;  ///< Verbosity level

  /// Use IndexVector for inliers and outliers since it is fast
  using IndexVector = FastVector<uint64_t>;
  ///< Slots in the factor graph corresponding to measurements that we know are inliers
  IndexVector knownInliers;
  ///< Slots in the factor graph corresponding to measurements that we know are outliers
  IndexVector knownOutliers;

  /// Set the maximum number of iterations in QGM (changing the max nr of iters might lead to less accurate solutions and is not recommended).
  void setMaxIterations(const size_t maxIter) {
    std::cout
        << "setMaxIterations: changing the max nr of iters might lead to less accurate solutions and is not recommended! "
        << std::endl;
    maxIterations = maxIter;
  }

  /// Set the maximum relative difference to stop iterating.
  void setRelativeCostTol(double value) {
    relativeCostTol = value;
  }

  /// Set the verbosity level.
  void setVerbosityQGM(const Verbosity value) {
    verbosity = value;
  }

  /** (Optional) Provide a vector of measurements that must be considered inliers. The enties in the vector
   * corresponds to the slots in the factor graph. For instance, if you have a nonlinear factor graph nfg,
   * and you provide  knownIn = {0, 2, 15}, QGM will not apply outlier rejection to nfg[0], nfg[2], and nfg[15].
   * This functionality is commonly used in SLAM when one may assume the odometry is outlier free, and
   * only apply QGM to prune outliers from the loop closures.
   * */
  void setKnownInliers(const IndexVector& knownIn) {
    for (size_t i = 0; i < knownIn.size(); i++){
      knownInliers.push_back(knownIn[i]);
    }
    std::sort(knownInliers.begin(), knownInliers.end());
  }

  /** (Optional) Provide a vector of measurements that must be considered outliers. The enties in the vector
   * corresponds to the slots in the factor graph. For instance, if you have a nonlinear factor graph nfg,
   * and you provide  knownOut = {0, 2, 15}, QGM will not apply outlier rejection to nfg[0], nfg[2], and nfg[15].
   * */
  void setKnownOutliers(const IndexVector& knownOut) {
    for (size_t i = 0; i < knownOut.size(); i++){
      knownOutliers.push_back(knownOut[i]);
    }
    std::sort(knownOutliers.begin(), knownOutliers.end());
  }

  /// Equals.
  bool equals(const QGMParams& other, double tol = 1e-9) const {
    return baseOptimizerParams.equals(other.baseOptimizerParams)
        && maxIterations == other.maxIterations
        && verbosity == other.verbosity && knownInliers == other.knownInliers
        && knownOutliers == other.knownOutliers;
  }

  /// Print.
  void print(const std::string& str) const {
    std::cout << str << "\n";
    std::cout << "maxIterations: " << maxIterations << "\n";
    std::cout << "relativeCostTol: " << relativeCostTol << "\n";
    std::cout << "verbosity: " << verbosity << "\n";
    for (size_t i = 0; i < knownInliers.size(); i++)
      std::cout << "knownInliers: " << knownInliers[i] << "\n";
    for (size_t i = 0; i < knownOutliers.size(); i++)
      std::cout << "knownOutliers: " << knownOutliers[i] << "\n";
    baseOptimizerParams.print("Base optimizer params: ");
  }
};

}
