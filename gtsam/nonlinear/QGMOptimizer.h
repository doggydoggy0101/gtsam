/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    QGMOptimizer.h
 * @brief   The QGMOptimizer class
 * @author  Bang-Shien Chen
 * 
 * modified by
 * @file    GncOptimizer.h
 * @brief   The GncOptimizer class
 * @author  Jingnan Shi
 * @author  Luca Carlone
 * @author  Frank Dellaert
 */

#pragma once

#include <gtsam/nonlinear/QGMParams.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/internal/ChiSquaredInverse.h>

namespace gtsam {
/*
 * Quantile of chi-squared distribution with given degrees of freedom at probability alpha.
 * Equivalent to chi2inv in Matlab.
 */
static double Chi2inv(const double alpha, const size_t dofs) {
  return internal::chi_squared_quantile(dofs, alpha);
}

/* ************************************************************************* */
template<class QGMParameters>
class QGMOptimizer {
 public:
  /// For each parameter, specify the corresponding optimizer: e.g., GaussNewtonParams -> GaussNewtonOptimizer.
  typedef typename QGMParameters::OptimizerType BaseOptimizer;

 private:
  NonlinearFactorGraph nfg_; ///< Original factor graph to be solved by QGM.
  Values state_; ///< Initial values to be used at each iteration by QGM.
  QGMParameters params_; ///< QGM parameters.
  Vector weights_;  ///< Weights associated to each factor in QGM (this could be a local variable in optimize, but it is useful to make it accessible from outside).
  Vector barcSq_;  ///< Inlier thresholds. A factor is considered an inlier if factor.error() < barcSq_[i] (where i is the position of the factor in the factor graph. Note that factor.error() whitens by the covariance.

 public:
  /// Constructor.
  QGMOptimizer(const NonlinearFactorGraph& graph, const Values& initialValues,
               const QGMParameters& params = QGMParameters())
      : state_(initialValues),
        params_(params) {

    // make sure all noiseModels are Gaussian or convert to Gaussian
    nfg_.resize(graph.size());
    for (size_t i = 0; i < graph.size(); i++) {
      if (graph[i]) {
        NoiseModelFactor::shared_ptr factor = graph.at<NoiseModelFactor>(i);
        auto robust =
            std::dynamic_pointer_cast<noiseModel::Robust>(factor->noiseModel());
        // if the factor has a robust loss, we remove the robust loss
        nfg_[i] = robust ? factor-> cloneWithNewNoiseModel(robust->noise()) : factor;
      }
    }

    // check that known inliers and outliers make sense:
    std::vector<size_t> inconsistentlySpecifiedWeights; // measurements the user has incorrectly specified
    // to be BOTH known inliers and known outliers
    std::set_intersection(params.knownInliers.begin(),params.knownInliers.end(),
                        params.knownOutliers.begin(),params.knownOutliers.end(),
                        std::inserter(inconsistentlySpecifiedWeights, inconsistentlySpecifiedWeights.begin()));
    if(inconsistentlySpecifiedWeights.size() > 0){ // if we have inconsistently specified weights, we throw an exception
      params.print("params\n");
      throw std::runtime_error("QGMOptimizer::constructor: the user has selected one or more measurements"
          " to be BOTH a known inlier and a known outlier.");
    }
    // check that known inliers are in the graph
    for (size_t i = 0; i < params.knownInliers.size(); i++){
      if( params.knownInliers[i] > nfg_.size()-1 ){ // outside graph
        throw std::runtime_error("QGMOptimizer::constructor: the user has selected one or more measurements"
                  "that are not in the factor graph to be known inliers.");
      }
    }
    // check that known outliers are in the graph
    for (size_t i = 0; i < params.knownOutliers.size(); i++){
      if( params.knownOutliers[i] > nfg_.size()-1 ){ // outside graph
        throw std::runtime_error("QGMOptimizer::constructor: the user has selected one or more measurements"
                  "that are not in the factor graph to be known outliers.");
      }
    }
    // initialize weights (if we don't have prior knowledge of inliers/outliers
    // the weights are all initialized to 1.
    weights_ = initializeWeightsFromKnownInliersAndOutliers();

    // set default barcSq_ (inlier threshold)
    double alpha = 0.99; // with this (default) probability, inlier residuals are smaller than barcSq_
    setInlierCostThresholdsAtProbability(alpha);
  }

  /** Set the maximum weighted residual error for an inlier (same for all factors). For a factor in the form f(x) = 0.5 * || r(x) ||^2_Omega,
   * the inlier threshold is the largest value of f(x) for the corresponding measurement to be considered an inlier.
   * In other words, an inlier at x is such that 0.5 * || r(x) ||^2_Omega <= barcSq.
   * Assuming an isotropic measurement covariance sigma^2 * Identity, the cost becomes: 0.5 * 1/sigma^2 || r(x) ||^2 <= barcSq.
   * Hence || r(x) ||^2 <= 2 * barcSq * sigma^2.
   * */
  void setInlierCostThresholds(const double inth) {
    barcSq_ = inth * Vector::Ones(nfg_.size());
  }

  /** Set the maximum weighted residual error for an inlier (one for each factor). For a factor in the form f(x) = 0.5 * || r(x) ||^2_Omega,
   * the inlier threshold is the largest value of f(x) for the corresponding measurement to be considered an inlier.
   * In other words, an inlier at x is such that 0.5 * || r(x) ||^2_Omega <= barcSq.
   * */
  void setInlierCostThresholds(const Vector& inthVec) {
    barcSq_ = inthVec;
  }

  /** Set the maximum weighted residual error threshold by specifying the probability
   * alpha that the inlier residuals are smaller than that threshold
   * */
  void setInlierCostThresholdsAtProbability(const double alpha) {
    barcSq_  = Vector::Ones(nfg_.size()); // initialize
    for (size_t k = 0; k < nfg_.size(); k++) {
      if (nfg_[k]) {
        barcSq_[k] = 0.5 * Chi2inv(alpha, nfg_[k]->dim()); // 0.5 derives from the error definition in gtsam
      }
    }
  }

  /** Set weights for each factor. This is typically not needed, but
   * provides an extra interface for the user to initialize the weightst
   * */
  void setWeights(const Vector w) {
    if (size_t(w.size()) != nfg_.size()) {
      throw std::runtime_error(
          "QGMOptimizer::setWeights: the number of specified weights"
          " does not match the size of the factor graph.");
    }
    weights_ = w;
  }

  /// Access a copy of the internal factor graph.
  const NonlinearFactorGraph& getFactors() const { return nfg_; }

  /// Access a copy of the internal values.
  const Values& getState() const { return state_; }

  /// Access a copy of the parameters.
  const QGMParameters& getParams() const { return params_;}

  /// Access a copy of the QGM weights.
  const Vector& getWeights() const { return weights_;}

  /// Get the inlier threshold.
  const Vector& getInlierCostThresholds() const {return barcSq_;}

  /// Equals.
  bool equals(const QGMOptimizer& other, double tol = 1e-9) const {
    return nfg_.equals(other.getFactors())
        && equal(weights_, other.getWeights())
        && params_.equals(other.getParams())
        && equal(barcSq_, other.getInlierCostThresholds());
  }

  Vector initializeWeightsFromKnownInliersAndOutliers() const{
    Vector weights = Vector::Ones(nfg_.size());
    for (size_t i = 0; i < params_.knownOutliers.size(); i++){
      weights[ params_.knownOutliers[i] ] = 0.0; // known to be outliers
    }
    return weights;
  }

  /// Compute optimal solution using QGM-based pose graph optimization.
  Values optimize() {
    NonlinearFactorGraph graph_initial = this->makeWeightedGraph(weights_);
    BaseOptimizer baseOptimizer(
        graph_initial, state_, params_.baseOptimizerParams);
    Values result = baseOptimizer.optimize();
    double prev_cost = graph_initial.error(result);
    double cost = 0.0;  // this will be updated in the main loop

    int nrUnknownInOrOut = nfg_.size() - ( params_.knownInliers.size() + params_.knownOutliers.size() );
    // ^^ number of measurements that are not known to be inliers or outliers 
    if (nrUnknownInOrOut == 0) { // no need to even call QGM in this case
      if (nrUnknownInOrOut==0 && params_.verbosity >= QGMParameters::Verbosity::SUMMARY) {
        std::cout << "QGM Optimizer stopped because all measurements are already known to be inliers or outliers"
                  << std::endl;
      }
      if (params_.verbosity >= QGMParameters::Verbosity::VALUES) {
        result.print("result\n");
      }
      return result;
    }

    size_t iter;
    for (iter = 0; iter < params_.maxIterations; iter++) {

      // display info
      if (params_.verbosity >= QGMParameters::Verbosity::WEIGHTS) {
        std::cout << "weights: " << weights_ << std::endl;
      }
      if (params_.verbosity >= QGMParameters::Verbosity::VALUES) {
        result.print("result\n");
      }
      // weights update
      weights_ = calculateWeights(result);

      // variable/values update
      NonlinearFactorGraph graph_iter = this->makeWeightedGraph(weights_);
      BaseOptimizer baseOptimizer_iter(
          graph_iter, state_, params_.baseOptimizerParams);
      result = baseOptimizer_iter.optimize();

      // stopping condition
      cost = graph_iter.error(result);
      if (checkCostConvergence(cost, prev_cost)) {
        break;
      }

      // get ready for next iteration
      prev_cost = cost;

      // display info
      if (params_.verbosity >= QGMParameters::Verbosity::VALUES) {
        std::cout << "previous cost: " << prev_cost << std::endl;
        std::cout << "current cost: " << cost << std::endl;
      }
    }
    // display info
    if (params_.verbosity >= QGMParameters::Verbosity::SUMMARY) {
      std::cout << "final iterations: " << iter << std::endl;
      std::cout << "previous cost: " << prev_cost << std::endl;
      std::cout << "current cost: " << cost << std::endl;
    }
    if (params_.verbosity >= QGMParameters::Verbosity::WEIGHTS) {
      std::cout << "final weights: " << weights_ << std::endl;
    }
    return result;
  }

  /// Check convergence of relative cost differences.
  bool checkCostConvergence(const double cost, const double prev_cost) const {
    bool costConverged = std::fabs(cost - prev_cost) / std::max(prev_cost, 1e-7)
        < params_.relativeCostTol;
    if (costConverged && params_.verbosity >= QGMParameters::Verbosity::SUMMARY){
      std::cout << "checkCostConvergence = true (prev. cost = " << prev_cost
                << ", curr. cost = " << cost << ")" << std::endl;
    }
    return costConverged;
  }

  /// Create a graph where each factor is weighted by the QGM weights.
  NonlinearFactorGraph makeWeightedGraph(const Vector& weights) const {
    // make sure all noiseModels are Gaussian or convert to Gaussian
    NonlinearFactorGraph newGraph;
    newGraph.resize(nfg_.size());
    for (size_t i = 0; i < nfg_.size(); i++) {
      if (nfg_[i]) {
        auto factor = nfg_.at<NoiseModelFactor>(i);
        auto noiseModel = std::dynamic_pointer_cast<noiseModel::Gaussian>(
            factor->noiseModel());
        if (noiseModel) {
          Matrix newInfo = weights[i] * noiseModel->information();
          auto newNoiseModel = noiseModel::Gaussian::Information(newInfo);
          newGraph[i] = factor->cloneWithNewNoiseModel(newNoiseModel);
        } else {
          throw std::runtime_error(
              "QGMOptimizer::makeWeightedGraph: unexpected non-Gaussian noise model.");
        }
      }
    }
    return newGraph;
  }

  /// Calculate QGM weights.
  Vector calculateWeights(const Values& currentEstimate) {
    Vector weights = initializeWeightsFromKnownInliersAndOutliers();

    // do not update the weights that the user has decided are known inliers
    std::vector<size_t> allWeights;
    for (size_t k = 0; k < nfg_.size(); k++) {
      allWeights.push_back(k);
    }
    std::vector<size_t> knownWeights;
    std::set_union(params_.knownInliers.begin(), params_.knownInliers.end(),
                   params_.knownOutliers.begin(), params_.knownOutliers.end(),
                   std::inserter(knownWeights, knownWeights.begin()));

    std::vector<size_t> unknownWeights;
    std::set_difference(allWeights.begin(), allWeights.end(),
                        knownWeights.begin(), knownWeights.end(),
                        std::inserter(unknownWeights, unknownWeights.begin()));

    // update weights of known inlier/outlier measurements
    for (size_t k : unknownWeights) {
      if (nfg_[k]) {
        double r2_k = nfg_[k]->error(currentEstimate); // squared residual
        weights[k] = 1 / std::pow((r2_k + barcSq_[k]), 2);
      }
    }
    return weights;
  }
};

}
