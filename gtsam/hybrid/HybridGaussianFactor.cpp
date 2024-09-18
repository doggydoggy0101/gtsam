/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file   HybridGaussianFactor.cpp
 * @brief  A set of Gaussian factors indexed by a set of discrete keys.
 * @author Fan Jiang
 * @author Varun Agrawal
 * @author Frank Dellaert
 * @date   Mar 12, 2022
 */

#include <gtsam/base/utilities.h>
#include <gtsam/discrete/DecisionTree-inl.h>
#include <gtsam/discrete/DecisionTree.h>
#include <gtsam/hybrid/HybridGaussianFactor.h>
#include <gtsam/hybrid/HybridValues.h>
#include <gtsam/linear/GaussianFactor.h>
#include <gtsam/linear/GaussianFactorGraph.h>

namespace gtsam {

/**
 * @brief Helper function to augment the [A|b] matrices in the factor components
 * with the normalizer values.
 * This is done by storing the normalizer value in
 * the `b` vector as an additional row.
 *
 * @param factors DecisionTree of GaussianFactors and arbitrary scalars.
 * Gaussian factor in factors.
 * @return HybridGaussianFactor::Factors
 */
HybridGaussianFactor::Factors augment(
    const HybridGaussianFactor::FactorValuePairs &factors) {
  // Find the minimum value so we can "proselytize" to positive values.
  // Done because we can't have sqrt of negative numbers.
  HybridGaussianFactor::Factors gaussianFactors;
  AlgebraicDecisionTree<Key> valueTree;
  std::tie(gaussianFactors, valueTree) = unzip(factors);

  // Normalize
  double min_value = valueTree.min();
  AlgebraicDecisionTree<Key> values =
      valueTree.apply([&min_value](double n) { return n - min_value; });

  // Finally, update the [A|b] matrices.
  auto update = [&values](const Assignment<Key> &assignment,
                          const HybridGaussianFactor::sharedFactor &gf) {
    auto jf = std::dynamic_pointer_cast<JacobianFactor>(gf);
    if (!jf) return gf;
    // If the log_normalizer is 0, do nothing
    if (values(assignment) == 0.0) return gf;

    GaussianFactorGraph gfg;
    gfg.push_back(jf);

    Vector c(1);
    c << std::sqrt(values(assignment));
    auto constantFactor = std::make_shared<JacobianFactor>(c);

    gfg.push_back(constantFactor);
    return std::dynamic_pointer_cast<GaussianFactor>(
        std::make_shared<JacobianFactor>(gfg));
  };
  return gaussianFactors.apply(update);
}

/* *******************************************************************************/
HybridGaussianFactor::HybridGaussianFactor(const KeyVector &continuousKeys,
                                           const DiscreteKeys &discreteKeys,
                                           const FactorValuePairs &factors)
    : Base(continuousKeys, discreteKeys), factors_(augment(factors)) {}

/* *******************************************************************************/
bool HybridGaussianFactor::equals(const HybridFactor &lf, double tol) const {
  const This *e = dynamic_cast<const This *>(&lf);
  if (e == nullptr) return false;

  // This will return false if either factors_ is empty or e->factors_ is empty,
  // but not if both are empty or both are not empty:
  if (factors_.empty() ^ e->factors_.empty()) return false;

  // Check the base and the factors:
  return Base::equals(*e, tol) &&
         factors_.equals(e->factors_,
                         [tol](const sharedFactor &f1, const sharedFactor &f2) {
                           return f1->equals(*f2, tol);
                         });
}

/* *******************************************************************************/
void HybridGaussianFactor::print(const std::string &s,
                                 const KeyFormatter &formatter) const {
  std::cout << (s.empty() ? "" : s + "\n");
  std::cout << "HybridGaussianFactor" << std::endl;
  HybridFactor::print("", formatter);
  std::cout << "{\n";
  if (factors_.empty()) {
    std::cout << "  empty" << std::endl;
  } else {
    factors_.print(
        "", [&](Key k) { return formatter(k); },
        [&](const sharedFactor &gf) -> std::string {
          RedirectCout rd;
          std::cout << ":\n";
          if (gf) {
            gf->print("", formatter);
            return rd.str();
          } else {
            return "nullptr";
          }
        });
  }
  std::cout << "}" << std::endl;
}

/* *******************************************************************************/
HybridGaussianFactor::sharedFactor HybridGaussianFactor::operator()(
    const DiscreteValues &assignment) const {
  return factors_(assignment);
}

/* *******************************************************************************/
GaussianFactorGraphTree HybridGaussianFactor::add(
    const GaussianFactorGraphTree &sum) const {
  using Y = GaussianFactorGraph;
  auto add = [](const Y &graph1, const Y &graph2) {
    auto result = graph1;
    result.push_back(graph2);
    return result;
  };
  const auto tree = asGaussianFactorGraphTree();
  return sum.empty() ? tree : sum.apply(tree, add);
}

/* *******************************************************************************/
GaussianFactorGraphTree HybridGaussianFactor::asGaussianFactorGraphTree()
    const {
  auto wrap = [](const sharedFactor &gf) { return GaussianFactorGraph{gf}; };
  return {factors_, wrap};
}

/* *******************************************************************************/
AlgebraicDecisionTree<Key> HybridGaussianFactor::errorTree(
    const VectorValues &continuousValues) const {
  // functor to convert from sharedFactor to double error value.
  auto errorFunc = [&continuousValues](const sharedFactor &gf) {
    return gf->error(continuousValues);
  };
  DecisionTree<Key, double> error_tree(factors_, errorFunc);
  return error_tree;
}

/* *******************************************************************************/
double HybridGaussianFactor::error(const HybridValues &values) const {
  const sharedFactor gf = factors_(values.discrete());
  return gf->error(values.continuous());
}

}  // namespace gtsam