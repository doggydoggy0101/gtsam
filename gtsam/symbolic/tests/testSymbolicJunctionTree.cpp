/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation, 
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    testJunctionTree.cpp
 * @brief   Unit tests for Junction Tree
 * @author  Kai Ni
 * @author  Frank Dellaert
 */

#include <CppUnitLite/TestHarness.h>
#include <gtsam/base/TestableAssertions.h>

#include <gtsam/symbolic/SymbolicFactorGraph.h>
#include <gtsam/symbolic/SymbolicEliminationTree.h>
#include <gtsam/symbolic/SymbolicJunctionTree.h>

#include <boost/assign/list_of.hpp>
using namespace boost::assign;

#include "symbolicExampleGraphs.h"

using namespace gtsam;
using namespace std;

/* ************************************************************************* *
 * 1 - 0 - 2 - 3
 * 2 3
 *   0 1 : 2
 ****************************************************************************/
TEST( JunctionTree, constructor )
{
  Ordering order; order += 0, 1, 2, 3;

  SymbolicJunctionTree actual(SymbolicEliminationTree(simpleChain, order));

  vector<Key> frontal1 = list_of(2)(3);
  vector<Key> frontal2 = list_of(0)(1);
  vector<Key> sep1;
  vector<Key> sep2 = list_of(2);
  EXPECT(assert_equal(frontal1, actual.roots().front()->keys));
  //EXPECT(assert_equal(sep1,     actual.roots().front()->separator));
  LONGS_EQUAL(1,                (long)actual.roots().front()->factors.size());
  EXPECT(assert_equal(frontal2, actual.roots().front()->children.front()->keys));
  //EXPECT(assert_equal(sep2,     actual.roots().front()->children.front()->separator));
  LONGS_EQUAL(2,                (long)actual.roots().front()->children.front()->factors.size());
  EXPECT(assert_equal(*simpleChain[2],   *actual.roots().front()->factors[0]));
  EXPECT(assert_equal(*simpleChain[0],   *actual.roots().front()->children.front()->factors[0]));
  EXPECT(assert_equal(*simpleChain[1],   *actual.roots().front()->children.front()->factors[1]));
}

/* ************************************************************************* */
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */