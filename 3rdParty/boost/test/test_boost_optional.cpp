//
// Created by dpantele on 4/23/15.
//

#include "gtest/gtest.h"
#include "boost/optional.hpp"

namespace {

TEST(BoostOptional, Test) {
  boost::optional<int> test;
  EXPECT_TRUE(!test);
  test = 10;
  EXPECT_EQ(10, *test);
}

}