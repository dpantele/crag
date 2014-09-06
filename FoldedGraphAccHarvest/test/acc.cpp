#include "gtest/gtest.h"
#include "acc.h"

#include <map>
#include <set>
#include <memory>
#include <random>

namespace crag {

namespace {

TEST(WordOp, Flip1) {
  Word a = {1, 2, 3};
  Flip(a.begin(), a.end());
  EXPECT_EQ(Word({3, 2, 1}), a);
}

TEST(WordOp, Flip2) {
  Word a = {1, 2, 3, 4};
  Flip(a.begin(), a.end());
  EXPECT_EQ(Word({4, 3, 2, 1}), a);
}

TEST(WordOp, LeftShift1) {
  Word a = {1, 2, 3};
  LeftShift(a.begin(), a.end());
  EXPECT_EQ(Word({2, 3, 1}), a);
}

TEST(WordOp, MinPermuation) {
  Word a = {2, 0, 1, 0, 0};
  PermuteToMin(&a);
  EXPECT_EQ(Word({0, 0, 2, 0, 1}), a);
}

TEST(WordOp, Inverse1) {
  Word a = {0, 2, 1};
  Invert(a.begin(), a.end());
  EXPECT_EQ(Word({0, 3, 1}), a);
}

TEST(WordOp, Inverse2) {
  Word a = {0, 2, 1, 3};
  Invert(a.begin(), a.end());
  EXPECT_EQ(Word({2, 0, 3, 1}), a);
}

TEST(WordOp, CyclicReduce1) {
  Word a = {0, 2, 2, 1};
  EXPECT_EQ(Word({2, 2}), CyclicReduce(a));
}

TEST(WordOp, CyclicReduce2) {
  Word a = {0, 2, 3, 1};
  EXPECT_EQ(Word({}), CyclicReduce(a));
}

TEST(WordOp, CyclicReduce3) {
  Word a = {0, 2, 3, 0};
  EXPECT_EQ(Word({0, 2, 3, 0}), CyclicReduce(a));
}

TEST(WordOp, CyclicReduce4) {
  Word a = {0, 2, 0};
  EXPECT_EQ(Word({0, 2, 0}), CyclicReduce(a));
}

TEST(WordOp, CyclicReduce5) {
  Word a = {0, 2, 1};
  EXPECT_EQ(Word({2u}), CyclicReduce(a));
}

TEST(WordOp, ReduceAndNormalize1) {
  EXPECT_EQ(std::vector<Word>({{0}, {0, 2, 1}, {2}}), ReduceAndNormalize({{0}, {1}, {0, 2, 1}, {2, 2, 1, 0, 3}}));
}

TEST(WordOp, Conjugate1) {
  EXPECT_EQ(Word({0, 0, 2, 1}), Conjugate({0u, 2u}, {1u}));
}

TEST(WordOp, Conjugate2) {
  EXPECT_EQ(Word({0u}), Conjugate({0u}, {0u, 0u}));
}

} //namespace

} //namespace crag
