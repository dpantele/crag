#include "gtest/gtest.h"
#include "compressed_word.h"

#include <map>
#include <set>
#include <memory>
#include <random>

namespace crag {

namespace {

TEST(CWord, Construct) {
  CWord w = {0, 1, 2, 3};
  EXPECT_EQ(0, w.size());
  w = {0, 2, 1, 3};
  EXPECT_EQ(0, w.GetFront());
  EXPECT_EQ(3, w.GetBack());
} 

TEST(CWord, PopFront) {
  CWord w = {0, 2, 1, 3};
  EXPECT_EQ(4, w.size());
  w.PopFront();
  EXPECT_EQ(CWord({2, 1, 3}), w);
  w.PopFront(2);
  EXPECT_EQ(CWord({3}), w);
} 

TEST(CWord, PopBack) {
  CWord w = {0, 2, 1, 3};
  EXPECT_EQ(4, w.size());
  w.PopBack();
  EXPECT_EQ(CWord({0, 2, 1}), w);
  EXPECT_NE(CWord({0, 2, 0}), w);
  w.PopBack(2);
  EXPECT_EQ(CWord({0}), w);
} 

TEST(CWord, Push) {
  CWord w = {0, 2, 1, 3};
  EXPECT_EQ(4, w.size());
  w.PushBack(2);
  EXPECT_EQ(CWord({0, 2, 1, 3, 2}), w);
  w.PushFront(1);
  EXPECT_EQ(CWord({2, 1, 3, 2 }), w);
} 

TEST(CWord, CyclicShift) {
  CWord w = {0, 2, 1, 3, 0, 2, 1, 3, 0, 3, 1, 2, 0, 3, 1, 2};
  w.CyclicLeftShift(2);
  EXPECT_EQ(CWord({1, 3, 0, 2, 1, 3, 0, 3, 1, 2, 0, 3, 1, 2, 0, 2}), w);
} 

TEST(CWord, Flip1) {
  CWord a = {1, 2, 3};
  a.Flip();
  EXPECT_EQ(CWord({3, 2, 1}), a);
}

TEST(CWord, Flip2) {
  CWord a = {1, 2, 3, 0};
  a.Flip();
  EXPECT_EQ(CWord({0, 3, 2, 1}), a);
}

TEST(CWord, Inverse1) {
  CWord a = {0, 2, 1};
  a.Invert();
  EXPECT_EQ(CWord({0, 3, 1}), a);
}

TEST(CWord, Inverse2) {
  CWord a = {0, 2, 1, 3};
  a.Invert();
  EXPECT_EQ(CWord({2, 0, 3, 1}), a);
}




} //namespace

} //namespace crag
