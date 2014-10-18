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
  EXPECT_EQ(CWord("x"), w);
} 

TEST(CWord, Push) {
  CWord w = {0, 2, 1, 3};
  EXPECT_EQ(4, w.size());
  w.PushBack(2);
  EXPECT_EQ(CWord({0, 2, 1, 3, 2}), w);
  w.PushFront(1);
  EXPECT_EQ(CWord({2, 1, 3, 2 }), w);
} 

inline unsigned int GetLabel(unsigned int i) {
  i %= (2 * CWord::kAlphabetSize);
  return i / CWord::kAlphabetSize + 2 * (i % CWord::kAlphabetSize);
}

TEST(CWord, PushFrontPopBack) {
  CWord w;
  for (auto i = 0u; i < CWord::kMaxLength; ++i) {
    w.PushFront(GetLabel(i));
    ASSERT_EQ(GetLabel(i), w.GetFront());
    ASSERT_EQ(GetLabel(0), w.GetBack());
  }
  auto front_l = GetLabel(CWord::kMaxLength - 1);
  for (auto i = 0u; i < CWord::kMaxLength; ++i) {
    ASSERT_EQ(GetLabel(i), w.GetBack());
    ASSERT_EQ(front_l, w.GetFront());
    w.PopBack();
  }
} 

TEST(CWord, LeftCyclicShift17) {
  CWord w;
  for (auto i = 0u; i < 55; ++i) {
    w.PushBack(GetLabel(i));
  }
  w.CyclicLeftShift(17);
  for (auto i = 0u; i < 55; ++i) {
    ASSERT_EQ(GetLabel((i + 17) % 55), w.GetFront());
    w.PopFront();
  }
} 

TEST(CWord, RightCyclicShift17) {
  CWord w;
  for (auto i = 0u; i < 55; ++i) {
    w.PushBack(GetLabel(i));
  }
  w.CyclicRightShift(17);
  for (auto i = 0u; i < 55; ++i) {
    ASSERT_EQ(GetLabel((i + 55 - 17) % 55), w.GetFront());
    w.PopFront();
  }
} 

TEST(CWord, Flip55) {
  CWord w;
  for (auto i = 0u; i < 55; ++i) {
    w.PushBack(GetLabel(i));
  }
  w.Flip();
  for (auto i = 0u; i < 55; ++i) {
    ASSERT_EQ(GetLabel(54 - i), w.GetFront());
    w.PopFront();
  }
} 

TEST(CWord, Flip30) {
  CWord w;
  for (auto i = 0u; i < 30; ++i) {
    w.PushBack(GetLabel(i));
  }
  w.Flip();
  for (auto i = 0u; i < 30; ++i) {
    ASSERT_EQ(GetLabel(29 - i), w.GetFront());
    w.PopFront();
  }
} 

TEST(CWord, Inverse55) {
  CWord w;
  for (auto i = 0u; i < 55; ++i) {
    w.PushBack(GetLabel(i));
  }
  w.Invert();
  for (auto i = 0u; i < 55; ++i) {
    ASSERT_EQ(GetLabel(54 - i) ^ 1, w.GetFront());
    w.PopFront();
  }
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

TEST(CWord, StringInput) {
  CWord a("xyxYXY");
  EXPECT_EQ(a.GetBack(), 3);
  a.PopBack();
  EXPECT_EQ(a.GetBack(), 1);
  a.PopBack();
  EXPECT_EQ(a.GetBack(), 3);
  a.PopBack();
  EXPECT_EQ(a.GetBack(), 0);
  a.PopBack();
  EXPECT_EQ(a.GetBack(), 2);
  a.PopBack();
  EXPECT_EQ(a.GetBack(), 0);
  a.PopBack();
}




} //namespace

} //namespace crag
