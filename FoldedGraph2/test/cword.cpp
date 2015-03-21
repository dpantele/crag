#include "gtest/gtest.h"
#include "compressed_word.h"

#include <chrono>
#include <map>
#include <memory>

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
    ASSERT_EQ(GetLabel(i), w.GetFront()) << "Iteration " << i;
    ASSERT_EQ(GetLabel(0), w.GetBack())  << "Iteration " << i;
  }
  auto front_l = GetLabel(CWord::kMaxLength - 1);
  for (auto i = 0u; i < CWord::kMaxLength; ++i) {
    ASSERT_EQ(GetLabel(i), w.GetBack());
    ASSERT_EQ(front_l, w.GetFront());
    w.PopBack();
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
  EXPECT_EQ(CWord({0, 3, 1}), a.Inverse());
}

TEST(CWord, Inverse2) {
  CWord a = {0, 2, 1, 3};
  a.Invert();
  EXPECT_EQ(CWord({2, 0, 3, 1}), a);
}

namespace naive_find_subword {

std::tuple<CWord::size_type, CWord::size_type, CWord::size_type> LongestCommonSubwordCyclic(CWord u, CWord v) {
  //this is actually not naive, but turns out that in general this is slower than naive approach
  std::vector<CWord::size_type> previous_longest_suffix(v.size() * 2 + 1);
  std::vector<CWord::size_type> current_longest_suffix(v.size() * 2 + 1);

  auto cur_u = u;
  auto cur_v = v;
  CWord::size_type max_common_length = 0;
  CWord::size_type u_common_begin = u.size();
  CWord::size_type v_common_begin = v.size();
  for(CWord::size_type i = 1u; i <= u.size() * 2; ++i) {
    std::swap(previous_longest_suffix, current_longest_suffix);
    for(CWord::size_type j = 1u; j <= v.size() * 2; ++j) {
      if (cur_u.GetFront() == cur_v.GetFront()) {
        current_longest_suffix[j] = previous_longest_suffix[j - 1] + 1;
        if (current_longest_suffix[j] > u.size()) {
          current_longest_suffix[j] = u.size();
        }
        if (current_longest_suffix[j] > v.size()) {
          current_longest_suffix[j] = v.size();
        }
        if (current_longest_suffix[j] >= max_common_length) {
          CWord::size_type cur_u_begin = (i - current_longest_suffix[j]) % u.size();
          CWord::size_type cur_v_begin = (j - current_longest_suffix[j]) % v.size();
          if (current_longest_suffix[j] > max_common_length || cur_u_begin < u_common_begin || (cur_u_begin == u_common_begin && cur_v_begin < v_common_begin)) {
            max_common_length = current_longest_suffix[j];
            u_common_begin = cur_u_begin;
            v_common_begin = cur_v_begin;
          }
        }
      } else {
        current_longest_suffix[j] = 0;
      }
      cur_v.PopFront();
      if (cur_v.Empty()) {
        cur_v = v;
      }
    }
    cur_u.PopFront();
    if (cur_u.Empty()) {
      cur_u = u;
    }
  }

  assert(cur_u == u && cur_v == v);
  if (max_common_length > u.size()) {
    max_common_length = u.size();
  }
  if (max_common_length > v.size()) {
    max_common_length = v.size();
  }

  return std::make_tuple(u_common_begin, v_common_begin, max_common_length);
}

}

TEST(LongestCommonSubwordCyclic, Example1) {
  CWord u("xy");
  CWord v("xY");
  EXPECT_EQ(std::make_tuple(0, 0, 1), LongestCommonSubwordCyclic(u, v));
  EXPECT_EQ(std::make_tuple(0, 0, 1), naive_find_subword::LongestCommonSubwordCyclic(u, v));
}

TEST(LongestCommonSubwordCyclic, Example2) {
  CWord u("xx");
  CWord v("xY");
  EXPECT_EQ(std::make_tuple(0, 0, 1), LongestCommonSubwordCyclic(u, v));
  EXPECT_EQ(std::make_tuple(0, 0, 1), naive_find_subword::LongestCommonSubwordCyclic(u, v));
}

TEST(LongestCommonSubwordCyclic, Example3) {
  CWord u("xYx");
  CWord v("xx");
  EXPECT_EQ(std::make_tuple(2, 0, 2), LongestCommonSubwordCyclic(u, v));
  EXPECT_EQ(std::make_tuple(2, 0, 2), naive_find_subword::LongestCommonSubwordCyclic(u, v));
}

TEST(LongestCommonSubwordCyclic, Example4) {
  CWord u("x");
  CWord v("xx");
  EXPECT_EQ(std::make_tuple(0, 0, 1), LongestCommonSubwordCyclic(u, v));
  EXPECT_EQ(std::make_tuple(0, 0, 1), naive_find_subword::LongestCommonSubwordCyclic(u, v));
}

TEST(LongestCommonSubwordCyclic, Example5) {
  CWord u("YXyyXyx");
  CWord v("yX");
  EXPECT_EQ(std::make_tuple(1, 1, 2), LongestCommonSubwordCyclic(u, v));
  EXPECT_EQ(std::make_tuple(1, 1, 2), naive_find_subword::LongestCommonSubwordCyclic(u, v));
}


TEST(LongestCommonSubwordCyclic, StressTest) {
  static const auto kDuration = std::chrono::seconds(10);
  std::mt19937_64 engine(17);
  RandomWord rw(30, 30);

  auto begin = std::chrono::steady_clock::now();

  std::chrono::high_resolution_clock::duration naive_op_duration{}, op_duration{};
  std::chrono::high_resolution_clock::time_point proc_begin;

  auto repeat = 0ull;
  while (std::chrono::steady_clock::now() - begin < kDuration) {
    ++repeat;
    auto u = rw(engine);
    auto v = rw(engine);

    proc_begin = std::chrono::high_resolution_clock::now();
    auto result = LongestCommonSubwordCyclic(u, v);
    op_duration += (std::chrono::high_resolution_clock::now() - proc_begin);

    proc_begin = std::chrono::high_resolution_clock::now();
    auto naive_result = naive_find_subword::LongestCommonSubwordCyclic(u, v);
    naive_op_duration += (std::chrono::high_resolution_clock::now() - proc_begin);
    ASSERT_EQ(naive_result, result) << u << ' ' << v;
  }
  std::cout << std::string(13, ' ') << repeat << " repeats" << std::endl;
  std::cout << std::string(13, ' ') 
    << std::chrono::duration_cast<std::chrono::milliseconds>(op_duration).count()
    << " vs "
    << std::chrono::duration_cast<std::chrono::milliseconds>(naive_op_duration).count()
    << std::endl;
  ASSERT_GT(repeat, 10000);
}



} //namespace

} //namespace crag
