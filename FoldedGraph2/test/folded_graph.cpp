#include "gtest/gtest.h"
#include "folded_graph2.h"

#include <map>
#include <set>
#include <memory>
#include <random>

namespace crag {

namespace {

TEST(FoldedGraph2, Trivial) {
  FoldedGraph2 g;
  EXPECT_EQ(0, g.vertex(g.root()).endpoint(0));
  EXPECT_EQ(0, g.vertex(g.root()).endpoint(1));
  EXPECT_EQ(0, g.vertex(g.root()).endpoint(2));
  EXPECT_EQ(0, g.vertex(g.root()).endpoint(3));
}

FoldedGraph2::Word AlphasToLabels(std::vector<int> alphas) {
  FoldedGraph2::Word result;
  result.reserve(alphas.size());
  for (auto alpha : alphas) {
    switch(alpha) {
      case 1:
        result.push_back(0);
        break;
      case -1:
        result.push_back(1);
        break;
      case 2:
        result.push_back(2);
        break;
      case -2:
        result.push_back(3);
        break;
      default:
        assert(false);
    }
  }
  return result;
}

TEST(FoldedGraph2, ReadWord) {
  FoldedGraph2 g;
  EXPECT_EQ(std::make_tuple(1, 0), g.ReadWord(AlphasToLabels({1, 2})));

  EXPECT_EQ(3, g.PushWord(AlphasToLabels({1, 2})));
  EXPECT_EQ(std::make_tuple(3, 2), g.ReadWord(AlphasToLabels({1, 2})));
  EXPECT_EQ(std::make_tuple(2, 1), g.ReadWord(AlphasToLabels({1, 1})));

  EXPECT_EQ(std::make_tuple(1, 2), g.ReadInverse(AlphasToLabels({1, 2}), 3));
  EXPECT_EQ(std::make_tuple(2, 1), g.ReadInverse(AlphasToLabels({2, 2}), 3));
}

TEST(FoldedGraph2, PushWord) {
  FoldedGraph2 g;
  EXPECT_EQ(3, g.PushWord(AlphasToLabels({1, 2})));

  EXPECT_EQ(2, g.vertex(1).endpoint(0));
  EXPECT_EQ(0, g.vertex(1).endpoint(2));
  EXPECT_EQ(0, g.vertex(1).endpoint(1));
  EXPECT_EQ(0, g.vertex(1).endpoint(3));

  EXPECT_EQ(0, g.vertex(2).endpoint(0));
  EXPECT_EQ(3, g.vertex(2).endpoint(2));
  EXPECT_EQ(1, g.vertex(2).endpoint(1));
  EXPECT_EQ(0, g.vertex(2).endpoint(3));

  EXPECT_EQ(0, g.vertex(3).endpoint(0));
  EXPECT_EQ(0, g.vertex(3).endpoint(2));
  EXPECT_EQ(0, g.vertex(3).endpoint(1));
  EXPECT_EQ(2, g.vertex(3).endpoint(3));

  EXPECT_EQ(3, g.PushWord(AlphasToLabels({1, 2})));
  EXPECT_EQ(4, g.PushWord(AlphasToLabels({1, 1})));

  EXPECT_EQ(2, g.vertex(1).endpoint(0));
  EXPECT_EQ(0, g.vertex(1).endpoint(2));
  EXPECT_EQ(0, g.vertex(1).endpoint(1));
  EXPECT_EQ(0, g.vertex(1).endpoint(3));

  EXPECT_EQ(4, g.vertex(2).endpoint(0));
  EXPECT_EQ(3, g.vertex(2).endpoint(2));
  EXPECT_EQ(1, g.vertex(2).endpoint(1));
  EXPECT_EQ(0, g.vertex(2).endpoint(3));

  EXPECT_EQ(0, g.vertex(3).endpoint(0));
  EXPECT_EQ(0, g.vertex(3).endpoint(2));
  EXPECT_EQ(0, g.vertex(3).endpoint(1));
  EXPECT_EQ(2, g.vertex(3).endpoint(3));

  EXPECT_EQ(0, g.vertex(4).endpoint(0));
  EXPECT_EQ(0, g.vertex(4).endpoint(2));
  EXPECT_EQ(2, g.vertex(4).endpoint(1));
  EXPECT_EQ(0, g.vertex(4).endpoint(3));

  EXPECT_EQ(5, g.PushWord(AlphasToLabels({1, 2}), 2));

  EXPECT_EQ(2, g.vertex(1).endpoint(0));
  EXPECT_EQ(0, g.vertex(1).endpoint(2));
  EXPECT_EQ(0, g.vertex(1).endpoint(1));
  EXPECT_EQ(0, g.vertex(1).endpoint(3));

  EXPECT_EQ(4, g.vertex(2).endpoint(0));
  EXPECT_EQ(3, g.vertex(2).endpoint(2));
  EXPECT_EQ(1, g.vertex(2).endpoint(1));
  EXPECT_EQ(0, g.vertex(2).endpoint(3));

  EXPECT_EQ(0, g.vertex(3).endpoint(0));
  EXPECT_EQ(0, g.vertex(3).endpoint(2));
  EXPECT_EQ(0, g.vertex(3).endpoint(1));
  EXPECT_EQ(2, g.vertex(3).endpoint(3));

  EXPECT_EQ(0, g.vertex(4).endpoint(0));
  EXPECT_EQ(5, g.vertex(4).endpoint(2));
  EXPECT_EQ(2, g.vertex(4).endpoint(1));
  EXPECT_EQ(0, g.vertex(4).endpoint(3));

  EXPECT_EQ(0, g.vertex(5).endpoint(0));
  EXPECT_EQ(0, g.vertex(5).endpoint(2));
  EXPECT_EQ(0, g.vertex(5).endpoint(1));
  EXPECT_EQ(4, g.vertex(5).endpoint(3));

}

TEST(FoldedGraph2, JoinVertices1) {
  FoldedGraph2 g;
  EXPECT_EQ(2, g.PushWord({0}));
  EXPECT_EQ(3, g.PushWord({2}));

  g.JoinVertices(2, 3);
  EXPECT_TRUE(g.Equal(2, 3));

  EXPECT_EQ(0, g.vertex(2).endpoint(0));
  EXPECT_EQ(0, g.vertex(2).endpoint(2));
  EXPECT_EQ(1, g.vertex(2).endpoint(1));
  EXPECT_EQ(1, g.vertex(2).endpoint(3));

  EXPECT_EQ(std::make_tuple(2, 1), g.ReadWord({0}));
  EXPECT_EQ(std::make_tuple(2, 1), g.ReadWord({2}));
}

TEST(FoldedGraph2, JoinVertices2) {
  FoldedGraph2 g;
  EXPECT_EQ(2, g.PushWord({0}));

  g.JoinVertices(1, 2);
  EXPECT_TRUE(g.Equal(2, 1));

  EXPECT_EQ(std::make_tuple(1, 1), g.ReadWord({0}));
  EXPECT_EQ(std::make_tuple(1, 0), g.ReadWord({2}));
  EXPECT_EQ(std::make_tuple(1, 1), g.ReadWord({0, 2}));
  EXPECT_EQ(std::make_tuple(1, 2), g.ReadWord({0, 0, 2, 0}));
}

TEST(FoldedGraph2, PushCycle1) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({0, 2}));

  EXPECT_EQ(std::make_tuple(1, 2), g.ReadWord({0, 2}));
  EXPECT_EQ(std::make_tuple(1, 4), g.ReadWord({0, 2, 0, 2}));
}

TEST(FoldedGraph2, PushCycle2) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({0, 2, 0}));
  EXPECT_EQ(true, g.PushCycle({0, 0}));

  EXPECT_EQ(std::make_tuple(1, 2), g.ReadWord({0, 0}));
  EXPECT_EQ(std::make_tuple(1, 3), g.ReadWord({0, 2, 0}));
  EXPECT_EQ(std::make_tuple(1, 4), g.ReadWord({0, 2, 2, 0}));
  EXPECT_EQ(std::make_tuple(1, 5), g.ReadWord({0, 2, 2, 2, 0}));
}

TEST(FoldedGraph2, PushCycle3) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({0, 2, 0}));
  EXPECT_EQ(std::make_tuple(1, 3), g.ReadWord({0, 2, 0}));

  EXPECT_EQ(true, g.PushCycle({0, 2, 0, 2, 0}));
  EXPECT_EQ(std::make_tuple(1, 2), g.ReadWord({2, 0}));
}

TEST(FoldedGraph2, PushCycle4) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({0, 2, 0}));
  EXPECT_EQ(false, g.PushCycle({0, 2, 0}));
  EXPECT_EQ(std::make_tuple(1, 3), g.ReadWord({0, 2, 0}));

  EXPECT_EQ(true, g.PushCycle({0, 2}));
  EXPECT_TRUE(g.Equal(1, 2));
  EXPECT_TRUE(g.Equal(2, 3));
}

TEST(FoldedGraph2, PushCycle5) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({1, 0, 0}));
}

TEST(FoldedGraph2, PushCycle6) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({0, 0, 0}));
  EXPECT_EQ(true, g.PushCycle({1, 1}));
}


//struct AWord {
//
//  const unsigned int kMaxLength;
//  FoldedGraph2::Word word_;
//  unsigned long state_;
//  unsigned long current_mask_;
//  bool in_initial_ = true;
//
//  AWord(unsigned int kMaxLength)
//    : kMaxLength(kMaxLength)
//    , word_({0, 0})
//    , state_(0)
//    , current_mask_(15)
//  { }
//
//  void next() {
//    ++state_;
//    in_initial_ = false;
//    if ((state_ & current_mask_) == 0) {
//      state_ = 0;
//      current_mask_ <<= 1;
//      current_mask_ ^= 1;
//      word_.emplace_back();
//      if (word_.size() > kMaxLength) {
//        word_ = {0, 0};
//        current_mask_ = 15;
//        in_initial_ = true;
//      }
//    }
//
//    for (size_t i = 0; i < word_.size() * 2; i += 2) {
//      unsigned char symbol = (state_ & (3 << i)) >> i;
//      if (symbol % 2 == 1) {
//        word_[i / 2] = -symbol / 2 - 1;
//      } else {
//        word_[i / 2] = symbol / 2 + 1;
//      }
//    }
//  }
//};

TEST(FoldedGraph2, Harvest1) {
  FoldedGraph2 g;
  g.PushCycle(AlphasToLabels({1, 1, 2, 1}));
  g.PushCycle(AlphasToLabels({1, 2, 2, -1}));

  auto words = g.Harvest(10, 1, 4);

  std::vector<FoldedGraph2::Word> correct = {
    { 0, 0, 2 },
    { 0, 0, 2, 0, 0, 0, 2 },
    { 0, 0, 2, 0, 0, 2, 2, 0, 2 },
    { 0, 0, 2, 0, 0, 2, 2, 1, 1 },
    { 0, 0, 2, 0, 0, 3, 3, 0, 2 },
    { 0, 0, 2, 0, 0, 3, 3, 1, 1 },
    { 0, 2, 2, 0, 2 },
    { 0, 2, 2, 0, 2, 0, 0, 0, 2 },
    { 0, 2, 2, 1, 1 },
    { 0, 2, 2, 1, 1, 3, 1, 1, 1 },
    { 0, 2, 2, 2, 2, 0, 2 },
    { 0, 2, 2, 2, 2, 1, 1 },
    { 0, 2, 2, 2, 2, 2, 2, 0, 2 },
    { 0, 2, 2, 2, 2, 2, 2, 1, 1 },
    { 0, 3, 3, 0, 2 },
    { 0, 3, 3, 0, 2, 0, 0, 0, 2 },
    { 0, 3, 3, 1, 1 },
    { 0, 3, 3, 1, 1, 3, 1, 1, 1 },
    { 0, 3, 3, 3, 3, 0, 2 },
    { 0, 3, 3, 3, 3, 1, 1 },
    { 0, 3, 3, 3, 3, 3, 3, 0, 2 },
    { 0, 3, 3, 3, 3, 3, 3, 1, 1 },
    { 1 },
    { 1, 3, 1, 1, 1 },
    { 1, 3, 1, 1, 1, 3, 1, 1, 1 },
    { 1, 3, 1, 2, 2, 0, 2 },
    { 1, 3, 1, 2, 2, 1, 1 },
    { 1, 3, 1, 2, 2, 2, 2, 0, 2 },
    { 1, 3, 1, 2, 2, 2, 2, 1, 1 },
    { 1, 3, 1, 3, 3, 0, 2 },
    { 1, 3, 1, 3, 3, 1, 1 },
    { 1, 3, 1, 3, 3, 3, 3, 0, 2 },
    { 1, 3, 1, 3, 3, 3, 3, 1, 1 }
  };

  EXPECT_EQ(correct, words);
}

//TEST(FoldedGraph2, PushCycleStress) {
//#ifdef NDEBUG
//  const unsigned int kMaxLength = 4;
//  const unsigned int kWordsNum = 3;
//#else
//  const unsigned int kMaxLength = 5;
//  const unsigned int kWordsNum = 2;
//#endif
//
//  std::vector<AWord> words;
//  words.reserve(kWordsNum);
//  for (size_t i = 0; i < kWordsNum; ++i) {
//    words.emplace_back(kMaxLength);
//  }
//
//  words.back().next();
//
//  while (!words.back().in_initial_) {
//    for (auto& word : words) {
//      word.next();
//      if (!word.in_initial_) {
//        break;
//      }
//    }
//
//    FoldedGraph2 g;
//
//    for (auto& word : words) {
//      g.PushCycle(word.word_);
//      ASSERT_EQ(std::make_tuple(1, word.word_.size()), g.ReadWord(word.word_));
//    }
//  }
//}

TEST(FoldedGraph2, PushCycleStressRandom) {
#ifdef NDEBUG
  static const unsigned int kRepeat = 100000;
#else
  static const unsigned int kRepeat = 1000000;
#endif
  static const unsigned int kWords = 3;
  std::mt19937_64 engine;
  std::uniform_int_distribution<> random_letter(0, 3);
  std::uniform_int_distribution<size_t> random_length(2, 10);

  for (auto i = 0u; i < kRepeat; ++i) {
    FoldedGraph2 g;
    for (auto j = 0u; j < kWords; ++j) {
      FoldedGraph2::Word w;
      size_t length = random_length(engine);
      w.reserve(length);
      while(w.size() < length) {
        w.push_back(random_letter(engine));
      }

      g.PushCycle(w);
      ASSERT_EQ(std::make_tuple(1, length), g.ReadWord(w));
    }
  }
}



} //namespace

} //namespace crag
