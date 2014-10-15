#include "gtest/gtest.h"
#include "acc.h"

#include <chrono>
#include <map>
#include <set>
#include <memory>
#include <random>

namespace crag {

namespace {

//TEST(WordOp, Flip1) {
//  Word a = {1, 2, 3};
//  Flip(a.begin(), a.end());
//  EXPECT_EQ(Word({3, 2, 1}), a);
//}
//
//TEST(WordOp, Flip2) {
//  Word a = {1, 2, 3, 4};
//  Flip(a.begin(), a.end());
//  EXPECT_EQ(Word({4, 3, 2, 1}), a);
//}
//
//TEST(WordOp, LeftShift1) {
//  Word a = {1, 2, 3};
//  LeftShift(a.begin(), a.end());
//  EXPECT_EQ(Word({2, 3, 1}), a);
//}

TEST(WordOp, MinPermuation) {
  Word a = {2, 0, 1, 0, 0};
  PermuteToMin(&a);
  EXPECT_EQ(Word({0, 0, 2, 0, 1}), a);
}

//TEST(WordOp, Inverse1) {
//  Word a = {0, 2, 1};
//  Invert(a.begin(), a.end());
//  EXPECT_EQ(Word({0, 3, 1}), a);
//}
//
//TEST(WordOp, Inverse2) {
//  Word a = {0, 2, 1, 3};
//  Invert(a.begin(), a.end());
//  EXPECT_EQ(Word({2, 0, 3, 1}), a);
//}

TEST(WordOp, CyclicReduce1) {
  Word a = {0, 2, 2, 1};
  EXPECT_EQ(Word({2, 2}), CyclicReduce(a));
}

TEST(WordOp, CyclicReduce2) {
  Word a = {0, 2, 3, 1};
  EXPECT_EQ(Word(), CyclicReduce(a));
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

TEST(WordOp, ReduceAndMinCycle1) {
  EXPECT_EQ(std::vector<Word>({{0}, {2}}), ReduceAndMinCycle({{0}, {1}, {0, 2, 1}, {2, 2, 1, 0, 3}}));
}

TEST(WordOp, Conjugate1) {
  EXPECT_EQ(Word({0, 0, 2, 1}), Conjugate({0u, 2u}, {1u}));
}

TEST(WordOp, Conjugate2) {
  EXPECT_EQ(Word({0u}), Conjugate({0u}, {0u, 0u}));
}

TEST(GenAllWords, L1) {
  EXPECT_EQ(std::set<Word>({Word(), Word({0u}), Word({1u}), Word({2u}), Word({3u})}), GenAllWords(1));
}

TEST(GenAllWords, L2) {
  EXPECT_EQ(std::set<Word>({
    {},
    {0u},
    {1u},
    {2u},
    {3u},
    {0u, 0u},
    {0u, 2u},
    {0u, 3u},
    {1u, 1u},
    {1u, 2u},
    {1u, 3u},
    {2u, 0u},
    {2u, 1u},
    {2u, 2u},
    {3u, 0u},
    {3u, 1u},
    {3u, 3u},
  }), GenAllWords(2));
}

namespace naive_get_canonical_pair {

Word Inverse(Word u) {
  u.Invert();
  return u;
}


std::pair<Word, Word> GetCanonicalPair(Word u, Word v) {
#define MORPHISM(X, Y) {Word(X), Inverse(Word(X)), Word(Y), Inverse(Word(Y))}
  static const Mapping morphisms[] = {
    MORPHISM("x", "y"),
    MORPHISM("x", "Y"),
    MORPHISM("X", "y"),
    MORPHISM("X", "Y"),
    MORPHISM("y", "x"),
    MORPHISM("y", "X"),
    MORPHISM("Y", "x"),
    MORPHISM("Y", "X"),
  };
#undef MORPHISM

  u = CyclicReduce(u);
  v = CyclicReduce(v);

  auto min_length_pairs = MinimizeTotalLength(u, v);
  assert(!min_length_pairs.empty());
  std::tie(u, v) = *min_length_pairs.begin();
  PermuteToMinWithInverse(&u);
  PermuteToMinWithInverse(&v);
  for (auto&& uv : min_length_pairs) {
    for(auto&& automorph : morphisms) {
      auto up = uv.first;
      auto vp = uv.second;
      ReduceMapAndMinCycle(automorph, &up);
      ReduceMapAndMinCycle(automorph, &vp);

      if (vp < up) {
        std::swap(up, vp);
      }

      if (up < u || (up == u && vp < v)) {
        u = up;
        v = vp;
      }
    }
  }
  return std::make_pair(u, v);
}


} //naive_get_canonical_pair

Word Inverse(Word u) {
  u.Invert();
  return u;
}

typedef std::map<std::pair<Word, Word>, std::vector<size_t>> Orbit;
typedef std::vector<Orbit::const_iterator> NewOrbitElements;

NewOrbitElements ProduceAutomorhicOrbit(const NewOrbitElements& just_added_elems, Orbit* full_orbit) {
#define MORPHISM(X, Y) {Word(X), Inverse(Word(X)), Word(Y), Inverse(Word(Y))}
  static const Mapping morphisms[] = {
    MORPHISM("x", "y"),   //0
    MORPHISM("x", "Y"),   //1
    MORPHISM("X", "y"),   //2
    MORPHISM("X", "Y"),   //3
    MORPHISM("y", "x"),   //4
    MORPHISM("y", "X"),   //5
    MORPHISM("Y", "x"),   //6
    MORPHISM("Y", "X"),   //7
    MORPHISM("yx", "y"),  //8
    MORPHISM("Yx", "y"),  //9
    MORPHISM("xy", "y"),  //10
    MORPHISM("xY", "y"),  //11
    MORPHISM("yxY", "y"), //12
    MORPHISM("Yxy", "y"), //13
    MORPHISM("x", "yx"),  //14
    MORPHISM("x", "yX"),  //15
    MORPHISM("x", "xy"),  //16
    MORPHISM("x", "Xy"),  //17
    MORPHISM("x", "Xyx"), //18
    MORPHISM("x", "xyX"), //19
  };
#undef MORPHISM
  NewOrbitElements new_orbit;
  for (auto&& just_added_elem : just_added_elems) {
    auto morphism_id = 0u;
    for (auto&& morhpism : morphisms) {
      try {
        auto images = std::make_pair(
          CyclicReduce(Map(just_added_elem->first.first, morhpism)),
          CyclicReduce(Map(just_added_elem->first.second, morhpism))
        );

        auto new_elem = full_orbit->emplace(images, just_added_elem->second);
        if (new_elem.second) {
          new_elem.first->second.push_back(morphism_id);
          new_orbit.push_back(new_elem.first);
        } 
      } catch(std::length_error&) {
        //it's ok
      }
      ++morphism_id;
    }
  }
  return new_orbit;
}


#define PAIR(X, Y) std::make_pair(Word(X), Word(Y))
static const std::pair<Word, Word> random_pairs[] = {
  PAIR("xyxYXY", "xxxYYYY"),
  PAIR("xyxyXYY", "xyxyxYYXy"),
  PAIR("xxYxyXXY", "xxyXXyxYY"),
  PAIR("xxyXyXY", "xxyXYxyXy"),
  PAIR("xyXYYxyy", "xxyXYXy"),
  PAIR("xyxYxYxyXY", "xyxYxyXYY"),
  PAIR("xxYXyxxYY", "xxYYxYYxy"),
  PAIR("xyyxY", "xxxxyXXYxy"),
  PAIR("xyxYXYYXy", "xyxYYYYXY"),
  PAIR("xxyXYXy", "xyyxYYXy"),
  PAIR("xxYYXyyXy", "xxYxyXXY"),
  PAIR("xyXyXYY", "xxxYYXy"),
  PAIR("xxxyyyxYY", "xxxxyyxY"),
  PAIR("xyxyxYxyXY", "xyxYxyyXY"),
  PAIR("xxyXYXyXY", "xxxxYxyxxy"),
  PAIR("xyxyyyyXy", "xxxyyyyXXy"),
  PAIR("xxYXYxyXy", "xxxyxxYxY"),
  PAIR("xyyyXYxYY", "xxxxyyyXYY"),
  PAIR("xxxxyXXYXY", "xxyyXYY"),
  PAIR("xxxyyXY", "xxyXYYXXy"),
  PAIR("xyXyxYxYXy", "xxyXYxYXy"),
};
#undef PAIR

TEST(GetCanonicalPair, Naive1) {
  auto initial = std::make_pair(Word("xyXXyxY"), Word("xxyXXYxy"));
  auto expected = std::make_pair(Word("xxyXy"), Word("xyXYYxyyyy"));
  EXPECT_EQ(expected, naive_get_canonical_pair::GetCanonicalPair(initial.first, initial.second)) << ::testing::PrintToString(MinimizeTotalLength(initial.first, initial.second));
}

TEST(GetCanonicalPair, Naive2) {
  auto initial = std::make_pair(Word("xyXXXyxY"), Word("xyXXYxy"));
  auto expected = std::make_pair(Word("xxyXy"), Word("xyXYYxyyyy"));
  EXPECT_EQ(expected, naive_get_canonical_pair::GetCanonicalPair(initial.first, initial.second)) << ::testing::PrintToString(MinimizeTotalLength(initial.first, initial.second));
}

TEST(GetCanonicalPair, Naive3) {
  auto initial = std::make_pair(Word("xYXXXYxy"), Word("xYXXyxY"));
  auto expected = std::make_pair(Word("xxyXy"), Word("xyXYYxyyyy"));
  EXPECT_EQ(expected, naive_get_canonical_pair::GetCanonicalPair(initial.first, initial.second)) << ::testing::PrintToString(MinimizeTotalLength(initial.first, initial.second));
}



size_t count_size_16_and_less(const Orbit& pairs) {
  size_t count = 0u;
  for (auto&& pair : pairs) {
    if (pair.first.first.size() <= 16 && pair.first.second.size() <= 16) {
      ++count;
    }
  }
  return count;
}

struct Stopwatch {
  using clock = std::chrono::high_resolution_clock;
  clock::duration duration_{0};
  clock::duration last_duration_{0};
  uint64_t iterations_count_ = 0u;

  typedef std::chrono::microseconds ResultUnits;

  struct Iteration {
    Stopwatch* stopwatch_;
    clock::duration duration_{0};
    clock::time_point last_click_;
    uint64_t clicks_count_{0};

    Iteration(Stopwatch* stopwatch)
      : stopwatch_(stopwatch)
    { }

    inline bool IsRunning() const {
      return clicks_count_ % 2 == 1;
    }

    void Click() {
      if (IsRunning()) {
        duration_ += (clock::now() - last_click_);
      } else {
        last_click_ = clock::now();
      }

      ++clicks_count_;
    }

    ~Iteration() {
      if (IsRunning()) {
        Click();
      }
      stopwatch_->Report(this);
    }
  };

  void Report(Iteration* iteration) {
    ++iterations_count_;
    duration_ += iteration->duration_;
    last_duration_ = iteration->duration_;
  }

  Iteration NewIter() {
    return Iteration(this);
  }

  long long last() {
    return std::chrono::duration_cast<ResultUnits>(last_duration_).count();
  }

  long long average() {
    return iterations_count_ ? std::chrono::duration_cast<ResultUnits>(duration_ / iterations_count_).count() : 0;
  }
};

TEST(GetCanonicalPair, Naive) {
  for (auto&& initial_pair : random_pairs) {
    Orbit pairs = {{initial_pair, {}}};
    NewOrbitElements next_elements = {pairs.begin()};

    Stopwatch produce_auto_orbit_time_watch;
    {
      auto produce_auto_orbit_time = produce_auto_orbit_time_watch.NewIter();

      auto i = 0;
      while (++i < 7) {
        produce_auto_orbit_time.Click();
        next_elements = ProduceAutomorhicOrbit(next_elements, &pairs);
        produce_auto_orbit_time.Click();
      }

      //auto old_size = 0u;
      //auto new_size = pairs.size();
      //while (old_size < new_size) {
      //  old_size = new_size;
      //  produce_auto_orbit_time.Click();
      //  next_elements = ProduceAutomorhicOrbit(next_elements, &pairs);
      //  produce_auto_orbit_time.Click();
      //  new_size = count_size_16_and_less(pairs);
      //  std::cout << ++i << ": " << new_size << std::endl;
      //}
    }

    std::cout << produce_auto_orbit_time_watch.average() << std::endl;

    std::cout << pairs.size() << std::endl;

    auto initial_pair_canonical = naive_get_canonical_pair::GetCanonicalPair(initial_pair.first, initial_pair.second);
    auto pair_count = 0u;
    for (auto&& pair : pairs) {
      ASSERT_EQ(initial_pair_canonical, naive_get_canonical_pair::GetCanonicalPair(pair.first.first, pair.first.second)) 
        << pair.first.first << ',' << pair.first.second << "; " 
        << initial_pair.first << ", " << initial_pair.second << "; "
        << ::testing::PrintToString(pair.second);
    }
/*
y->yXX
Y->xxY

x  y  y  x  Y, x  x  x  x  y  X  X  Y  x  y
x  yXXyXXxxxY, x  x  x  x  yXXX  X  xxYx  yXX
xyXXyxY, xxyXXYxy


*/

  }
}

} //namespace

} //namespace crag
