//
// Created by dpantele on 5/17/15.
//

#include <gtest/gtest.h>
#include "gtest/gtest.h"

#include "VanKampen.h"
#include "Word.h"

namespace crag { namespace free_group { namespace {

using Word = Word<int>;

TEST(VanKampen, Test1) {
  Word a = {1, 2, 1, -1};

  VanKampen<Word> d;
  d.Push(a);
  d.Push(Word({1, 1}));
  d.Push(Word({1, 2, -1}));
}

TEST(VanKampen, Test2) {
  const Word r = {1, 2, -1, -2};

  VanKampen<Word> d;
  d.Push(r.Conjugate(Word({-1, -1})));
  d.Push(r.Conjugate(Word({1, -2})));
  d.Push(r.Conjugate(Word({1, 1, -2})));
  d.Push(r.Conjugate(Word({1})));
  d.Push(r.Conjugate(Word({1, 1})));

  d.Fold();
  d.PrintTo(&std::cout);
}

TEST(VanKampen, Test3) {
  VanKampen<Word> d;
  d.Push(Word({1, 2}));
  d.Push(Word({1}));

  d.Fold();
  d.PrintTo(&std::cout);
}

TEST(VanKampen, Test4) {
  VanKampen<Word> d;
  d.Push(Word({1, 2}));
  d.Push(Word({-1}));
  d.PrintTo(&std::cout);

  d.Fold();
  d.PrintTo(&std::cout);
}

TEST(VanKampen, Test5) {
  VanKampen<Word> d;
  d.Push(Word({2, 1}));
  d.Push(Word({-1}));
//  d.Push(Word({2}));
//  d.Push(Word({1}));
//  d.Push(Word({-1, 2}));
  d.PrintTo(&std::cout);

  d.Fold();
  d.PrintTo(&std::cout);
}




}}}