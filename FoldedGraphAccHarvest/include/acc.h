#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>

#include "folded_graph2.h"

namespace crag {

typedef FoldedGraph2::Word Word;

typedef std::vector<std::pair<Word, Word>> AllPairs;

//! Flips sequence
template<typename InputIterator>
void Flip(InputIterator first, InputIterator last) {
  while (first != last) {
    --last;
    if (first != last) {
      std::swap(*first, *last);
      ++first;
    }
  }
}

//! Shifts sequnece one to the left
template<typename InputIterator>
void LeftShift(InputIterator first, InputIterator last) {
  if (first == last) {
    return;
  }
  auto front = std::move(*first);
  auto second = std::next(first);
  for (; second != last; ++second, ++first) {
    std::swap(*first, *second);
  }
  *first = std::move(front);
}


//! Shifts word by @ref shift elements to the left
/**
 * Shift(2, {1, 2, 3, 4, 5}) is {3, 4, 5, 1, 2}
 */
//void Shift(size_t shift, Word* w);

Word CyclicReduce(Word w);

//! Returns the minimal cyclic permutation of w
void PermuteToMin(Word* w);

//! Build the inverse of @ref w
template<typename InputIterator>
void Invert(InputIterator first, InputIterator last) {
  while (first != last) {
    --last;
    if (first != last) {
      std::swap(*first, *last);
      *first = FoldedGraph2::Inverse(*first);
      *last = FoldedGraph2::Inverse(*last);
      ++first;
    } else {
      *first = FoldedGraph2::Inverse(*first);
    }
  }
}

std::vector<Word> ReduceAndMinCycle(std::vector<Word>);

std::pair<Word, Word> GetCanonicalPair(const char* u, const char* v);
std::pair<Word, Word> GetCanonicalPair(Word u, Word v);
void GetCanonicalPairs(Word* u, std::vector<Word>* vs);

//! w -> s^(-1) w s
Word Conjugate(Word w, Word s);

uint64_t ipow(uint64_t base, uint64_t exp);

std::set<Word> GenAllWords(unsigned int max_length);
} //namespae crag