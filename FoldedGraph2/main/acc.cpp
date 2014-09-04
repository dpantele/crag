#include <algorithm>
#include <cassert>
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

//! Shifts word by @ref shift elements to the left
/**
 * Shift(2, {1, 2, 3, 4, 5}) is {3, 4, 5, 1, 2}
 */
void Shift(size_t shift, Word* w) {
  shift %= w->size();
  Flip(w->begin(), w->begin() + shift);
  Flip(w->begin() + shift, w->end());
  Flip(w->begin(), w->end());
}

inline size_t LabelToIndex(FoldedGraph2::Label l) {
  switch(l) {
    case 1:
      return 3 * 4;
    case -1:
      return 2 * 4;
    case 2:
      return 1 * 4;
    case -2:
      return 0 * 4;
    default:
      assert(false);
      return ~0u;
  }
}

inline size_t InvLabelToIndex(FoldedGraph2::Label l) {
  switch(l) {
    case 1:
      return 0 * 4;
    case -1:
      return 1 * 4;
    case 2:
      return 2 * 4;
    case -2:
      return 3 * 4;
    default:
      assert(false);
      return ~0u;
  }
}

//! Returns the minimal cyclic permutation of w
void PermuteToMin(Word* w) {

}

int main() {
  AllPairs pairs;
  AllPairs::value_type initial = {{1, 2, 2}, {2, 1, 1, 2, 2, 2, 2}};
  pairs.emplace_back(initial);
  pairs.emplace_back()



}

}  // namespace crag
