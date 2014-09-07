#include "acc.h"

namespace crag {

//void Shift(size_t shift, Word* w) {
//  shift %= w->size();
//  Flip(w->begin(), w->begin() + shift);
//  Flip(w->begin() + shift, w->end());
//  Flip(w->begin(), w->end());
//}
//
Word CyclicReduce(Word w) {
  while (w.size() > 1) {
    if (w.GetFront() == FoldedGraph2::Inverse(w.GetBack())) {
      w.PopBack();
      w.PopFront();
    } else {
      break;
    }
  }
  return w;
}

void PermuteToMin(Word* w) {
  if (w->Empty()) return;

  auto current_permutation = *w;
  for (auto i = 0u; i < w->size() - 1; ++i) {
    current_permutation.CyclicLeftShift();
    if (current_permutation < *w) {
      *w = current_permutation;
    }
  }
}

std::vector<Word> ReduceAndNormalize(std::vector<Word> words) {
  for (auto& word : words) {
    word = CyclicReduce(word);
    auto inverse = word;
    inverse.Invert();
    PermuteToMin(&word);
    PermuteToMin(&inverse);
    if (inverse < word) {
      word = std::move(inverse);
    }
  }

  std::sort(words.begin(), words.end());
  auto end = std::unique(words.begin(), words.end());
  words.erase(end, words.end());

  return words;
}

Word Conjugate(Word w, Word s) {
  auto s_inv = s;
  s_inv.Invert();
  
  while (!s_inv.Empty() && !w.Empty() && s_inv.GetBack() == FoldedGraph2::Inverse(w.GetFront())) {
    s_inv.PopBack();
    w.PopFront();
  }

  while(!s_inv.Empty()) {
    if (w.size() == Word::kMaxLength) {
      return Word{};
    }
    w.PushFront(s_inv.GetBack());
    s_inv.PopBack();
  }

  while (!s.Empty() && !w.Empty() && s.GetFront() == FoldedGraph2::Inverse(w.GetBack())) {
    s.PopFront();
    w.PopBack();
  }

  while(!s.Empty()) {
    if (w.size() == Word::kMaxLength) {
      return Word{};
    }
    w.PushBack(s.GetFront());
    s.PopFront();
  }

  return w;
}

} //namespace crag
