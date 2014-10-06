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

void PermuteToMinWithInverse(Word* w) {
  assert(w->GetFront() != FoldedGraph2::Inverse(w->GetBack()));

  auto w_inv = *w;
  w_inv.Invert();
  PermuteToMin(w);
  PermuteToMin(&w_inv);
  if (w_inv < *w) {
    *w = w_inv;
  }
}

void ReduceAndMinCycle(Word* word) {
  *word = CyclicReduce(*word);
  PermuteToMinWithInverse(word);
}

std::vector<Word> ReduceAndMinCycle(std::vector<Word> words) {
  for (auto& word : words) {
    ReduceAndMinCycle(&word);
  }

  std::sort(words.begin(), words.end());
  auto end = std::unique(words.begin(), words.end());
  words.erase(end, words.end());

  return words;
}

std::pair<Word, Word> GetCanonicalPair(Word u, Word v) {
  PermuteToMinWithInverse(&u);
  PermuteToMinWithInverse(&v);
  return std::make_pair(u, v);
}

std::pair<Word, Word> GetCanonicalPair(const char* u_string, const char* v_string) {
  return GetCanonicalPair(Word(u_string), Word(v_string));
}

void GetCanonicalPairs(Word* u, std::vector<Word>* vs) {
  PermuteToMinWithInverse(u);
  for (auto& v : *vs) {
    PermuteToMinWithInverse(&v);
  }

  std::sort(vs->begin(), vs->end());
  auto end = std::unique(vs->begin(), vs->end());
  vs->erase(end, vs->end());
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

uint64_t ipow(uint64_t base, uint64_t exp) {
  uint64_t result = 1;
  while (exp) {
    if (exp % 2) {
        result *= base;
    }
    exp /= 2;
    base *= base;
  }

  return result;
}

struct IterationVector {
  typedef unsigned int Symbol;

  std::vector<Symbol> current_;
  Symbol max_symbol_;
  uint64_t total_combinations_;

  explicit operator bool() const {
    return total_combinations_ != 0;
  }

  IterationVector(size_t length, Symbol max_symbol)
    : current_(length)
    , max_symbol_(max_symbol)
    , total_combinations_(ipow(max_symbol_, length))
  { }

  IterationVector& operator++() {
    if (total_combinations_ == 0) {
      return *this;
    }
    --total_combinations_;
    for (auto& num : current_) {
      ++num;
      if (num == max_symbol_) {
        num = 0;
      } else {
        break;
      }
    }
    return *this;
  }
};

std::set<Word> GenAllWords(unsigned int max_length) {
  std::set<Word> result = {{}};
  for (auto i = 1u; i <= max_length; ++i) {
    IterationVector v(i, Word::kAlphabetSize * 2);
    while (v) {
      Word w;
      for (const auto& letter : v.current_) {
        w.PushFront(letter);
      }
      result.emplace(w);
      ++v;
    }
  }
  return result;
}

} //namespace crag
