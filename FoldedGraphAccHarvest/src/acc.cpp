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

unsigned int GetMaxEqualLetterLength(Word* w) {
  auto max_equal_letter_length = 0u;
  auto current_equal_letter_length = 1u;
  for (auto shift = 0u; shift < w->size(); ++shift) {
    w->CyclicLeftShift();
    if (w->GetBack() == w->GetFront()) {
      ++current_equal_letter_length;
    } else {
      current_equal_letter_length = 1u;      
    }
    max_equal_letter_length = std::max(max_equal_letter_length, current_equal_letter_length);
  }

  return max_equal_letter_length;
}

Word Map(Word w, const std::array<unsigned int, 2 * Word::kAlphabetSize>& mapping) {
  Word result;
  while(!w.Empty()) {
    result.PushBack(mapping[w.GetFront()]);
    w.PopFront();
  }
  return result;
}

typedef std::array<unsigned int, 2 * Word::kAlphabetSize> Mapping;

Mapping MapToMin(Word* w) {
  Mapping result = {0, 1, 2, 3};
  if (w->Empty()) {
    return result;
  }

  auto count = 0u;
  while (w->GetBack() == w->GetFront() && count < w->size()) {
    w->CyclicRightShift();
    ++count;
  }

  if (count == w->size()) {
    *w = Word(w->size(), 0);
    result[w->GetFront()] = 0;
    result[FoldedGraph2::Inverse(w->GetFront())] = 1;
    if (w->GetFront() >= 2) {
      result[0] = 2;
      result[1] = 3;
    }
    return result;
  }

  auto candidate = *w;

  auto max_equal_letter_length = GetMaxEqualLetterLength(w);
  
  auto current_letter_length = 1u;

  for (auto shift = 0u; shift < w->size(); ++shift) {
    w->CyclicLeftShift();
    if (w->GetBack() == w->GetFront()) {
      ++current_letter_length;
    } else {
      if (current_letter_length == max_equal_letter_length) {
        auto y_preimage = w->GetFront();
        w->CyclicRightShift(max_equal_letter_length);
        auto x_preimage = w->GetFront();

        assert(x_preimage != y_preimage && 
          x_preimage != FoldedGraph2::Inverse(y_preimage));

        Mapping mapping = {0};
        mapping[x_preimage] = 0;
        mapping[FoldedGraph2::Inverse(x_preimage)] = 1;
        mapping[y_preimage] = 2;
        mapping[FoldedGraph2::Inverse(y_preimage)] = 3;

        auto new_candidate = Map(*w, mapping);
        if (new_candidate < candidate) {
          candidate = new_candidate;
          result = mapping;
        }
        w->CyclicLeftShift(max_equal_letter_length);
      }
      current_letter_length = 1u;
    }
  }
  assert(candidate.size() == w->size());
  assert(candidate.GetBack() != (candidate.GetFront() ^ 1));
  assert(candidate.GetBack() != candidate.GetFront());
  *w = candidate;
  return result;
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

void ReduceAndNormalize(Word* word) {
  *word = CyclicReduce(*word);
  auto inverse = *word;
  inverse.Invert();
  PermuteToMin(word);
  PermuteToMin(&inverse);
  if (inverse < *word) {
    *word = std::move(inverse);
  }
}

std::vector<Word> ReduceAndNormalize(std::vector<Word> words) {
  for (auto& word : words) {
    ReduceAndNormalize(&word);
  }

  std::sort(words.begin(), words.end());
  auto end = std::unique(words.begin(), words.end());
  words.erase(end, words.end());

  return words;
}

void ReduceAndNormalize(const Mapping& map, Word* word) {
  *word = CyclicReduce(*word);
  *word = Map(*word, map);
  auto inverse = *word;
  inverse.Invert();
  PermuteToMin(word);
  PermuteToMin(&inverse);
  if (inverse < *word) {
    *word = std::move(inverse);
  }
}

Mapping MapToMinWithInverse(Word* w) {
  *w = CyclicReduce(*w);
  auto mapping = MapToMin(w);
  auto w_inv = *w;
  w_inv.Invert();
  auto mapping_inv = MapToMin(&w_inv);
  if (w_inv < *w) {
    *w = w_inv;
    mapping = mapping_inv;
  }
  return mapping;
}

std::pair<Word, Word> GetCanonicalPair(const char* u_string, const char* v_string) {
  return GetCanonicalPair(Word(u_string), Word(v_string));
}

std::pair<Word, Word> GetCanonicalPair(Word u, Word v) {
  auto mapping = MapToMinWithInverse(&u);
  
  v = CyclicReduce(v);
  v = Map(v, mapping);
  PermuteToMin(&v);

  return std::make_pair(v, u);
}

void ReduceAndNormalize(Word* w, std::vector<Word>* words) {
  auto mapping = MapToMinWithInverse(w);
  for (auto& word : *words) {
    ReduceAndNormalize(mapping, &word);
  }

  std::sort(words->begin(), words->end());
  auto end = std::unique(words->begin(), words->end());
  words->erase(end, words->end());
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
