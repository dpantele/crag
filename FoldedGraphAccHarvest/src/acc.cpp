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

Word Map(Word w, const Mapping& mapping) {
  Word result;
  while(!w.Empty()) {
    result.PushBack(mapping[w.GetFront()]);
    w.PopFront();
  }
  return result;
}

Mapping MapToMin(Word* w) {
  Mapping result = {Word(1, 0), Word(1, 1), Word(1, 2), Word(1, 3)};
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
    result[w->GetFront()] = Word(1, 0);
    result[FoldedGraph2::Inverse(w->GetFront())] = Word(1, 1);
    if (w->GetFront() >= 2) {
      result[0] = Word(1, 2);
      result[1] = Word(1, 3);
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

        Mapping mapping;
        mapping[x_preimage] = Word(1, 0);
        mapping[FoldedGraph2::Inverse(x_preimage)] = Word(1, 1);
        mapping[y_preimage] = Word(1, 2);
        mapping[FoldedGraph2::Inverse(y_preimage)] = Word(1, 3);

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
  //PrintTo(*w, &std::cout);
  //std::cout << std::endl;
  for (auto i = 0u; i < w->size() - 1; ++i) {
    current_permutation.CyclicLeftShift();
    if (current_permutation < *w) {
      *w = current_permutation;
    }
  }
}

void PermuteToMinWithInverse(Word* w) {
  assert(w->Empty() || w->GetFront() != FoldedGraph2::Inverse(w->GetBack()));

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

Mapping MapToMinWithInverse(Word* w) {
  *w = CyclicReduce(*w);
  auto w_inv = *w;
  w_inv.Invert();

  auto mapping = MapToMin(w);
  auto mapping_inv = MapToMin(&w_inv);
  if (w_inv < *w) {
    *w = w_inv;
    mapping = mapping_inv;
  }
  return mapping;
}

void ReduceMapAndMinCycle(const Mapping& mapping, Word* w) {
  *w = Map(*w, mapping);
  PrintTo(*w, &std::cout);
  std::cout << std::endl;
  *w = CyclicReduce(*w);
  PermuteToMinWithInverse(w);
}

Word Inverse(Word w) {
  w.Invert();
  return w;
}

std::set<std::pair<Word, Word>> MinimizeTotalLength(Word u, Word v, size_t max_length) {
#define MORPHISM(X, Y) {Word(X), Inverse(Word(X)), Word(Y), Inverse(Word(Y))}
  static const std::array<Word, 2 * FoldedGraph2::kAlphabetSize> morphisms[] = {
    MORPHISM("yx", "y"),
    MORPHISM("Yx", "y"),
    MORPHISM("xy", "y"),
    MORPHISM("xY", "y"),
    MORPHISM("yxY", "y"),
    MORPHISM("Yxy", "y"),
    MORPHISM("x", "yx"),
    MORPHISM("x", "yX"),
    MORPHISM("x", "xy"),
    MORPHISM("x", "Xy"),
    MORPHISM("x", "Xyx"),
    MORPHISM("x", "xyX"),
  };
#undef MORPHISM

  std::set<std::pair<Word, Word>> min_length_pairs = {std::make_pair(u, v)};

  bool progress = true;
  while (progress) {
    progress = false;
    for (auto&& morphism : morphisms) {
      try {
        auto u_image = CyclicReduce(Map(u, morphism));
        auto v_image = CyclicReduce(Map(v, morphism));
        if (max_length == 0 || (u_image.size() <= max_length && v_image.size() <= max_length)) {
          if (u_image.size() + v_image.size() < u.size() + v.size()) {
            u = u_image;
            v = v_image;
            progress = true;
            min_length_pairs.clear();
            min_length_pairs.emplace(u_image, v_image);
            break;
          } else if (u_image.size() + v_image.size() == u.size() + v.size()) {
            min_length_pairs.emplace(u_image, v_image);
          }
        }
      } catch(std::length_error&) {
        continue;
      }
    }
  }

  //now apply automorphisms to all pairs trying find all same-length pairs
  std::vector<std::pair<Word, Word>> pairs_to_check(min_length_pairs.begin(), min_length_pairs.end());

  while (!pairs_to_check.empty()) {
    std::vector<std::pair<Word, Word>> new_pairs;
    for (auto&& morphism : morphisms) {
      for (auto&& pair : pairs_to_check) {
        try {
          auto u_image = CyclicReduce(Map(pair.first, morphism));
          auto v_image = CyclicReduce(Map(pair.second, morphism));
          if (max_length == 0 || (u_image.size() <= max_length && v_image.size() <= max_length)) {
            assert(u_image.size() + v_image.size() >= pair.first.size() + pair.second.size());
            if (u_image.size() + v_image.size() == pair.first.size() + pair.second.size()) {
              if (min_length_pairs.emplace(u_image, v_image).second) {
                new_pairs.emplace_back(u_image, v_image);
              }
            }
          }
        } catch(std::length_error&) {
          continue;
        }
      }
    }
    pairs_to_check = std::move(new_pairs);
  }
  return min_length_pairs;
}

std::pair<Word, Word> GetCanonicalPair(Word u, Word v, size_t max_length) {
  u = CyclicReduce(u);
  v = CyclicReduce(v);

  const auto min_length_pairs = MinimizeTotalLength(u, v, max_length);
  assert(!min_length_pairs.empty());

  //since min_length_pairs is a map, the first element has the shortest u
  auto min_size = min_length_pairs.begin()->first.size();

  for (auto&& pair : min_length_pairs) {
    min_size = std::min(min_size, pair.second.size());
  }

  std::tie(u, v) = *min_length_pairs.begin();
  PermuteToMinWithInverse(&u);
  PermuteToMinWithInverse(&v);
  for (auto&& uv : min_length_pairs) {
    if (uv.first.size() != min_size && uv.second.size() != min_size) {
      continue;
    }
    auto up = uv.first;
    auto vp = uv.second;
    if (up.size() > vp.size()) {
      std::swap(up, vp);
    }

    std::cout << "!!!!!!!!\n";
    PrintTo(up, &std::cout);
    std::cout << std::endl;
    auto u_min_mapping = MapToMinWithInverse(&up);
    std::cout << "!!!!!!!!\n";
    PrintTo(up, &std::cout);
    std::cout << std::endl;

    if (vp.size() == up.size()) {
      auto v_min_mapping = MapToMinWithInverse(&vp);

      if (vp < up) {
        up = vp;
        vp = uv.first;
        ReduceMapAndMinCycle(v_min_mapping, &vp);
      } else {
        vp = uv.second;
        ReduceMapAndMinCycle(u_min_mapping, &vp);
      }
    } else {
      ReduceMapAndMinCycle(u_min_mapping, &vp);
    } 
    
    if (up < u || (up == u && vp < v)) {
      PrintTo(vp, &std::cout);
      vp.Invert();
      PrintTo(vp, &std::cout);
      vp.Invert();
      u = up;
      v = vp;
    }
  }
  return std::make_pair(u, v);
}

//std::pair<Word, Word> GetCanonicalPair(Word u, Word v, size_t max_length) {
//  MapToMinWithInverse(&u);
//  MapToMinWithInverse(&v);
//  if (v < u) {
//    std::swap(u, v);
//  }
//  auto min_length_pairs = MinimizeTotalLength(u, v, max_length);
//
//  for (auto&& uv : min_length_pairs) {
//    auto u_min = uv.first;
//    auto u_min_mapping = MapToMinWithInverse(&u_min);
//
//    auto v_min = uv.second;
//    auto v_min_mapping = MapToMinWithInverse(&v_min);
//
//    Word u_mapped = uv.first;
//    Word v_mapped = uv.second;
//    if (v_min < u_min) {
//      u_mapped = v_min;
//      v_mapped = uv.first;
//      ReduceMapAndMinCycle(v_min_mapping, &v_mapped);
//    } else {
//      u_mapped = u_min;
//      v_mapped = uv.second;
//      ReduceMapAndMinCycle(u_min_mapping, &v_mapped);
//    }
//    if (u_mapped < u || (u_mapped == u && v_mapped < v)) {
//      u = u_mapped;
//      v = v_mapped;
//    }
//  }
//  return std::make_pair(u, v);
//}

std::pair<Word, Word> GetCanonicalPair(const char* u_string, const char* v_string, size_t max_length) {
  return GetCanonicalPair(Word(u_string), Word(v_string));
}

void GetCanonicalPairs(Word* u, std::vector<Word>* vs) {
  auto mapping = MapToMinWithInverse(u);
  for (auto& v : *vs) {
    ReduceMapAndMinCycle(mapping, &v);
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
