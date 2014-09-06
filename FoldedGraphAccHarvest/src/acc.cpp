#include "acc.h"

namespace crag {

void Shift(size_t shift, Word* w) {
  shift %= w->size();
  Flip(w->begin(), w->begin() + shift);
  Flip(w->begin() + shift, w->end());
  Flip(w->begin(), w->end());
}

Word CyclicReduce(const Word& w) {
  auto first = w.begin();
  auto last = w.end();
  while (first != last) {
    --last;
    if (first == last || FoldedGraph2::Inverse(*last) != *first) {
      return Word(first, ++last);
    }
    ++first;
  }
  return {};
}

void PermuteToMin(Word* w) {
  if (w->empty()) return;

  auto current_permutation = *w;
  for (auto i = 0u; i < w->size() - 1; ++i) {
    LeftShift(current_permutation.begin(), current_permutation.end());
    if (*w > current_permutation) {
      *w = current_permutation;
    }
  }
}

std::vector<Word> ReduceAndNormalize(std::vector<Word> words) {
  for (auto& word : words) {
    word = CyclicReduce(word);
    auto inverse = word;
    Invert(inverse.begin(), inverse.end());
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

Word Conjugate(Word w, const Word& s) {
  auto result = s;
  Invert(result.begin(), result.end());
  result.insert(result.end(), w.begin(), w.end());
  result.insert(result.end(), s.begin(), s.end());
  w = std::move(result);
  result.clear();
  for (const auto& letter : w) {
    if (result.empty()) {
      result.push_back(letter);
      continue;
    }
    if (FoldedGraph2::Inverse(letter) == result.back()) {
      result.pop_back();
    } else {
      result.push_back(letter);
    }
  }
  return result;
}

} //namespace crag
