//
// Created by dpantele on 5/17/15.
//

#ifndef CRAG_WORD_H
#define CRAG_WORD_H

#include <deque>
#include <stdexcept>
#include "../../../../../../usr/include/c++/4.9/bits/stl_deque.h"
#include "../../general/include/ConfigFile.h"

namespace crag { namespace free_group {

//! Element of a free group
template <typename LetterT>
class Word {
 public:
  using Letter = LetterT;

  Word() { }

  Word(std::initializer_list<Letter> i) {
    for (auto&& l : i) {
      PushBack(l);
    }
  }

  bool IsTrivial() const {
    return word_.empty();
  }

  void PushBack(Letter l) {
//    if (!IsTrivial() && word_.back() == -l) {
//      PopBack();
//    }
    word_.push_back(l);
  }

  void PushFront(Letter l) {
//    if (!IsTrivial() && word_.front() == -l) {
//      PopFront();
//    }
    word_.push_front(l);
  }

  Letter PopFront() {
    auto front_l = front();
    word_.pop_front();
    return front_l;
  }

  Letter PopBack() {
    auto back_l = back();
    word_.pop_back();
    return back_l;
  }

  Letter front() const {
    if (IsTrivial()) {
      throw std::runtime_error("Word is empty");
    }
    auto front = word_.front();
    return front;
  }

  Letter back() const {
    if (IsTrivial()) {
      throw std::runtime_error("Word is empty");
    }
    auto back = word_.back();
    return back;
  }

  //! Last letter becomes first
  void ShiftRight() {
    if (IsTrivial()) {
      return;
    }
    PushFront(back());
    PopBack();
  }

  //! First letter becomes last
  void ShiftLeft() {
    if (IsTrivial()) {
      return;
    }
    PushBack(front());
    PopFront();
  }

  Word Inverse() const {
    Word result_;
    for (auto&& l : word_) {
      result_.PushFront(-l);
    }
    return result_;
  }

  Word& ConcatenateWith(const Word& other) {
    for (auto&& l : other.word_) {
      PushBack(l);
    }
    return *this;
  }

  Word& ConjugateWith(const Word& other) {
    ConcatenateWith(other);
    *this = other.Inverse() * std::move(*this);
    return *this;
  }

  Word Conjugate(const Word& other) const {
    return Word(*this).ConjugateWith(other);
  }

  Word operator*(const Word& other) const& {
    Word copy = *this;
    copy.ConcatenateWith(other);
    return copy;
  }

  Word operator*(const Word& other) && {
    ConcatenateWith(other);
    return *this;
  }

  Word& Reduce() {
    if (word_.empty()) {
      return *this;
    }

    auto current_letter = std::next(word_.begin());
    auto current_end = current_letter;
    while (current_letter != word_.end()) {
      if (current_end == word_.begin()
        || -*current_letter != *std::prev(current_end)
      ) {
        //just push the letter
        *current_end = *current_letter;
        ++current_end;
      } else {
        // otherwise we 'pop' the last letter
        --current_end;
      }
      ++current_letter;
    }
    word_.erase(current_end, word_.end());
    return *this;
  }

  Word Reduced() const {
    Word copy = *this;
    return copy.Reduce();
  }

  void Invert() {
    *this = Inverse();
  }

  typename std::deque<Letter>::const_iterator begin() const {
    return word_.begin();
  }

  typename std::deque<Letter>::const_iterator end() const {
    return word_.end();
  }

  size_t length() const {
    return word_.size();
  }

  void PrintTo(std::ostream* out) {
    for (auto&& l : word_) {
      *out << l;
    }
  }

 private:
  std::deque<Letter> word_;
};

} //free_group
} //crag


#endif //CRAG_WORD_H
