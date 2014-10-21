/**
* \file 
* @brief Short word of 2-alphabet stored in a single 64-bit field
*/

#pragma once
#ifndef CRAG_COMPRESSED_WORD_H_
#define CRAG_COMPRESSED_WORD_H_

#include <cstdint>
#include <iostream>
#include <random>

namespace crag {

class CWord;

//! Stores word of lenght up to 32
class CWord {
public:
  typedef unsigned short size_type;
  static const unsigned short kAlphabetSize = 2;
  static const size_type kMaxLength = 32;

  CWord()
    : size_(0)
    , letters_(0)
  { }

  CWord(std::initializer_list<unsigned int> letters)
    : size_(0)
    , letters_(0)
  {
    for (auto letter : letters) {
      assert(letter < kAlphabetSize * 2);
      if (!Empty() && GetBack() == (letter ^ 1)) {
        PopBack();
      } else {
        if (size_ == kMaxLength) {
          throw std::length_error("Length of CWord is limited by 32");
        }
        letters_ <<= kLetterShift;
        letters_ |= letter;
        ++size_;
      }
    }
    assert((letters_ >> (kLetterShift * size_)) == 0);
  }

  CWord(size_t count, unsigned int letter) 
    : size_(0)
    , letters_(0)
  {
    while (count > 0) {
      --count;
      PushBack(letter);
    }
    assert((letters_ >> (kLetterShift * size_)) == 0);
  }

  explicit CWord(const std::string& letters) 
    : CWord(letters.c_str())
  { }

  explicit CWord(const char* letters) 
    : size_(0)
    , letters_(0)
  {
    for (; *letters != 0; ++letters) {
      auto sym = *letters;
      switch (sym) {
        case 'x':
          PushBack(0);
          break;
        case 'X':
          PushBack(1);
          break;
        case 'y':
          PushBack(2);
          break;
        case 'Y':
          PushBack(3);
          break;
        default:
          throw std::invalid_argument("Only x, y, X, Y are allowed to be passed to CWord constructor");
      }
    }
    assert((letters_ >> (kLetterShift * size_)) == 0);
  }

  bool Empty() const {
    return size_ == 0;
  }

  void Clear() {
    size_ = 0;
    letters_ = 0;
  }

  size_type size() const {
    return size_;
  }

  void PushBack(uint64_t letter) {
    assert(letter < kAlphabetSize * 2);
    if(Empty() || (letter ^ 1) != GetBack()) {
      if (size_ == kMaxLength) {
        throw std::length_error("Length of CWord is limited by 32");
      }
      letters_ <<= kLetterShift;
      letters_ |= letter;
      ++size_;
    } else {
      PopBack();
    }
    assert(size_ == kMaxLength || (letters_ >> (kLetterShift * size_)) == 0);
  }

  void PushBack(CWord w) {
    while (!w.Empty()) {
      PushBack(w.GetFront());
      w.PopFront();
    }
    assert(size_ == kMaxLength || (letters_ >> (kLetterShift * size_)) == 0);
  }

  void PushFront(uint64_t letter) {
    assert(letter < kAlphabetSize * 2);
    if(Empty() || (letter ^ 1) != GetFront()) {
      if (size_ == kMaxLength) {
        throw std::length_error("Length of CWord is limited by 32");
      }
      letters_ |= (letter << (kLetterShift * size_));
      ++size_;
    } else {
      PopFront();
    }
    assert(size_ == kMaxLength || (letters_ >> (kLetterShift * size_)) == 0);
  }

  void PopBack() {
    assert(!Empty());
    --size_;
    letters_ >>= kLetterShift;
    assert(size_ == kMaxLength || (letters_ >> (kLetterShift * size_)) == 0);
  }

  void PopBack(size_type count) {
    assert(size() >= count);
    size_ -= count;
    letters_ >>= (kLetterShift * count);
    assert(size_ == kMaxLength || (letters_ >> (kLetterShift * size_)) == 0);
  }

  void PopFront() {
    assert(!Empty());
    --size_;
    ZeroUnused();
  }

  void PopFront(size_type count) {
    assert(size() >= count);
    size_ -= count;
    ZeroUnused();
  }

  void Flip() {
    // swap consecutive pairs
    letters_ = ((letters_ >> 2) & 0x3333333333333333ull) | ((letters_ & 0x3333333333333333ull) << 2);
    // swap nibbles ... 
    letters_ = ((letters_ >> 4) & 0x0F0F0F0F0F0F0F0Full) | ((letters_ & 0x0F0F0F0F0F0F0F0Full) << 4);
    // swap bytes
    letters_ = ((letters_ >> 8) & 0x00FF00FF00FF00FFull) | ((letters_ & 0x00FF00FF00FF00FFull) << 8);
    // swap 2-byte long pairs
    letters_ = ( letters_ >> 16 & 0x0000FFFF0000FFFFull) | ((letters_ & 0x0000FFFF0000FFFFull) << 16);  
    // swap 4-byte long pairs
    letters_ = ( letters_ >> 32                        ) | ((letters_                        ) << 32);  
    // shift significant part back to the right
    letters_ >>= (sizeof(letters_) * 8 - size_ * kLetterShift);
  }

  void Invert() {
    static const uint64_t kInvertMask = 0x5555555555555555ull;
    Flip();
    letters_ ^= kInvertMask;
    ZeroUnused();
  }

  unsigned int GetFront() const {
    assert(!Empty());
    return static_cast<unsigned int>(letters_ >> (kLetterShift * (size_ - 1))); 
  }

  unsigned int GetBack() const {
    assert(!Empty());
    return static_cast<unsigned int>(letters_ & kLetterMask);
  }

  bool operator < (const CWord& other) const {
    return size_ == other.size_ ?  letters_ < other.letters_ : size_ < other.size_;
  }

  bool operator == (const CWord& other) const {
    return letters_ == other.letters_ && size_ == other.size();
  }

  bool operator != (const CWord& other) const {
    return !(*this == other);
  }

  void CyclicLeftShift(size_type shift = 1) {
    shift %= size_;
    auto left_part = (letters_ >> (kLetterShift * (size_ - shift)));
    letters_ <<= kLetterShift * shift;
    letters_ |= left_part;
    ZeroUnused();
  }

  void CyclicRightShift(size_type shift = 1) {
    CyclicLeftShift(size_ - shift);
  }

  void CyclicReduce() {
    while(!Empty() && (GetFront() == (GetBack() ^ 1))) {
      PopFront();
      PopBack();
    }
  }

private:
  size_type size_;
  uint64_t letters_;

  static const uint64_t kLetterMask = 3;
  static const uint64_t kLetterShift = 2;
  static const uint64_t kFullMask = ~uint64_t{0};

  void ZeroUnused() {
    if (size_) {
      letters_ &= (kFullMask >> (sizeof(kFullMask) * 8 - kLetterShift * size_));
    } else {
      letters_ = 0;
    }
    assert(size_ == kMaxLength || (letters_ >> (kLetterShift * size_)) == 0);
  }

};

inline void PrintWord(const CWord& w1, std::ostream* out) {
  auto w = w1;
  while (!w.Empty()) {
    *out << static_cast<char>((w.GetFront() % 2 ? 'X' : 'x') + w.GetFront() / 2);
    w.PopFront();
  }
}

//gtest debugging print
inline void PrintTo(const CWord& w, ::std::ostream* out) {
  *out << w.size() << ": ";
  PrintWord(w, out);
}

inline std::ostream& operator<<(std::ostream& out, const CWord& w) {
  PrintTo(w, &out);
  return out;
}

class RandomWord {
public:
  RandomWord(size_t min_size, size_t max_size)
    : random_letter_(0, 2 * CWord::kAlphabetSize - 1)
    , random_length_(min_size, max_size)
  { }

  template<class RandomEngine>
  CWord operator()(RandomEngine& engine) {
    CWord w;
    size_t length = random_length_(engine);
    while(w.size() < length) {
      w.PushBack(random_letter_(engine));
    }
    return w;
  }

private:
  std::uniform_int_distribution<> random_letter_;
  std::uniform_int_distribution<size_t> random_length_;


};


}

#endif //CRAG_COMPREESED_WORDS_H_