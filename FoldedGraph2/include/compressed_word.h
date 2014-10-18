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

//! Stores word of lenght up to 64
class CWord {
public:
  typedef uint32_t size_type;
  static const int kAlphabetSize = 2;
  static const size_type kMaxLength = 64;

  CWord()
    : size_(0)
    , letters_(0, 0)
  { }

  CWord(std::initializer_list<unsigned int> letters)
    : size_(0)
    , letters_(0, 0)
  {
    for (auto letter : letters) {
      assert(letter < kAlphabetSize * 2);
      if (!Empty() && GetBack() == (letter ^ 1)) {
        PopBack();
      } else {
        if (size_ == kMaxLength) {
          throw std::length_error("Length of CWord is limited by 32");
        }
        letters_.first <<= kLetterShift;
        letters_.first |= (letters_.second >> 62);
        letters_.second <<= kLetterShift;
        letters_.second |= letter;
        ++size_;
      }
    }
    assert(CheckZeroed());
  }

  CWord(size_t count, unsigned int letter) 
    : size_(0)
    , letters_(0, 0)
  {
    while (count > 0) {
      --count;
      PushBack(letter);
    }
    assert(CheckZeroed());
  }

  explicit CWord(const std::string& letters) 
    : CWord(letters.c_str())
  { }

  explicit CWord(const char* letters) 
    : size_(0)
    , letters_(0, 0)
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
    assert(CheckZeroed());
  }

  bool Empty() const {
    return size_ == 0;
  }

  void Clear() {
    size_ = 0;
    letters_.first = 0;
    letters_.second = 0;
  }

  uint32_t size() const {
    return size_;
  }

  void PushBack(uint64_t letter) {
    assert(letter < kAlphabetSize * 2);
    if(Empty() || (letter ^ 1) != GetBack()) {
      if (size_ == kMaxLength) {
        throw std::length_error("Length of CWord is limited by 32");
      }
      letters_.first <<= kLetterShift;
      letters_.first |= (letters_.second >> 62);
      letters_.second <<= kLetterShift;
      letters_.second |= letter;
      ++size_;
    } else {
      PopBack();
    }
    assert(CheckZeroed());
  }

  void PushBack(CWord w) {
    while (!w.Empty()) {
      PushBack(w.GetFront());
      w.PopFront();
    }
    assert(CheckZeroed());
  }

  void PushFront(uint64_t letter) {
    assert(letter < kAlphabetSize * 2);
    if(Empty() || (letter ^ 1) != GetFront()) {
      if (size_ == kMaxLength) {
        throw std::length_error("Length of CWord is limited by 32");
      }
      if (size_ >= 32) {
        letters_.first |= (letter << (kLetterShift * (size_ - 32)));
      } else {
        letters_.second |= (letter << (kLetterShift * size_ ));
      }
      ++size_;
    } else {
      PopFront();
    }
    assert(CheckZeroed());
  }

  void PopBack() {
    assert(!Empty());
    RightLetterwiseShift(1);
    --size_;
    assert(CheckZeroed());
  }

  void PopBack(size_type count) {
    assert(size() >= count);
    if (count == 0) {
      return;
    }

    RightLetterwiseShift(count);
    size_ -= count;
    assert(CheckZeroed());
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
    letters_.second = Flip(letters_.second);

    // shift significant part back to the right
    if (size_ <= 32) {
      letters_.second >>= (64 - size_ * kLetterShift);
    } else {
      letters_.first = Flip(letters_.first);
      std::swap(letters_.first, letters_.second);
      RightLetterwiseShift(kMaxLength - size_);
    }
  }

  void Invert() {
    static const uint64_t kInvertMask = 0x5555555555555555ull;
    Flip();
    letters_.first ^= kInvertMask;
    letters_.second ^= kInvertMask;
    ZeroUnused();
  }

  unsigned int GetFront() const {
    assert(!Empty());
    if (size_ > 32) {
      return static_cast<unsigned int>(letters_.first >> (kLetterShift * (size_ - 33))); 
    } else {
      return static_cast<unsigned int>(letters_.second >> (kLetterShift * (size_ - 1))); 
    }
  }

  unsigned int GetBack() const {
    assert(!Empty());
    return static_cast<unsigned int>(letters_.second & kLetterMask);
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
    if (shift == 0) {
      return;
    }
    if (shift < 0) {
      shift += size_;
    }
    if (size_ <= 32) {
      auto left_part = (letters_.second >> (kLetterShift * (size_ - shift)));
      letters_.second <<= kLetterShift * shift;
      letters_.second |= left_part;
      ZeroUnused();
    }
    if (size_ > 32) {
      if (shift > 32) {
        CyclicLeftShift(32);
        CyclicLeftShift(shift - 32);
      } else if (shift == 32) {
        std::swap(letters_.first, letters_.second);
        letters_.second <<= (64 - size_) * kLetterShift;
        letters_.second |= (letters_.first >> (size_ - 32) * kLetterShift);
        ZeroUnused();
      } else {
        auto first_left_part = (letters_.first >> (kLetterShift * (size_ - shift - 32)));
        letters_.first <<= (kLetterShift * shift);
        letters_.first |= (letters_.second >> (kLetterShift * (32 - shift)));
        letters_.second <<= kLetterShift * shift;
        letters_.second |= first_left_part;
        ZeroUnused();
      }
    }
    assert(CheckZeroed());
  }

  /*
  
   1122, 334455
1  2233, 445511
2  3344, 551122
3  4455, 112233
4  5511, 223344
  */

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
  std::pair<uint64_t, uint64_t> letters_;

  static const uint64_t kLetterMask = 3;
  static const uint64_t kLetterShift = 2;
  static const uint64_t kFullMask = ~uint64_t{0};

  void ZeroUnused() {
    if (size_ > 32) {
      letters_.first &= (kFullMask >> (sizeof(kFullMask) * 8 - kLetterShift * (size_ - 32)));
    } else if (size_) {
      letters_.first = 0;
      letters_.second &= (kFullMask >> (sizeof(kFullMask) * 8 - kLetterShift * size_));
    } else {
      letters_.first = 0;
      letters_.second = 0;
    }
    assert(CheckZeroed());
  }

  bool CheckZeroed() const {
    if (size_ == 64) {
      return true;
    }
    if (size_ > 32) {
      return (letters_.first >> (kLetterShift * (size_ - 32))) == 0;
    }
    if (size_ == 32) {
      return letters_.first == 0;
    }
    return (letters_.second >> (kLetterShift * size_)) == 0;
  }

  static uint64_t Flip(uint64_t letters) {
    // swap consecutive pairs
    letters = ((letters >> 2) & 0x3333333333333333ull) | ((letters & 0x3333333333333333ull) << 2);
    // swap nibbles ... 
    letters = ((letters >> 4) & 0x0F0F0F0F0F0F0F0Full) | ((letters & 0x0F0F0F0F0F0F0F0Full) << 4);
    // swap bytes
    letters = ((letters >> 8) & 0x00FF00FF00FF00FFull) | ((letters & 0x00FF00FF00FF00FFull) << 8);
    // swap 2-byte long pairs
    letters = ( letters >> 16 & 0x0000FFFF0000FFFFull) | ((letters & 0x0000FFFF0000FFFFull) << 16);  
    // swap 4-byte long pairs
    letters = ( letters >> 32                        ) | ((letters                        ) << 32);  
    return letters;
  }
  
  void RightLetterwiseShift(size_type count) {
    if (count >= 32) {
      letters_.second = letters_.first;
      letters_.first = 0;
      letters_.second >>= (kLetterShift * (count - 32));
    } else {
      letters_.second >>= (kLetterShift * count);
      letters_.second |= ((letters_.first & (kFullMask >> (32 - count) * kLetterShift)) << (32 - count) * kLetterShift);
      letters_.first >>= (kLetterShift * count);
    }
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