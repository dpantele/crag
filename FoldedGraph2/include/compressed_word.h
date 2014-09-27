/**
* \file 
* @brief Short word of 2-alphabet stored in a single 64-bit field
*/

#include <cstdint>
#include <iostream>

namespace crag {

class CWord;

//! Stores word of lenght up to 16
class CWord {
public:
  typedef uint32_t size_type;
  static const int kAlphabetSize = 2;
  static const size_type kMaxLength = 16;

  CWord()
    : size_(0)
    , letters_(0)
  { }

  CWord(std::initializer_list<unsigned int> letters)
    : size_(0)
    , letters_(0)
  {
    assert(letters.size() <= kMaxLength);
    for (auto letter : letters) {
      assert(letter < kAlphabetSize * 2);
      if (!Empty() && GetBack() == (letter ^ 1)) {
        PopBack();
      } else {
        letters_ <<= kLetterShift;
        letters_ |= letter;
        ++size_;
      }
    }
  }

  explicit CWord(const std::string& letters) 
    : size_(0)
    , letters_(0)
  {
    for (auto&& sym : letters) {
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
          throw std::runtime_error("Only x, y, X, Y are allowed to be passed to CWord constructor");
      }
    }
  }

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
          throw std::runtime_error("Only x, y, X, Y are allowed to be passed to CWord constructor");
      }
    }
  }


  bool Empty() const {
    return size_ == 0;
  }

  void Clear() {
    size_ = 0;
    letters_ = 0;
  }

  uint32_t size() const {
    return size_;
  }

  void PushBack(unsigned int letter) {
    assert(letter < kAlphabetSize * 2);
    assert(size_ < kMaxLength);
    if(Empty() || (letter ^ 1) != GetBack()) {
      letters_ <<= kLetterShift;
      letters_ |= letter;
      ++size_;
    } else {
      PopBack();
    }
  }

  void PushFront(unsigned int letter) {
    assert(letter < kAlphabetSize * 2);
    assert(size_ < kMaxLength);
    if(Empty() || (letter ^ 1) != GetFront()) {
      letters_ |= (letter << (kLetterShift * size_));
      ++size_;
    } else {
      PopFront();
    }
  }

  void PopBack() {
    assert(!Empty());
    --size_;
    letters_ >>= kLetterShift;
  }

  void PopBack(size_type count) {
    assert(size() >= count);
    size_ -= count;
    letters_ >>= (kLetterShift * count);
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
    letters_ = ((letters_ >> 2) & 0x33333333u) | ((letters_ & 0x33333333u) << 2);
    // swap nibbles ... 
    letters_ = ((letters_ >> 4) & 0x0F0F0F0Fu) | ((letters_ & 0x0F0F0F0Fu) << 4);
    // swap bytes
    letters_ = ((letters_ >> 8) & 0x00FF00FFu) | ((letters_ & 0x00FF00FFu) << 8);
    // swap 2-byte long pairs
    letters_ = ( letters_ >> 16              ) | ( letters_               << 16);  
    // shift significant part back to the right
    letters_ >>= (sizeof(letters_) * 8 - size_ * kLetterShift);
  }

  void Invert() {
    static const uint32_t kInvertMask = 0x55555555u;
    Flip();
    letters_ ^= kInvertMask;
    ZeroUnused();
  }

  unsigned int GetFront() const {
    assert(!Empty());
    return letters_ >> (kLetterShift * (size_ - 1)); 
  }

  unsigned int GetBack() const {
    assert(!Empty());
    return letters_ & 3;
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
    ZeroUnused();
    letters_ |= left_part;
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
  uint32_t letters_;

  static const uint32_t kLetterMask = 3;
  static const uint32_t kLetterShift = 2;
  static const uint32_t kFullMask = ~uint32_t{0};

  void ZeroUnused() {
    if (size_) {
      letters_ &= (kFullMask >> (sizeof(kFullMask) * 8 - kLetterShift * size_));
    } else {
      letters_ = 0;
    }
  }

public:
  
};

inline void PrintTo(const CWord& w1, ::std::ostream* os) {
  auto w = w1;
  *os << w.size() << ": ";
  while (!w.Empty()) {
    *os << static_cast<char>((w.GetFront() % 2 ? 'X' : 'x') + w.GetFront() / 2);
    w.PopFront();
  }
}

inline std::ostream& operator<<(std::ostream& out, const CWord& w) {
  PrintTo(w, &out);
  return out;
}


}