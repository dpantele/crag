/**
* @file
*
* @brief Definition of a basic letter for CWord
*/
#include <cassert>
#include <ostream>
#include <stdexcept>

#pragma once
#ifndef _CRAG_XY_LETTER_H_
#define _CRAG_XY_LETTER_H_

//! Small store which stores a x, y, X or Y
/**
* Mainly was created to have normal print and invert functions
*/
class XYLetter {
 public:
  static const unsigned short kAlphabetSize = 2;

  /*implicit*/ XYLetter(unsigned short letter)
    : letter_(letter)
  {
    assert(letter < 2 * kAlphabetSize);
  }

  /*implicit*/ XYLetter(int letter)
    : letter_(static_cast<unsigned short>(letter))
  { }

  /*implicit*/ XYLetter(unsigned int letter)
    : letter_(static_cast<unsigned short>(letter))
  { }

  explicit XYLetter(char letter)
    : letter_(CharToShort(letter))
  { }

  inline char AsChar() const {
    switch(letter_) {
      case 0:
        return 'x';
      case 1:
        return 'X';
      case 2:
        return 'y';
      case 3:
        return 'Y';
      default:
        break;
    }

    assert(letter_ < 4);
    return '\0';
  }

  inline unsigned int AsInt() const {
    return letter_;
  }

  inline XYLetter Inverse() const {
    return XYLetter(static_cast<unsigned short>(letter_ ^ 1));
  }

 private:
  unsigned short letter_;

  static inline unsigned short CharToShort(char letter) {
    switch(letter) {
      case 'x':
        return 0;
      case 'X':
        return 1;
      case 'y':
        return 2;
      case 'Y':
        return 3;
      default:
        throw std::invalid_argument("Only x, y, X, Y can be transformed into XYLetter");
    }
  }
};

//! Comparison operator to utilize implicit constructors
inline bool operator==(const XYLetter& lhs, const XYLetter& rhs) {
  return lhs.AsInt() == rhs.AsInt();
}

//! Comparison operator to utilize implicit constructors
inline bool operator!=(const XYLetter& lhs, const XYLetter& rhs) {
  return lhs.AsInt() != rhs.AsInt();
}

inline std::ostream& operator<<(std::ostream& out, const XYLetter& l) {
  return out << l.AsChar();
}

#endif //_CRAG_XY_LETTER_H_
