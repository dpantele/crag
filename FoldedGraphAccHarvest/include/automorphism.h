//
// Created by dpantele on 8/27/15.
//

#ifndef CRAG_AUTOMORPHISM_H
#define CRAG_AUTOMORPHISM_H

#include "compressed_word.h"

namespace crag {

class Automorhpism {
 public:
  Automorhpism(CWord x_image, CWord y_image)
      : images_({x_image, x_image, y_image, y_image})
  {
    images_[1].Invert();
    images_[3].Invert();
  }

  static Automorhpism x_map(CWord x_image) {
    return Automorhpism(x_image, CWord("y"));
  }

  static Automorhpism y_map(CWord y_image) {
    return Automorhpism(CWord("x"), y_image);
  }

  CWord Apply(CWord w) const {
    CWord result;
    while (!w.Empty()) {
      result.PushBack(images_[w.GetFront()]);
      w.PopFront();
    }
    result.CyclicReduce();
    return result;
  }

  bool operator<(const Automorhpism& other) const {
    for (auto i = 0u; i + 1 < images_.size(); ++i) {
      if (images_[i] != other.images_[i]) {
        return images_[i] < other.images_[i];
      }
    }
    return images_.back() < other.images_.back();
  }
 private:
  std::array<CWord, 4> images_;
};

}

#endif //CRAG_AUTOMORPHISM_H
