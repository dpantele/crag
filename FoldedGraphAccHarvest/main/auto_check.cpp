#include "acc.h"

#include <bitset>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <string>

using namespace crag;

using Label = FoldedGraph2::Label;

char LabelToChar(Label l) {
  return l % 2 ? 'X' + l / 2 : 'x' + l/2;
}

int main(int argc, const char *argv[]) {
  std::ifstream words("h12_c3_unproc_words.txt");

  std::set<std::pair<Word, Word>> all_words;

  std::pair<std::string, std::string> next_pair;

  words.ignore(64, '\n');

  while (words) {
    std::string appeared_in;
    words >> appeared_in;
    words >> next_pair.first;
    next_pair.first.pop_back();
    words >> next_pair.second;
    if (words) {
      all_words.emplace(Word(next_pair.first), Word(next_pair.second));
    }
  }

  std::pair<Word, Word> x_y_images[] = {
    {Word("xy"), Word("y")},
  };

  for (auto&& x_y_image : x_y_images) {
    Word x_image;
    Word y_image;

    std::tie(x_image, y_image) = x_y_image;

    std::array<Word, 4> images = {
      x_image,
      x_image,
      y_image,
      y_image,
    }; 

    images[1].Invert();
    images[3].Invert();

    for (Label l = 0; l < 4; ++l) {
      std::cout << LabelToChar(l) << " -> ";
      PrintTo(images[l], &std::cout);
      std::cout << ", ";
    }
    std::cout << std::endl;

    auto TransformWord = [&images](Word w) -> Word {
      Word result;
      while (!w.Empty()) {
        result.PushBack(images[w.GetFront()]);
        w.PopFront();
      }
      return result;
    };

    auto success_count = 0ull;
    auto overflow_count = 0ull;

    std::pair<Word, Word> transformed;
    for (auto&& word_pair : all_words) {
      try {
        transformed = std::make_pair(
            TransformWord(word_pair.first),
            TransformWord(word_pair.second)
        );
      } catch (const std::length_error&) {
        ++overflow_count;
        continue;
      }

      if (transformed.first.size() > 12 || transformed.second.size() > 12) {
        ++overflow_count;
        continue;
      }

      auto normalized = GetCanonicalPair(transformed.first, transformed.second);
      std::swap(normalized.first, normalized.second);
      if(all_words.count(normalized)) {
        ++success_count;
      } else {
        //std::cout << "\n=========\n";
        //PrintTo(word_pair.first, &std::cout);
        //std::cout << ", ";
        //PrintTo(word_pair.second, &std::cout);
        //std::cout << std::endl;
        //PrintTo(transformed.first, &std::cout);
        //std::cout << ", ";
        //PrintTo(transformed.second, &std::cout);
        //std::cout << std::endl;
        //PrintTo(normalized.first, &std::cout);
        //std::cout << ", ";
        //PrintTo(normalized.second, &std::cout);
        //std::cout << std::endl;
      }
    }

    std::cout << "Total: " << all_words.size() << ", success: " << success_count << "\n\n";
    std::cout << "Overflow: " << overflow_count << "\n\n";
  }
  return 0;
}
