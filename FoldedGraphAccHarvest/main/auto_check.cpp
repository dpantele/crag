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
    words >> next_pair.first;
    next_pair.first.pop_back();
    words >> next_pair.second;
    if (words) {
      all_words.emplace(Word(next_pair.first), Word(next_pair.second));
    }
  }

  std::pair<Label, Label> x_y_images[] = {
    {2, 0},
    {3, 0},
    {2, 1},
    {1, 2},
    {0, 3},
  };

  for (auto&& x_y_image : x_y_images) {
  
    Label x_image;
    Label y_image;

    std::tie(x_image, y_image) = x_y_image;

    std::array<Label, 4> images = {
      x_image,
      FoldedGraph2::Inverse(x_image),
      y_image,
      FoldedGraph2::Inverse(y_image),
    }; 

    for (Label l = 0; l < 4; ++l) {
      std::cout << LabelToChar(l) << " -> " << LabelToChar(images[l]) << ", ";
    }
    std::cout << std::endl;

    auto TransformWord = [&images](Word w) -> Word {
      Word result;
      while (!w.Empty()) {
        result.PushBack(images[w.GetFront()]);
        w.PopFront();
      }
      ReduceAndNormalize(&result);
      return result;
    };

    auto success_count = 0ull;

    for (auto&& word_pair : all_words) {
      auto transformed = std::make_pair(
          TransformWord(word_pair.first),
          TransformWord(word_pair.second)
      );
      if(all_words.count(transformed)) {
        ++success_count;
      } else {
        PrintTo(word_pair.first, &std::cout);
        std::cout << ", ";
        PrintTo(word_pair.second, &std::cout);
        std::cout << std::endl;
        PrintTo(transformed.first, &std::cout);
        std::cout << ", ";
        PrintTo(transformed.second, &std::cout);
        std::cout << std::endl;
      }
    }

    std::cout << "Total: " << all_words.size() << ", success: " << success_count << "\n\n";
  }
  return 0;
}
