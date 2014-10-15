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
  std::ifstream normalized_words_in("aut2_unproc_words.txt");

  std::set<std::pair<Word, Word>> all_words;

  std::pair<std::string, std::string> next_pair;

  normalized_words_in.ignore(64, '\n');

  while (normalized_words_in) {
    std::string appeared_in;
    normalized_words_in >> appeared_in;
    normalized_words_in >> next_pair.first;
    next_pair.first.pop_back();
    normalized_words_in >> next_pair.second;
    if (normalized_words_in) {
      all_words.emplace(Word(next_pair.first), Word(next_pair.second));
    }
  }

  std::cout << "Loaded " << all_words.size() << " pairs" << std::endl;

  std::ifstream checked_words_in("h12_unproc_words.txt");

  auto total_count = 0u;
  auto half_success_count = 0u;
  auto success_count = 0u;

  checked_words_in.ignore(64, '\n');

  while (checked_words_in) {
    std::string appeared_in;
    checked_words_in >> appeared_in;
    checked_words_in >> next_pair.first;
    next_pair.first.pop_back();
    checked_words_in >> next_pair.second;
    if (checked_words_in) {
      ++total_count;
      auto checked_pair = GetCanonicalPair(next_pair.first.c_str(), next_pair.second.c_str());
      auto swapped = GetCanonicalPair(next_pair.second.c_str(), next_pair.first.c_str());
      if (all_words.count(checked_pair)) {
        ++success_count;
      } else if (all_words.count(checked_pair) || all_words.count(swapped)) {
        ++half_success_count;
      } else {
        
      }
    }
  }

  std::cout << "Total: " << total_count << std::endl;
  std::cout << "Success: " << success_count << std::endl;
  std::cout << "With swapped: " << half_success_count << std::endl;
  return 0;
}
