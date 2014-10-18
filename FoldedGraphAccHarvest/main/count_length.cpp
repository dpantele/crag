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
#include <map>

using namespace crag;

using Label = FoldedGraph2::Label;

char LabelToChar(Label l) {
  return l % 2 ? 'X' + l / 2 : 'x' + l/2;
}

int main(int argc, const char *argv[]) {
  std::string prefix = "h16";

  std::ifstream words_in(prefix + "_unproc_words.txt");

  std::map<size_t, uint64_t> length;

  std::pair<std::string, std::string> next_pair;

  words_in.ignore(64, '\n');

  while (words_in) {
    std::string appeared_in;
    words_in >> appeared_in;
    words_in >> next_pair.first;
    next_pair.first.pop_back();
    words_in >> next_pair.second;
    if (words_in) {
      ++length[next_pair.first.size() + next_pair.second.size()];
    }
  }

  std::ofstream out(prefix + "_count.txt");

  for (auto&& count : length) {
    out << std::setw(6) << count.first << " ";
  }
  out << std::endl;

  for (auto&& count : length) {
    out << std::setw(6) << count.second << " ";
  }
  out << std::endl;
  return 0;
}
