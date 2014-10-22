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
#include <cstdio>

using namespace crag;

using Label = FoldedGraph2::Label;

char LabelToChar(Label l) {
  return l % 2 ? 'X' + l / 2 : 'x' + l/2;
}

#ifdef _MSC_VER
#define snprintf _snprintf_s
#endif 

int main(int argc, const char *argv[]) {
  static const char filename_template[] = "tb_c2_%s_unproc_words.txt";
  static const char out_filename[] = "tb_c2_lengths_count.txt";
  static const size_t filename_size = sizeof(filename_template) + 5;
  static const char * parameters[] = {"15", "16", "17", "18", "19", "20", "21", "22", "23"};//, "24"};

  std::set<size_t> total_lengths;
  std::map<std::string, std::map<size_t, uint64_t>> length_count;
  for (auto&& parameter : parameters) {
    char filename[filename_size] = {};
    snprintf(filename, filename_size, filename_template, parameter);
    std::ifstream words_in(filename);

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
        total_lengths.insert(next_pair.first.size() + next_pair.second.size());
      }
    }
    length_count.emplace(parameter, std::move(length));
  }

  std::ofstream out(out_filename);

  out << "  ";

  for (auto&& count : total_lengths) {
    out << std::setw(6) << count << " ";
  }
  out << std::endl;

  for (auto&& iteration : length_count) {
    out << iteration.first;
    for (auto&& count : iteration.second) {
      out << std::setw(6) << count.second << " ";
    }
    out << std::endl;
  }
  return 0;
}
