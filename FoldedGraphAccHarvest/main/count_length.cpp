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
#include <memory>
#include <cstdio>
#include <cstring>

using namespace crag;

using Label = FoldedGraph2::Label;

char LabelToChar(Label l) {
  return l % 2 ? 'X' + l / 2 : 'x' + l/2;
}

#ifdef _MSC_VER
#define snprintf(BUFFER, SIZE, ...) _snprintf_s(BUFFER, SIZE, _TRUNCATE, __VA_ARGS__)
#endif 

int main(int argc, const char *argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: input_filename_pattern output_filename first_length last_length\nExample: h3_%d_unproc_words.txt h3_lengths_count.txt 10 16\n";
    exit(-1);
  }
  static const char* filename_template = argv[1];
  static const char* out_filename = argv[2];
  static const size_t filename_size = strlen(filename_template) + 5;
  
  static const auto first_length = std::stoi(argv[3]);
  static const auto last_length = std::stoi(argv[4]);

  std::set<size_t> total_lengths;
  std::map<std::string, std::map<size_t, uint64_t>> length_count;
  for (auto parameter = first_length; parameter <= last_length; ++parameter) {
    std::unique_ptr<char[]> filename(new char[filename_size]);
    snprintf(filename.get(), filename_size, filename_template, parameter);
    std::ifstream words_in(filename.get());

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
    length_count.emplace(std::to_string(parameter), std::move(length));
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
