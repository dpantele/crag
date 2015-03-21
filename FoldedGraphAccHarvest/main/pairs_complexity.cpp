#include "acc.h"
#include "folded_graph2.h"

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

unsigned int PairComplexity(const Word& u, const Word& v) {
  FoldedGraph2 g;

  auto IsSingleVertex = [&g]() -> bool {
    static const CWord x = {0};
    static const CWord y = {2};

    auto x_trace = g.ReadWord(x);
    if (std::get<0>(x_trace) != g.root() || std::get<1>(x_trace) != 1) {
      return false;
    }

    auto y_trace = g.ReadWord(y);
    if (std::get<0>(y_trace) != g.root() || std::get<1>(y_trace) != 1) {
      return false;
    }

    return true;
  };

  auto complexity = 0u;
  for (; !IsSingleVertex(); ++complexity) {
    auto g_size = g.size();
    g.CompleteWith(u, g_size);
    g.CompleteWith(v, g_size);
  }

  return complexity;
}

int main(int argc, const char *argv[]) {
  std::vector<std::tuple<unsigned int, Word, Word>> complexity;
  std::map<unsigned int, size_t> complexity_count;
  size_t total = 0;

  std::ifstream words_in("h2_14_unproc_words.txt");

  std::pair<std::string, std::string> next_pair;

  words_in.ignore(64, '\n');

  auto PrintStats = [&total, &complexity_count](std::ostream& out) {
    out << "Total: " << total << "\n";
    for (auto&& count : complexity_count) {
      out << count.first << ": " << count.second << "\n";
    }
    out << std::flush;
  };

  auto SaveStats = [&complexity, &PrintStats]() {
    std::ofstream words_out("h2_14_pairs_complexity.txt");
    PrintStats(words_out);
    words_out << "\n";

    std::sort(complexity.begin(), complexity.end(),
      [](const std::tuple<unsigned int, Word, Word>& lhs, const std::tuple<unsigned int, Word, Word>& rhs) {

        if (std::get<0>(lhs) != std::get<0>(rhs)) {
          return std::get<0>(lhs) < std::get<0>(rhs);
        }
        return std::get<1>(lhs).size() + std::get<2>(lhs).size() < std::get<1>(rhs).size() + std::get<2>(rhs).size();
      });

    for (auto&& pair : complexity) {
      words_out << std::get<0>(pair) << ", " << std::get<1>(pair) << ", " << std::get<2>(pair) << "\n";
    }
  };

  while (words_in) {
    std::string appeared_in;
    words_in >> appeared_in;
    words_in >> next_pair.first;
    next_pair.first.pop_back();
    words_in >> next_pair.second;

    if (words_in) {
      ++total;
      Word u(next_pair.first);
      Word v(next_pair.second);
      auto pair_complexity = PairComplexity(u, v);
      ++complexity_count[pair_complexity];
      complexity.emplace_back(pair_complexity, u, v);

      if (total % 100 == 0) {
        PrintStats(std::cout);
      }

      if (total % 5000 == 0) {
        SaveStats();
      }
    }
  }

  PrintStats(std::cout);
  SaveStats();


  return 0;
}

