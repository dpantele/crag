#include "acc.h"

#include <bitset>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <iomanip>
#include <list>
#include <string>

using namespace crag;

void PrintWord(const Word& w, std::ostream* out) {
  PrintTo(w, out);
}

std::vector<std::string> Split(const std::string& str, char split = ':') {
  std::vector<std::string> result;
  decltype(str.find(split)) last_split = 0;
  while (last_split < str.size()) {
    auto new_split = str.find(split, last_split);
    result.push_back(str.substr(last_split, new_split - last_split));
    if (new_split == std::string::npos) {
      break;
    }
    last_split = new_split + 1;
  }
  return result;
}

struct Stopwatch {
  using clock = std::chrono::high_resolution_clock;
  clock::duration duration_;
  clock::duration last_duration_;
  clock::time_point last_click_;
  uint64_t clicks_count_ = 0u;

  void click() {
    if (clicks_count_ % 2) {
      last_duration_ = (clock::now() - last_click_);
      duration_ += last_duration_;
    } else {
      last_click_ = clock::now();
    }

    ++clicks_count_;
  }
};

std::ostream& operator<<(std::ostream& out, const Stopwatch& s) {
  return out
    << std::chrono::duration_cast<std::chrono::microseconds>(s.last_duration_).count() << "mcs ("
    << std::chrono::duration_cast<std::chrono::microseconds>(s.duration_ / (s.clicks_count_ / 2)).count() << " avms)";
}
 
int main(int argc, const char *argv[]) {
  size_t max_harvest_length = Word::kMaxLength;
  std::vector<uint16_t> complete_count(Word::kMaxLength + 1, 2);
  std::pair<Word, Word> initial = {Word("xyxYXY"), Word("xxxYYYY")};
  std::pair<Word, Word> required = {Word("x"), Word("y")};

  for (int argi = 1; argi < argc; ++argi) {
    auto arg = Split(argv[argi]);
    if (arg.empty()) {
      continue;
    }

    if (arg.front() == "maxhl") {
      max_harvest_length = std::stoi(arg[1]);
    } else if (arg.front() == "comp") {
      auto up_to = std::stoi(arg[1]);
      auto count = std::stoi(arg[2]);
      if (complete_count.size() < up_to) {
        complete_count.resize(up_to, count);
      }
      for (size_t i = 0; i < up_to; ++i) {
        complete_count[i] = count; 
      }
    } else if (arg.front() == "init") {
      initial.first = Word(arg[1]);
      initial.second = Word(arg[2]);
    }
  }

  auto normalized = ReduceAndNormalize({initial.first, initial.second});
  initial = {normalized[0], normalized[1]};

  std::cout << "Configuration: " << std::endl;
  std::cout << "Max length for harvest: " << max_harvest_length << std::endl;
  std::cout << "Initial words:          ";
  PrintWord(initial.first, &std::cout);
  std::cout << " | ";
  PrintWord(initial.second, &std::cout);
  std::cout << "\nHow many times graph is completed with v: " << std::endl;
  std::cout << "Length of v: ";
  for (auto i = 0u; i < complete_count.size(); ++i) {
    std::cout << std::setw(2) << i << ' ';
  }
  std::cout << "\nTimes:       ";
  for (auto i = 0u; i < complete_count.size(); ++i) {
    std::cout << std::setw(2) << complete_count[i] << ' ';
  }
  std::cout << "\n\n";

  std::set<std::pair<Word, Word>> unprocessed_pairs = {initial};
  std::set<std::pair<Word, Word>> all_pairs = {initial};

  int counter = 0;
  Stopwatch folding_time;
  Stopwatch harvest_time;
  Stopwatch normalize_time;

  while (!all_pairs.count(required) && !unprocessed_pairs.empty()) {
    ++counter;
    std::cout << std::left << std::setw(6) << counter << ", ";
    std::cout.width(9);
    std::cout << unprocessed_pairs.size() << ", ";
    std::cout.width(9);
    std::cout << all_pairs.size() << ", ";
    // std::cout << std::endl;

    Word u, v;
    auto next_pair = *unprocessed_pairs.begin();
    unprocessed_pairs.erase(unprocessed_pairs.begin());
    std::tie(v, u) = std::move(next_pair);
    //std::cout << "u = ";
    std::cout << "  (";
    PrintWord(u, &std::cout);
    std::cout << " | ";
    // std::cout << "\nv = ";
    PrintWord(v, &std::cout);
    std::cout << ")";
    auto exists = all_pairs.emplace(u, v);
    if (exists.second) {
      unprocessed_pairs.emplace(*exists.first);
    }

    folding_time.click();
    FoldedGraph2 g;
    g.PushCycle(u, g.root(), 1);

    for (auto i = 0u; i < complete_count[v.size()]; ++i) {
      g.CompleteWith(v);
    }
    folding_time.click();

    harvest_time.click();
    auto eq_u = g.Harvest(max_harvest_length, g.root());
    harvest_time.click();

    normalize_time.click();
    eq_u = ReduceAndNormalize(eq_u);
    normalize_time.click();

    std::cout << "\nfold " << folding_time;
    std::cout << "\nharvest " << harvest_time;
    std::cout << "\nnormalize " << normalize_time;

    std::bitset<Word::kMaxLength> available_sizes;
    for (auto u_p = eq_u.begin(); u_p != eq_u.end(); ++u_p) {
      auto exists = all_pairs.emplace(*u_p, v);
      if (exists.second) {
        if (u_p->size() > 0) {
          available_sizes.set(u_p->size() - 1);
        }
        unprocessed_pairs.emplace(*exists.first);
      }
    }
    std::cout << ",\tsz: ";
    for (auto sz = 0u; sz < Word::kMaxLength; ++sz) {
      if (available_sizes[sz]) {
        std::cout << sz + 1 << ",";
      }
    }
    std::cout << std::endl;
  }

  return 0;
}
