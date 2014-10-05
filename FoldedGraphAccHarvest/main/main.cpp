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
  clock::duration duration_{0};
  clock::duration last_duration_{0};
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

  typedef std::chrono::microseconds ResultUnits;

  long long last() {
    return std::chrono::duration_cast<ResultUnits>(last_duration_).count();
  }

  long long average() {
    return std::chrono::duration_cast<ResultUnits>(duration_ / (clicks_count_ / 2)).count();
  }
};

int main(int argc, const char *argv[]) {
  size_t max_harvest_length = Word::kMaxLength;
  std::vector<uint16_t> complete_count(Word::kMaxLength + 1, 2);
  std::pair<Word, Word> initial = {Word("xyxYXY"), Word("xxxYYYY")};
  std::pair<Word, Word> required = {Word("x"), Word("y")};

  std::vector<Word> ws = {Word("yxyXYX"), Word("yyyXXXX")};
  ws = ReduceAndNormalize(ws);
  PrintTo(ws[0], &std::cout);
  std::cout << std::endl;
  PrintTo(ws[1], &std::cout);
  std::cout << std::endl;

  std::ostream* out = &std::cout;

  std::ofstream stats_out;
  std::ofstream unproc_words;
  std::ofstream proc_words;
  std::ofstream estats_out;

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
    } else if (arg.front() == "out") {
      stats_out.open(arg[1] + ".txt");
      out = &stats_out;

      proc_words.open(arg[1] + "_proc_words.txt");
      proc_words << "iter, u, v\n";
      unproc_words.open(arg[1] + "_unproc_words.txt");
      unproc_words << "u, v\n";
      estats_out.open(arg[1] + "_ext.txt");
    }
  }

  *out << std::left ;
  *out << std::setw(7) << "iter";
  *out << std::setw(9) << ", unproc";
  *out << std::setw(9) << ", total";
  *out << std::setw(4) << ", u";
  *out << std::setw(4) << ", v";
  *out <<                 ", min new size";
  *out <<                 ", folding time";
  *out <<                 ", harvest time";
  *out << std::endl;

  auto normalized = ReduceAndNormalize({initial.first, initial.second});
  initial = {normalized[0], normalized[1]};

  if (estats_out.is_open()) {
    estats_out << "Configuration: " << std::endl;
    estats_out << "Max length for harvest: " << max_harvest_length << std::endl;
    estats_out << "Initial words:          ";
    PrintTo(initial.first, &estats_out);
    estats_out << " | ";
    PrintTo(initial.second, &estats_out);
    estats_out << "\nHow many times graph is completed with v: " << std::endl;
    estats_out << "Length of v: ";
    for (auto i = 0u; i < complete_count.size(); ++i) {
      estats_out << std::setw(2) << i << ' ';
    }
    estats_out << "\nTimes:       ";
    for (auto i = 0u; i < complete_count.size(); ++i) {
      estats_out << std::setw(2) << complete_count[i] << ' ';
    }
    estats_out << "\n\n";
    estats_out << "iteration, unprocessed count, total count, ";
    estats_out << "graph size, ";
    estats_out << "non-zero edges before reweight, ";
    estats_out << "non-zero edges after reweight, ";
    estats_out << "harvest result size, ";
    estats_out << "after normalize size, ";
    estats_out << "u size, v size, new sizes, ";
    estats_out << "folding time, average folding time, ";
    estats_out << "reweight time, average reweight time, ";
    estats_out << "harvest time, average harvest time, ";
    estats_out << "normalize time, average normalize time, ";
    estats_out << "\n";
  }

  std::set<std::pair<Word, Word>> unprocessed_pairs = {initial};
  std::set<std::pair<Word, Word>> all_pairs = {initial};

  int counter = 0;
  Stopwatch folding_time;
  Stopwatch harvest_time;
  Stopwatch normalize_time;
  Stopwatch reweight_time;

  while (!all_pairs.count(required) && !unprocessed_pairs.empty()) {
    ++counter;
    *out << std::left << std::setw(7) << counter << ", ";
    *out << std::right << std::setw(7) << unprocessed_pairs.size() << ", ";
    *out << std::right << std::setw(7) << all_pairs.size() << ", ";

    if (estats_out.is_open()) {
      estats_out << counter << ", ";
      estats_out << unprocessed_pairs.size() << ", ";
      estats_out << all_pairs.size() << ", ";
    }

    Word u, v;
    auto next_pair = *unprocessed_pairs.begin();
    unprocessed_pairs.erase(unprocessed_pairs.begin());
    std::tie(v, u) = std::move(next_pair);
    *out << std::setw(2) << u.size() << ", ";
    *out << std::setw(2) << v.size() << ", ";

    if (proc_words.is_open()) {
      proc_words << counter << ", ";
      PrintWord(u, &proc_words);
      proc_words << ", ";
      PrintWord(v, &proc_words);
      proc_words << "\n";
    }

    auto exists = all_pairs.emplace(u, v);
    if (exists.second) {
      unprocessed_pairs.emplace(*exists.first);

      if (unproc_words.is_open()) {
        PrintWord(v, &unproc_words);
        unproc_words << ", ";
        PrintWord(u, &unproc_words);
        unproc_words << "\n";
      }
    }

    folding_time.click();
    FoldedGraph2 g;
    g.PushCycle(u, g.root(), 1);

    for (auto i = 0u; i < complete_count[v.size()]; ++i) {
      g.CompleteWith(v);
    }
    folding_time.click();

    if (estats_out.is_open()) {
      estats_out << g.size() << ", ";
      estats_out << g.CountNontrivialEdges() << ", ";
    }

    reweight_time.click();
    g.Reweight();
    reweight_time.click();

    if (estats_out.is_open()) {
      estats_out << g.CountNontrivialEdges() << ", ";
    }

    harvest_time.click();
    auto eq_u = g.Harvest(max_harvest_length, g.root());
    harvest_time.click();

    if (estats_out.is_open()) {
      estats_out << eq_u.size() << ", ";
    }

    normalize_time.click();
    eq_u = ReduceAndNormalize(eq_u);
    normalize_time.click();

    if (estats_out.is_open()) {
      estats_out << eq_u.size() << ", ";
      estats_out << u.size() << ", ";
      estats_out << v.size() << ", ";
    }

    std::bitset<Word::kMaxLength> available_sizes;
    for (auto u_p = eq_u.begin(); u_p != eq_u.end(); ++u_p) {
      auto exists = all_pairs.emplace(*u_p, v);
      if (exists.second) {
        if (u_p->size() > 0) {
          available_sizes.set(u_p->size() - 1);
        }
        unprocessed_pairs.emplace(*exists.first);
        if (unproc_words.is_open()) {
          PrintWord(v, &unproc_words);
          unproc_words << ", ";
          PrintWord(*u_p, &unproc_words);
          unproc_words << "\n";
        }
      }
    }

    bool is_first = true;
    for (auto sz = 0u; sz < Word::kMaxLength; ++sz) {
      if (available_sizes[sz]) {
        if (is_first) {
          *out << sz + 1;
        }
        if (estats_out.is_open()) {
          estats_out << sz + 1 << "; ";
        }
        is_first = false;
      }
    }

    *out << ", " << folding_time.last() << ", ";
    *out << harvest_time.last() << std::endl;

    if (estats_out.is_open()) {
      estats_out << ", ";
      estats_out << folding_time.last() << ", ";
      estats_out << folding_time.average() << ", ";
      estats_out << reweight_time.last() << ", ";
      estats_out << reweight_time.average() << ", ";
      estats_out << harvest_time.last() << ", ";
      estats_out << harvest_time.average() << ", ";
      estats_out << normalize_time.last() << ", ";
      estats_out << normalize_time.average() << "\n";
    }

    out->flush();
    estats_out.flush();
    proc_words.flush();
    unproc_words.flush();
  }

  return 0;
}
