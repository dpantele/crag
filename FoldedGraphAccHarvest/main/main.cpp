#include "acc.h"

#include <bitset>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <sstream>
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
  uint64_t iterations_count_ = 0u;

  typedef std::chrono::microseconds ResultUnits;

  struct Iteration {
    Stopwatch* stopwatch_;
    clock::duration duration_{0};
    clock::time_point last_click_;
    uint64_t clicks_count_{0};

    Iteration(Stopwatch* stopwatch)
      : stopwatch_(stopwatch)
    { }

    inline bool IsRunning() const {
      return clicks_count_ % 2 == 1;
    }

    void Click() {
      if (IsRunning()) {
        duration_ += (clock::now() - last_click_);
      } else {
        last_click_ = clock::now();
      }

      ++clicks_count_;
    }

    void Reset() {
      if (IsRunning()) {
        Click();
      }
      if (clicks_count_) {
        stopwatch_->Report(this);
        clicks_count_ = 0;
        duration_ = clock::duration(0);
      }
    }

    ~Iteration() {
      Reset();
    }
  };

  void Report(Iteration* iteration) {
    ++iterations_count_;
    duration_ += iteration->duration_;
    last_duration_ = iteration->duration_;
  }

  Iteration NewIter() {
    return Iteration(this);
  }

  long long last() {
    return std::chrono::duration_cast<ResultUnits>(last_duration_).count();
  }

  long long average() {
    return iterations_count_ ? std::chrono::duration_cast<ResultUnits>(duration_ / iterations_count_).count() : 0;
  }

  uint64_t iterations() const {
    return iterations_count_;
  }

  long long total_time() const {
    return duration_.count();
  }
};

template<typename U, typename V>
std::pair<V, U> Swapped(std::pair<U, V> p) {
  return std::pair<V, U>(std::move(p.second), std::move(p.first));
}

//structures and algoritms for main below:

struct Stopwatches {
  typedef Stopwatch Field;
  std::array<Field, 4> data_;

  Field& folding() {
    return data_[0];
  }
  Field& harvest() {
    return data_[1];
  }
  Field& normalize() {
    return data_[2];
  }
  Field& reweight() {
    return data_[3];
  }

  long long total() const {
    long long total = 0;
    for (auto&& timer : data_) {
      total += timer.total_time();
    }
    return total;
  }
};


struct StopwatchesIteration {
  typedef Stopwatch::Iteration Field;
  std::array<Field, 4> data_;

  Field& folding() {
    return data_[0];
  }
  Field& harvest() {
    return data_[1];
  }
  Field& normalize() {
    return data_[2];
  }
  Field& reweight() {
    return data_[3];
  }

  StopwatchesIteration(Stopwatches* s)
    : data_({s->data_[0].NewIter(), s->data_[1].NewIter(), s->data_[2].NewIter(), s->data_[3].NewIter()})
  { }

  void Reset() {
    for (auto&& timer : data_) {
      timer.Reset();
    }
  }
};

struct Result {
  Word u_;
  Word v_;

  StopwatchesIteration time_;

  std::vector<std::pair<Word, Word>> new_pairs_;

  std::ostringstream estats_out_;

  Result(Word u, Word v, Stopwatches* s)
    : u_(std::move(u))
    , v_(std::move(v))
    , time_(s)
  { }
};

struct Parameters {
  size_t max_harvest_length = Word::kMaxLength;
  std::vector<uint16_t> complete_count;
  std::pair<std::string, std::string> initial_strings;
  std::pair<std::string, std::string> required_strings;

  Parameters()
    : complete_count(Word::kMaxLength + 1, 2)
    , initial_strings("xyxYXY", "xxxYYYY")
    , required_strings("x", "y")
  { }

};

void ProceedPair(const Word& u, const Word& v, const Parameters& p, Result* result) {
  auto& time = result->time_;

  time.folding().Click();
  FoldedGraph2 g;
  g.PushCycle(u, g.root(), 1);

  for (auto i = 0u; i < p.complete_count[v.size()]; ++i) {
    g.CompleteWith(v);
  }
  time.folding().Click();

  result->estats_out_ << g.size() << ", ";
  result->estats_out_ << g.CountNontrivialEdges() << ", ";

  time.reweight().Click();
  g.Reweight();
  time.reweight().Click();

  result->estats_out_ << g.CountNontrivialEdges() << ", ";

  time.harvest().Click();
  auto eq_u = g.Harvest(p.max_harvest_length, g.root());
  time.harvest().Click();

  result->estats_out_ << eq_u.size() << ", ";

  for (auto u_p = eq_u.begin(); u_p != eq_u.end(); ++u_p) {
    time.normalize().Click();
    auto new_pair = GetCanonicalPair(v, *u_p);
    time.normalize().Click();
    result->new_pairs_.push_back(std::move(new_pair));  
  }

  std::sort(result->new_pairs_.begin(), result->new_pairs_.end());
  auto new_end = std::unique(result->new_pairs_.begin(), result->new_pairs_.end());
  result->new_pairs_.erase(new_end, result->new_pairs_.end());
}

int main(int argc, const char *argv[]) {
  Parameters p;
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
      p.max_harvest_length = std::stoi(arg[1]);
    } else if (arg.front() == "comp") {
      auto up_to = std::stoi(arg[1]);
      auto count = std::stoi(arg[2]);
      if (p.complete_count.size() < up_to) {
        p.complete_count.resize(up_to, count);
      }
      for (size_t i = 0; i < up_to; ++i) {
        p.complete_count[i] = count; 
      }
    } else if (arg.front() == "init") {
      p.initial_strings.first = arg[1];
      p.initial_strings.second = arg[2];
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
  *out <<                 ", folding time";
  *out <<                 ", harvest time";
  *out << std::endl;

  if (estats_out.is_open()) {
    estats_out << "Configuration: " << std::endl;
    estats_out << "Max length for harvest: " << p.max_harvest_length << std::endl;
    estats_out << "Pairs are also normalized using automorphisms: yes";
    estats_out << "Initial words:          " << p.initial_strings.first << " | " << p.initial_strings.second << std::endl;
    estats_out << "\nHow many times graph is completed with v: " << std::endl;
    estats_out << "Length of v: ";
    for (auto i = 0u; i < p.complete_count.size(); ++i) {
      estats_out << std::setw(2) << i << ' ';
    }
    estats_out << "\nTimes:       ";
    for (auto i = 0u; i < p.complete_count.size(); ++i) {
      estats_out << std::setw(2) << p.complete_count[i] << ' ';
    }
    estats_out << "\n\n";
    estats_out << "iteration, unprocessed count, total count, ";
    estats_out << "graph size, ";
    estats_out << "non-zero edges before reweight, ";
    estats_out << "non-zero edges after reweight, ";
    estats_out << "harvest result size, ";
    estats_out << "after normalize size, ";
    estats_out << "u size, v size, ";
    estats_out << "folding time, average folding time, ";
    estats_out << "reweight time, average reweight time, ";
    estats_out << "harvest time, average harvest time, ";
    estats_out << "normalize time, average normalize time, ";
    estats_out << "\n";
  }

  auto initial = Swapped(GetCanonicalPair(p.initial_strings.first.c_str(), p.initial_strings.second.c_str()));
  auto required = Swapped(GetCanonicalPair(p.required_strings.first.c_str(), p.required_strings.second.c_str()));

  std::set<std::pair<Word, Word>> unprocessed_pairs = {};
  std::set<std::pair<Word, Word>> all_pairs = {};

  int counter = 0;

  auto AddPair = [&unprocessed_pairs, &all_pairs, &unproc_words, &counter](const std::pair<Word, Word>& pair) -> bool {
    auto exists = all_pairs.insert(pair);
    if (exists.second) {
      unprocessed_pairs.insert(pair);
      if (unproc_words.is_open()) {
        unproc_words << counter << ", ";
        PrintWord(pair.second, &unproc_words);
        unproc_words << ", ";
        PrintWord(pair.first, &unproc_words);
        unproc_words << "\n";
      }
      return true;
    }
    return false;
  };

  AddPair(initial);
  AddPair(Swapped(initial));

  Stopwatches total_time;
  Stopwatch overall_time;
  auto overall_time_stamp = overall_time.NewIter();
  overall_time_stamp.Click();

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

    Result result(u, v, &total_time);
    ProceedPair(u, v, p, &result);

    if (estats_out.is_open()) {
      estats_out << result.estats_out_.str();
    }

    auto pairs_count_before = all_pairs.size();

    for (auto&& new_pair : result.new_pairs_) {
      if (AddPair(Swapped(new_pair))) {
        AddPair(new_pair);
      }
    }

    if (estats_out.is_open()) {
      estats_out << all_pairs.size() - pairs_count_before << ", ";
      estats_out << u.size() << ", ";
      estats_out << v.size() << ", ";
    }

    result.time_.Reset();

    *out << total_time.folding().last() << ", ";
    *out << total_time.harvest().last() << std::endl;

    if (estats_out.is_open()) {
      estats_out << total_time.folding().last() << ", ";
      estats_out << total_time.folding().average() << ", ";
      estats_out << total_time.reweight().last() << ", ";
      estats_out << total_time.reweight().average() << ", ";
      estats_out << total_time.harvest().last() << ", ";
      estats_out << total_time.harvest().average() << ", ";
      estats_out << total_time.normalize().last() << ", ";
      estats_out << total_time.normalize().average() << "\n";
    }

    out->flush();
    estats_out.flush();
    proc_words.flush();
    unproc_words.flush();
  }
  overall_time_stamp.Click();
  overall_time_stamp.Reset();
  auto average_other_time = (overall_time.total_time() - total_time.total()) / total_time.folding().iterations();
  *out << "total time: " << overall_time.total_time() << std::endl;
  *out << "total rest(av): " << average_other_time << std::endl;
  *out << "total iterations: " << total_time.folding().iterations() << std::endl;


  return 0;
}
