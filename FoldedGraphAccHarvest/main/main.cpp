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

    ~Iteration() {
      if (IsRunning()) {
        Click();
      }
      stopwatch_->Report(this);
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
};

template<typename U, typename V>
std::pair<V, U> Swapped(std::pair<U, V> p) {
  return std::pair<V, U>(std::move(p.second), std::move(p.first));
}

Word Inverse(Word u) {
  u.Invert();
  return u;
}

#define MORPHISM(X, Y) {Word(X), Inverse(Word(X)), Word(Y), Inverse(Word(Y))}
static const Mapping length_changing_autos[] = {
  MORPHISM("yx", "y"),  //8
  MORPHISM("Yx", "y"),  //9
  MORPHISM("xy", "y"),  //10
  MORPHISM("xY", "y"),  //11
  MORPHISM("yxY", "y"), //12
  MORPHISM("Yxy", "y"), //13
  MORPHISM("x", "yx"),  //14
  MORPHISM("x", "yX"),  //15
  MORPHISM("x", "xy"),  //16
  MORPHISM("x", "Xy"),  //17
  MORPHISM("x", "Xyx"), //18
  MORPHISM("x", "xyX"), //19
};
#undef MORPHISM


int main(int argc, const char *argv[]) {
  size_t max_harvest_length = Word::kMaxLength;
  std::vector<uint16_t> complete_count(Word::kMaxLength + 1, 2);
  auto initial_strings = std::pair<std::string, std::string>("xyxYXY", "xxxYYYY");
  auto required_strings = std::pair<std::string, std::string>("x", "y");

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
      initial_strings.first = arg[1];
      initial_strings.second = arg[2];
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

  if (estats_out.is_open()) {
    estats_out << "Configuration: " << std::endl;
    estats_out << "Max length for harvest: " << max_harvest_length << std::endl;
    estats_out << "Pairs are also normalized using automorphisms: yes";
    estats_out << "Initial words:          " << initial_strings.first << " | " << initial_strings.second << std::endl;
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

  auto initial = GetCanonicalPair(initial_strings.first.c_str(), initial_strings.second.c_str());
  auto required = GetCanonicalPair(required_strings.first.c_str(), required_strings.second.c_str());

  std::set<std::pair<Word, Word>> unprocessed_pairs = {};
  std::set<std::pair<Word, Word>> all_pairs = {};

  int counter = 0;

  std::vector<std::set<std::pair<Word, Word>>::const_iterator> unprocessed_elements;
  std::vector<std::set<std::pair<Word, Word>>::const_iterator> new_unprocessed_elements;

  //wrong naming - actually it just used as something what adds pairs to 
  auto AddPairToAll = [&all_pairs, &new_unprocessed_elements](const std::pair<Word, Word>& pair) -> bool {
    auto insert_result = all_pairs.emplace(pair.first, pair.second);
    if (insert_result.second) {
      new_unprocessed_elements.push_back(insert_result.first);
      return true;
    }
    return false;
  };

  auto AddPairAsUnprocessed = [&AddPairToAll, &unproc_words, &unprocessed_pairs, &counter](const std::pair<Word, Word>& pair) {
    unprocessed_pairs.insert(pair);
    if (unproc_words.is_open()) {
      unproc_words << counter << ", ";
      PrintWord(pair.second, &unproc_words);
      unproc_words << ", ";
      PrintWord(pair.first, &unproc_words);
      unproc_words << "\n";
    }
  };

  auto AddPair = [
    &unprocessed_elements,
    &new_unprocessed_elements, 
    &AddPairAsUnprocessed, 
    &AddPairToAll,
    max_harvest_length
  ](const std::pair<Word, Word>& pair) -> bool {
    assert(pair == GetCanonicalPair(pair.first, pair.second));
    unprocessed_elements.clear();
    new_unprocessed_elements.clear();
    if (AddPairToAll(pair)) {
      AddPairAsUnprocessed(pair);
      AddPairToAll(Swapped(pair));
      AddPairAsUnprocessed(Swapped(pair));
      
      while (!new_unprocessed_elements.empty()) {
        std::swap(unprocessed_elements, new_unprocessed_elements);
        new_unprocessed_elements.clear();
        for (auto&& just_added_elem : unprocessed_elements) {
          for (auto&& morhpism : length_changing_autos) {
            try {
              auto images = std::make_pair(
                CyclicReduce(Map(just_added_elem->first, morhpism)),
                CyclicReduce(Map(just_added_elem->second, morhpism))
              );

              if (images.first.size() > max_harvest_length || images.second.size() > max_harvest_length) {
                continue;
              }

              images = GetOrbitCanonicalPair(images.first, images.second);

              if(AddPairToAll(images)) {
                AddPairToAll(Swapped(images));
              }
            } catch(std::length_error&) {
              //it's ok
            }
          }
        }
      }
      return true;
    }
    return false;
  };

  AddPair(initial);

  Stopwatch folding_time_total;
  Stopwatch harvest_time_total;
  Stopwatch normalize_time_total;
  Stopwatch reweight_time_total;

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

    //auto exists = all_pairs.insert(Swapped(GetCanonicalPair(v, u)));
    //if (exists.second) {
    //  unprocessed_pairs.emplace(*exists.first);

    //  if (unproc_words.is_open()) {
    //    unproc_words << counter << ", ";
    //    PrintWord(exists.first->second, &unproc_words);
    //    unproc_words << ", ";
    //    PrintWord(exists.first->first, &unproc_words);
    //    unproc_words << "\n";
    //  }
    //}

    auto folding_time = folding_time_total.NewIter();
    folding_time.Click();
    FoldedGraph2 g;
    g.PushCycle(u, g.root(), 1);

    for (auto i = 0u; i < complete_count[v.size()]; ++i) {
      g.CompleteWith(v);
    }
    folding_time.Click();

    if (estats_out.is_open()) {
      estats_out << g.size() << ", ";
      estats_out << g.CountNontrivialEdges() << ", ";
    }

    auto reweight_time = reweight_time_total.NewIter();
    reweight_time.Click();
    g.Reweight();
    reweight_time.Click();

    if (estats_out.is_open()) {
      estats_out << g.CountNontrivialEdges() << ", ";
    }

    auto harvest_time = harvest_time_total.NewIter();
    harvest_time.Click();
    auto eq_u = g.Harvest(max_harvest_length, g.root());
    harvest_time.Click();

    if (estats_out.is_open()) {
      estats_out << eq_u.size() << ", ";
    }

    auto normalize_time = normalize_time_total.NewIter();

    auto pairs_count_before = all_pairs.size();

    std::bitset<Word::kMaxLength> available_sizes;
    for (auto u_p = eq_u.begin(); u_p != eq_u.end(); ++u_p) {

      normalize_time.Click();
      auto new_pair = GetOrbitCanonicalPair(v, *u_p);
      normalize_time.Click();

      if (all_pairs.count(new_pair) != 0) {
        continue;
      }

      normalize_time.Click();
      auto new_pair = GetCanonicalPair(new_pair.first, new_pair.second);
      normalize_time.Click();
  
      if (AddPair(new_pair)) {
        if (u_p->size() > 0) {
          available_sizes.set(u_p->size() - 1);
        }
      }
    }

    if (estats_out.is_open()) {
      estats_out << all_pairs.size() - pairs_count_before << ", ";
      estats_out << u.size() << ", ";
      estats_out << v.size() << ", ";
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

    *out << ", " << folding_time_total.last() << ", ";
    *out << harvest_time_total.last() << std::endl;

    if (estats_out.is_open()) {
      estats_out << ", ";
      estats_out << folding_time_total.last() << ", ";
      estats_out << folding_time_total.average() << ", ";
      estats_out << reweight_time_total.last() << ", ";
      estats_out << reweight_time_total.average() << ", ";
      estats_out << harvest_time_total.last() << ", ";
      estats_out << harvest_time_total.average() << ", ";
      estats_out << normalize_time_total.last() << ", ";
      estats_out << normalize_time_total.average() << "\n";
    }

    out->flush();
    estats_out.flush();
    proc_words.flush();
    unproc_words.flush();
  }

  return 0;
}
