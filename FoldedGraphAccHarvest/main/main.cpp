#include "acc.h"

#include <bitset>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <list>

using namespace crag;

void PrintWord(const Word& w, std::ostream* out) {
  PrintTo(w, out);
}

int main(int argc, char *argv[]) {
  std::vector<Word> conjugators = {
    {},
    {0u},
    {1u},
    {2u},
    {3u},
    {0u, 0u},
    {0u, 2u},
    {0u, 3u},
    {1u, 1u},
    {1u, 2u},
    {1u, 3u},
    {2u, 0u},
    {2u, 1u},
    {2u, 2u},
    {3u, 0u},
    {3u, 1u},
    {3u, 3u},
  };

  std::pair<Word, Word> initial = { { 0u, 2u, 0u, 3u, 1u, 3u }, { 0u, 0u, 0u, 3u, 3u, 3u, 3u } };
  // std::pair<Word, Word> initial = { { 0u, 2u, 0u, 3u, 1u, 3u }, { 0u, 0u, 3u, 3u, 3u } };
  std::pair<Word, Word> required = { { 0u }, { 2u } };

  auto normalized = ReduceAndNormalize({initial.first, initial.second});
  initial = {normalized[0], normalized[1]};

  std::set<std::pair<Word, Word>> unprocessed_pairs = {initial};
  std::set<std::pair<Word, Word>> all_pairs = {initial};

  int counter = 0;

  size_t max_harvest_length = 16;
  if (argc > 0) {
    max_harvest_length = std::strtol(argv[1], nullptr, 10);
  }
  while (!all_pairs.count(required) && !unprocessed_pairs.empty()) {
    ++counter;
    std::cout << counter << ": ";
    std::cout.width(9);
    std::cout << unprocessed_pairs.size();
    std::cout.width(9);
    std::cout << all_pairs.size();
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

    std::vector<Word> new_u;

    std::set<Word> upp_s;
    
    for (size_t i = 0; i < u.size(); ++i) {
      u.CyclicLeftShift();
      for(const auto& s : conjugators) {
        //std::cout << "s=";
        //PrintWord(s, &std::cout);
        upp_s.insert(Conjugate(u, s));
        //std::cout << "\nupp=";
        //PrintWord(upp, &std::cout);
        //std::cout << std::endl;
      }
    }

    for (auto&& upp : upp_s) {
      if (upp.Empty()) {
        continue;
      }
      // PrintWord(upp, &std::cout);
      // std::cout << std::endl;

      FoldedGraph2 g;
      auto end = g.PushWord(upp);
      g.CompleteWith(v);
      g.CompleteWith(v);
      if (v.size() <= 7) g.CompleteWith(v);
      if (v.size() <= 5) g.CompleteWith(v);
      auto eq_u = g.Harvest(max_harvest_length, g.root(), end);
      eq_u = ReduceAndNormalize(eq_u);
      new_u.reserve(new_u.size() + eq_u.size());
      for(auto& up : eq_u) {
        new_u.push_back(std::move(up));
        // PrintWord(up, &std::cout);
        // std::cout << std::endl;
      }
    }

    std::sort(new_u.begin(), new_u.end());
    auto new_u_end = std::unique(new_u.begin(), new_u.end());
    std::bitset<Word::kMaxLength> available_sizes;
    for (auto u_p = new_u.begin(); u_p != new_u_end; ++u_p) {
      auto exists = all_pairs.emplace(*u_p, v);
      if (exists.second) {
        if (u_p->size() > 0) {
          available_sizes.set(u_p->size() - 1);
        }
        unprocessed_pairs.emplace(*exists.first);
      }
    }
    std::cout << "\tsz: ";
    for (auto sz = 0u; sz < Word::kMaxLength; ++sz) {
      if (available_sizes[sz]) {
        std::cout << sz + 1 << ",";
      }
    }
    std::cout << std::endl;
  }

  return 0;
}
