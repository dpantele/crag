#include "acc.h"

#include <deque>
#include <iostream>
#include <list>

using namespace crag;

void PrintWord(const Word& w, std::ostream* out) {
  for (auto&& l : w) {
    (*out) << l;
  }
}

int main() {
  std::vector<Word> conjugators = {
    {},
    {0u},
    {1u},
    {2u},
    {3u},
    {0u, 0u},
    {0u, 1u},
    {0u, 2u},
    {0u, 3u},
    {1u, 0u},
    {1u, 1u},
    {1u, 2u},
    {1u, 3u},
    {2u, 0u},
    {2u, 1u},
    {2u, 2u},
    {2u, 3u},
    {3u, 0u},
    {3u, 1u},
    {3u, 2u},
    {3u, 3u},
  };

  std::pair<Word, Word> initial = {{0u, 2u, 0u, 3u, 1u, 3u}, {0u, 0u, 0u, 2u, 2u, 2u, 2u}};
  std::pair<Word, Word> required = {{0u}, {2u}};

  std::set<std::pair<Word, Word>> unprocessed_pairs = {initial};
  std::set<std::pair<Word, Word>> all_pairs = {initial};

  while (!all_pairs.count(required) && !unprocessed_pairs.empty()) {
    std::cout << "all pairs: " << all_pairs.size() << std::endl;
    std::cout << "unporcessed: " << unprocessed_pairs.size() << std::endl;

    Word u, v;
    auto next_pair = *unprocessed_pairs.begin();
    unprocessed_pairs.erase(unprocessed_pairs.begin());
    std::tie(u, v) = std::move(next_pair);
    std::cout << "u = ";
    PrintWord(u, &std::cout);
    std::cout << "\nv = ";
    PrintWord(v, &std::cout);
    std::cout << std::endl;
    auto exists = all_pairs.emplace(v, u);
    if (exists.second) {
      if (v.size() < u.size())
        std::cout << "Adding" << std::endl;
      unprocessed_pairs.emplace(*exists.first);
    }

    std::vector<Word> new_u;

    std::set<Word> upp_s;
    
    for (size_t i = 0; i < u.size(); ++i) {
      LeftShift(u.begin(), u.end());
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
      FoldedGraph2 g;
      auto end = g.PushWord(upp);
      g.CompleteWith(v);
      g.CompleteWith(v);
      auto eq_u = g.Harvest(12, g.root(), end);
      ReduceAndNormalize(eq_u);
      new_u.reserve(new_u.size() + eq_u.size());
      for(auto& up : eq_u) {
        new_u.push_back(std::move(up));
      }
    }
    std::sort(new_u.begin(), new_u.end());
    auto new_u_end = std::unique(new_u.begin(), new_u.end());
    auto u_min_size = ~0u;
    for (auto u_p = new_u.begin(); u_p != new_u_end; ++u_p) {
      auto exists = all_pairs.emplace(v, *u_p);
      if (exists.second) {
        if (u_min_size > u_p->size()) {
          u_min_size = u_p->size();
        }
        unprocessed_pairs.emplace(*exists.first);
      }
    }
    std::cout << "Min size: " << u_min_size << std::endl;
    if (u_min_size <= 4) {
      std::cout << "done" << std::endl;
      return 0;
    }
  }

  return 0;
}