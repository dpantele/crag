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
#include <unordered_map>

using namespace crag;

using Label = FoldedGraph2::Label;

char LabelToChar(Label l) {
  return l % 2 ? 'X' + l / 2 : 'x' + l/2;
}

//! Simple reduction of *u in <x, y | v>
bool reduce(Word v, Word* u) {
  Word::size_type common_part_length, common_u_begin, common_v_begin;
  std::tie(common_u_begin, common_v_begin, common_part_length) = LongestCommonSubwordCyclic(*u, v);
  if (common_part_length > v.size() / 2) {
    auto u_copy = *u;
    u_copy.CyclicLeftShift(common_u_begin); //u = as
    auto v_copy = v; 
    v_copy.CyclicLeftShift(common_v_begin); //v = at
    v_copy.PopFront(common_part_length);    //v = t
    v_copy.Invert();                        //v = T
    u_copy.PopFront(common_part_length);    //u = s
    v_copy.PushBack(u_copy);                //v = Ts
    *u = v_copy;
    return true;
  }
  return false;
}

int main(int argc, const char *argv[]) {
  std::ifstream normalized_words_in("h2_16_proc_words.txt");

  std::set<std::pair<Word, Word>> all_words;

  std::pair<std::string, std::string> next_pair;

  normalized_words_in.ignore(64, '\n');

  while (normalized_words_in) {
    std::string appeared_in;
    normalized_words_in >> appeared_in;
    normalized_words_in >> next_pair.first;
    next_pair.first.pop_back();
    normalized_words_in >> next_pair.second;
    if (normalized_words_in) {
      all_words.emplace(Word(next_pair.first), Word(next_pair.second));
    }
  }

  std::cout << "Loaded " << all_words.size() << " pairs" << std::endl;

  std::map<Word, std::set<Word>> conj_classes;
  int reduced_pairs_count = 0;

  for (auto&& pair : all_words) {
    Word u, v, v_inv;
    std::tie(u, v) = pair;
    v_inv = v;
    v_inv.Invert();
    while(reduce(v, &u) || reduce(v_inv, &u)) { }
    auto result = GetCanonicalPair(u, v);
    auto exists = conj_classes[result.second].emplace(result.first);
    if (exists.second) {
      ++reduced_pairs_count;
    }
  }


  std::cout << "Total: " << all_words.size() << std::endl;
  std::cout << "After reduction: " << reduced_pairs_count << std::endl;
  std::cout << "Different v: " << conj_classes.size() << std::endl;

  std::map<int, int> count_sizes;

  for (auto&& conj_class : conj_classes) {
    ++count_sizes[conj_class.second.size()];
  }

  for (auto&& size : count_sizes) {
    std::cout << size.first << ' ' << size.second << std::endl;
  }

  return 0;
}
