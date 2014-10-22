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

using Label = FoldedGraph2::Label;

char LabelToChar(Label l) {
  return l % 2 ? 'X' + l / 2 : 'x' + l/2;
}

template<typename T, size_t kSetSize>
class RandomSetPicker {
public:
  RandomSetPicker()
    : random_kicker_(0, kSetSize - 1) {
    current_set_.reserve(kSetSize);
  }

  template<typename RandomEngine>
  void Push(T value, RandomEngine& random_engine) {
    ++elements_count_;
    if (elements_count_ <= kSetSize) {
      current_set_.push_back(std::move(value));
      return;
    }
    if (random_number_(random_engine) * elements_count_ <= kSetSize) {
      //we pick new element with the probability k/n
      //so kick one of the previous with probability 1/k
      std::swap(current_set_[random_kicker_(random_engine)], current_set_.back());
      current_set_.back() = std::move(value);
    }
  }

  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;

  iterator begin() {
    return current_set_.begin();
  }

  iterator end() {
    return current_set_.end();
  }

  const_iterator begin() const {
    return current_set_.begin();
  }

  const_iterator end() const {
    return current_set_.end();
  }
private:
  std::vector<T> current_set_;
  std::uniform_real_distribution<double> random_number_;
  std::uniform_int_distribution<size_t> random_kicker_;
  size_t elements_count_ = 0;

};

int main(int argc, const char *argv[]) {
  std::ifstream words_in("diff.txt");
  std::random_device rd;
  std::mt19937_64 random_engine(rd());

  RandomSetPicker<std::pair<Word, Word>, 20> words;

  std::pair<std::string, std::string> next_pair;

  words_in.ignore(64, '\n');

  while (words_in) {
    std::string appeared_in;
    words_in >> appeared_in;
    words_in >> next_pair.first;
    next_pair.first.pop_back();
    words_in >> next_pair.second;
    if (words_in) {
      words.Push(std::make_pair(Word(next_pair.first), Word(next_pair.second)), random_engine);
    }
  }

  for (auto&& pair : words) {
    std::cout << "PAIR(\"";
    PrintWord(pair.first, &std::cout);
    std::cout << "\", \"";
    PrintWord(pair.second, &std::cout);
    std::cout << "\"),\n";
  }

  return 0;
}

