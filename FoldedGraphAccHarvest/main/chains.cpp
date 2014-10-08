#include "acc.h"

#include <fstream>
#include <map>
#include <string>
#include <utility>

using namespace crag;

using Label = FoldedGraph2::Label;

char LabelToChar(Label l) {
  return l % 2 ? 'X' + l / 2 : 'x' + l/2;
}

int main(int argc, const char *argv[]) {
  std::string prefix = "h14_c3";
  std::ifstream generated_words_in(prefix + "_unproc_words.txt");

  std::map<std::pair<Word, Word>, unsigned int> generated_words;

  std::pair<std::string, std::string> next_pair;

  generated_words_in.ignore(64, '\n');
  while (generated_words_in) {
    unsigned int appeared_in;
    generated_words_in >> appeared_in;
    generated_words_in.ignore(8, ' ');
    generated_words_in >> next_pair.first;
    next_pair.first.pop_back();
    generated_words_in >> next_pair.second;
    if (generated_words_in) {
      generated_words.emplace(
        std::make_pair(Word(next_pair.first), Word(next_pair.second)), 
        appeared_in
      );
    }
  }

  std::vector<std::pair<Word, Word>> generating_words;
  generating_words.emplace_back(Word(), Word());

  std::ifstream generating_words_in(prefix + "_proc_words.txt");
  generating_words_in.ignore(64, '\n');

  while (generating_words_in) {
    unsigned int appeared_in;
    generating_words_in >> appeared_in;
    generating_words_in.ignore(8, ' ');
    generating_words_in >> next_pair.first;
    next_pair.first.pop_back();
    generating_words_in >> next_pair.second;
    if (generating_words_in) {
      generating_words.emplace_back(
        std::make_pair(Word(next_pair.first), Word(next_pair.second))
      );
    }
  }


  //x->X, y->Y, x<->y, x->xy, x->yx, y->yx, y->xy
  std::pair<Word, Word> x_y_images[] = {
    {Word("X"), Word("y")},
    {Word("x"), Word("Y")},
    {Word("y"), Word("x")},
    {Word("xy"), Word("y")},
    {Word("yx"), Word("y")},
    {Word("x"), Word("yx")},
    {Word("x"), Word("xy")},
  };

  auto checked_pair = generating_words[1];//GetCanonicalPair("xyXYxYxyyXy", "xyxYxYXyXy");

  for (auto&& x_y_image : x_y_images) {
    Word x_image;
    Word y_image;

    std::tie(x_image, y_image) = x_y_image;

    std::array<Word, 4> images = {
      x_image,
      x_image,
      y_image,
      y_image,
    }; 

    images[1].Invert();
    images[3].Invert();

    for (Label l = 0; l < 4; ++l) {
      std::cout << "$" << LabelToChar(l) << " \\rightarrow ";
      PrintWord(images[l], &std::cout);
      std::cout << "$, ";
    }
    std::cout << "\n\n";

    auto TransformWord = [&images](Word w) -> Word {
      Word result;
      while (!w.Empty()) {
        result.PushBack(images[w.GetFront()]);
        w.PopFront();
      }
      return result;
    };

    auto TransformPair = [&TransformWord](const std::pair<Word, Word>& words) {
      return GetCanonicalPair(TransformWord(words.first), TransformWord(words.second));
    };

    auto current_words = TransformPair(checked_pair);

    auto next = generated_words.find(current_words);
    if (next == generated_words.end()) {
      std::cout << " missing" << std::endl;
      continue;
    }
    std::cout << "\\[\\begin{split}\n";
    while (next->second != 0) {
      std::cout << "\\langle";
      PrintTo(current_words.first, &std::cout);
      std::cout << " \\mid ";
      PrintTo(current_words.second, &std::cout);
      std::cout << "\\rangle\\\\ \n&\\leftarrow ";
      auto next_index = next->second;
      current_words = generating_words[next_index];
      next = generated_words.find(current_words);
    }
    std::cout << "\\langle";
    PrintTo(current_words.first, &std::cout);
    std::cout << " \\mid ";
    PrintTo(current_words.second, &std::cout);
    std::cout << "\\rangle\n";
    std::cout << "\\end{split}\\]\n\n";
  }
  return 0;
}
