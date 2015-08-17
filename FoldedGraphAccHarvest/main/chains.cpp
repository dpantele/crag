#include "acc.h"

#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <sstream>

using namespace crag;

using Label = FoldedGraph2::Label;

char LabelToChar(Label l) {
  return l % 2 ? 'X' + l / 2 : 'x' + l/2;
}

template<typename U, typename V>
std::pair<V, U> Swapped(std::pair<U, V> p) {
  return std::pair<V, U>(std::move(p.second), std::move(p.first));
}

int NormalFormExtendedBasic(Word* u) {
  auto u_min = *u;
  int min_shift = 0;

  int current_shift = 1;
  auto current_u = *u;
  current_u.CyclicLeftShift();

  while (current_u != *u) {
    if (current_u < u_min) {
      u_min = current_u;
      min_shift = current_shift;
    }

    current_u.CyclicLeftShift();
    ++current_shift;
  }


  *u = u_min;
  return min_shift;
}

std::tuple<int, bool, Word> NormalFormExtended(const Word& u) {
  auto u_min = u;
  auto min_shift = NormalFormExtendedBasic(&u_min);

  Word u_inv = u;
  u_inv.Invert();
  auto min_shift_inv = NormalFormExtendedBasic(&u_inv);

  if (u_inv < u_min) {
    return std::make_tuple(min_shift_inv, true, u_inv);
  } else {
    return std::make_tuple(min_shift, false, u_min);
  }
}

void PrintPair(const Word& u, const Word& v) {
  std::cout << "\\langle";
  PrintTo(u, &std::cout);
  std::cout << " \\mid ";
  PrintTo(v, &std::cout);
  std::cout << "\\rangle";
}

std::pair<Word, bool> GetCanConj(Word u) {
  auto normal_form = NormalFormExtended(u);
  if (std::get<1>(normal_form)) {
    u.Invert();
  }

  auto c = u;
  c.PopBack(c.size() - std::get<0>(normal_form));

  if (Conjugate(u, c) != std::get<2>(normal_form)) {
    std::cout << "\nu: ";
    PrintWord(u, &std::cout);
    std::cout << "\n Shift: " << std::get<0>(normal_form);
    std::cout << "\n New: ";
    PrintWord(Conjugate(u, c), &std::cout);
    std::cout << "\n Conj: ";
    PrintWord(c, &std::cout);
    std::cout << "\n Old: ";
    PrintWord(std::get<2>(normal_form), &std::cout);
    std::cout << std::endl;
    throw std::runtime_error("Incorrect shifting conjugate");
  }

  return std::make_pair(c, std::get<1>(normal_form));
}

std::pair<Word, Word> ExtendedGetCanoncial(Word u, Word v) {
  u = CyclicReduce(u);
  v = CyclicReduce(v);

  Word u_can, v_can;
  int u_shift, v_shift;
  bool u_inv, v_inv;

  std::tie(u_shift, u_inv, u_can) = NormalFormExtended(u);
  std::tie(v_shift, v_inv, v_can) = NormalFormExtended(v);

  auto old_can = GetCanonicalPair(u, v);

  if (old_can.first != u_can || old_can.second != v_can) {
    std::cout << " New: ";
    PrintPair(u_can, v_can);
    std::cout << " Old: ";
    PrintPair(old_can.first, old_can.second);
    std::cout << std::flush;
    throw std::runtime_error("New and old normal forms are not equal");
  }

  if (u_shift != 0) {
    std::cout << " u_{" << u_shift << "} ";
  }

  if (u_inv) {
    std::cout << " u_{inv} ";
  }

  if (v_shift != 0) {
    std::cout << " v_{" << v_shift << "} ";
  }

  if (v_inv) {
    std::cout << " v_{inv} " << " ";
  }

  return std::make_pair(u_can, v_can);
};

std::string ToString(const Word& w) {
  std::stringstream out;
  PrintWord(w, &out);
  return out.str();
}

int main(int argc, const char *argv[]) {
  std::string prefix = "h13_no_auto_ak4";
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
  checked_pair = std::make_pair(checked_pair.second, checked_pair.first);

  for (auto&& x_y_image : x_y_images) {
    Word x_image;
    Word y_image;

    std::tie(x_image, y_image) = x_y_image;

    std::cout << "\\begin{table}[h]\n"
                   "\\centering\n"
                   "\\begin{tabularx}{\\textwidth}{|X|X|X|X|}\n"
                   "\\hline\n"
                   "$u$ & $v$ & $c$ & $u'$\\\\\n"
                   "\\hline";

    std::array<Word, 4> images = {
      x_image,
      x_image,
      y_image,
      y_image,
    }; 

    images[1].Invert();
    images[3].Invert();

    auto TransformWord = [&images](Word w) -> Word {
      Word result;
      while (!w.Empty()) {
        result.PushBack(images[w.GetFront()]);
        w.PopFront();
      }
      return result;
    };

    auto u0 = CyclicReduce(TransformWord(checked_pair.first));
    auto u_c = GetCanConj(u0);

//    std::cout << "$(";
//    PrintWord(u_c.first, &std::cout);
//    std::cout << ")^{-1}";
//    std::cout << "(";
//    PrintWord(u0, &std::cout);
//    std::cout << ")";
    if (u_c.second) {
//      std::cout << "^{-1}";
      u0.Invert();
    }
//    std::cout << "(";
//    PrintWord(u_c.first, &std::cout);
//    std::cout << ") = ";
    u0 = Conjugate(u0, u_c.first);
//    PrintWord(u0, &std::cout);
//    std::cout << "$, \n";


    auto v0 = CyclicReduce(TransformWord(checked_pair.second));
    auto v_c = GetCanConj(v0);

//    std::cout << "$(";
//    PrintWord(v_c.first, &std::cout);
//    std::cout << ")^{-1}";
//    std::cout << "(";
//    PrintWord(v0, &std::cout);
//    std::cout << ")";
    if (v_c.second) {
//      std::cout << "^{-1}";
      v0.Invert();
    }
//    std::cout << "(";
//    PrintWord(v_c.first, &std::cout);
//    std::cout << ") = ";
    v0 = Conjugate(v0, v_c.first);
//    PrintWord(v0, &std::cout);
//    std::cout << "$, \n";

    std::cout << "$"<< ToString(u0) << "$ & $" << ToString(v0) << "$ & & \\\\\n";

    auto current_words = std::make_pair(u0, v0);

    if (current_words != GetCanonicalPair(TransformWord(checked_pair.first), TransformWord(checked_pair.second))) {
      throw std::runtime_error("Error");
    }

    auto next = generated_words.find(current_words);
    if (next == generated_words.end()) {
      std::cout << " missing" << std::endl;
      continue;
    }
//    std::cout << "\\begin{align*}\n";
    while (next->second != 0) {
//      std::cout << "\\,&\\langle";
//      PrintTo(current_words.first, &std::cout);
//      std::cout << " \\mid ";
//      PrintTo(current_words.second, &std::cout);
//      std::cout << "\\rangle";
      auto next_index = next->second;
      auto previous_words = current_words;
      current_words = generating_words[next_index];
      next = generated_words.find(current_words);

//      if (Swapped(current_words) == previous_words) {
//        std::cout << "\\;&\\text{swapped}";
//      } else {
      if (Swapped(current_words) != previous_words) {
        std::cout << "$"
                  << ToString(current_words.first)
                  << "$ & $"
                  << ToString(current_words.second)
                  << "$ & ";
        FoldedGraph2 g;
        g.PushCycle(current_words.first, g.root(), 1);
        for (auto i = 0u; i < 3; ++i) {
          g.CompleteWith(current_words.second);
        }

        for (auto v = g.root(); v < g.size(); ++v) {
          auto v_harvest = g.Harvest(previous_words.second.size(), v, v, 1);
          for (auto u : v_harvest) {
            if (GetCanonicalPair(current_words.second, u) == previous_words) {

              auto conjugator = g.GetPathFromRoot(v);
              auto path = g.ReadWord(u, v);
              if (std::get<0>(path) != v) {
                throw std::runtime_error("Word was not read");
              }

              if (std::get<1>(path) != u.size()) {
                throw std::runtime_error("Word was not read");
              }

              if (std::get<2>(path) != 1 && std::get<2>(path) != -1) {
                throw std::runtime_error("Weight is not identity");
              }

              auto cycl_conj = GetCanConj(u);
              //conjugator.PushBack(cycl_conj.first);

//              if (conjugator.Empty() && cycl_conj.first.Empty()) {
//                std::cout << "\\;&\\text{no conjugator}";
//              } else {
//                std::cout << "\\;&c = ";
                std::cout << "$";
                PrintWord(conjugator, &std::cout);
                if (conjugator.Empty()) {
                  std::cout << "\\varepsilon";
                }
                std::cout << " \\cdot ";
                PrintWord(cycl_conj.first, &std::cout);
                if (cycl_conj.first.Empty()) {
                  std::cout << "\\varepsilon";
                }
              std::cout << "$ & ";

//              }
              v = g.size();

//              std::cout << ", u' = ";
              if (cycl_conj.second) {
                u.Invert();
              }

//              std::cout << std::endl;
//              PrintWord(u, &std::cout);
//              std::cout << std::endl;
//              PrintWord(cycl_conj.first, &std::cout);
              u = Conjugate(u, cycl_conj.first);

              //here I got the normal form
              if (u != previous_words.second) {
                std::cout << std::endl;
                PrintWord(u, &std::cout);
                std::cout << std::endl;
                PrintPair(previous_words.first, previous_words.second);
                std::cout << std::endl;
                throw std::runtime_error("Error");
              }

              if (std::get<2>(path) == -1) {
                u.Invert();
              }

              if (cycl_conj.second) {
                u.Invert();
              }

              std::cout << "$";
              PrintWord(u, &std::cout);
              std::cout << "$\\\\\n";
              break;
            }
          }
        }
      }
    }

    std::cout << "\\hline\n"
                   "\\end{tabularx}\n"
                   "\\captionsetup{width=\\textwidth, justification=centering}\n"
                   "\\caption{ACM-moves for ";

    if (x_image != Word("x")) {
      std::cout << "$x \\rightarrow ";
      PrintWord(x_image, &std::cout);
      std::cout << "$ ";
    }

    if (y_image != Word("y")) {
      std::cout << " $y \\rightarrow ";
      PrintWord(y_image, &std::cout);
      std::cout << "$";
    }

    std::cout << ",\n\\\\";

    std::cout << "$\\varphi(xyxYXY) = ";
    if (!u_c.first.Empty()) {
      std::cout << "\\left(";
      PrintWord(u_c.first, &std::cout);
      std::cout << "\\right)";
    }
    if (u_c.second) {
      std::cout << "\\left(";
    }
    PrintWord(u0, &std::cout);
    if (u_c.second) {
      std::cout << "\\right)^{-1}";
    }
    if (!u_c.first.Empty()) {
      std::cout << "\\left(";
      PrintWord(u_c.first, &std::cout);
      std::cout << "\\right)^{-1}";
    }
    std::cout << "$, \n";

    std::cout << "$\\varphi(xxxYYYY) = ";
    if (!v_c.first.Empty()) {
      std::cout << "\\left(";
      PrintWord(v_c.first, &std::cout);
      std::cout << "\\right)";
    }
    if (v_c.second) {
      std::cout << "\\left(";
    }
    PrintWord(v0, &std::cout);
    if (v_c.second) {
      std::cout << "\\right)^{-1}";
    }
    if (!v_c.first.Empty()) {
      std::cout << "\\left(";
      PrintWord(v_c.first, &std::cout);
      std::cout << "\\right)^{-1}";
    }
    std::cout << "$\\\\\n";

    std::cout << "}\n"
                 "\\end{table}";



  }
  return 0;
}
