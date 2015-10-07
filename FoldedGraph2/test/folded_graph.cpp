#include "gtest/gtest.h"
#include "folded_graph2.h"

#include <algorithm>
#include <chrono>
#include <deque>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <memory>
#include <random>
#include <utility>

#include "acc.h"

namespace crag {

namespace {

using Label = FoldedGraph2::Label;
using Word = FoldedGraph2::Word;
using Weight = FoldedGraph2::Weight;
using Vertex = FoldedGraph2::Vertex;

namespace naive_graph_folding {
  

Weight Gcd(Weight a, Weight b) {
  while (b != 0) {
    auto c = a % b;
    a = b;
    b = c;
  }

  return a;
}

Label Inverse(Label l) {
  return l ^ 1;
}

struct Graph {
  struct Edge {
    std::tuple<Label, Vertex, Weight> data_;

    Label label() const {
      return std::get<0>(data_);
    }

    Vertex to() const {
      return std::get<1>(data_);
    }

    Weight w() const {
      return std::get<2>(data_);
    }

    Edge(Label label, Vertex to, Weight w)
      : data_(label, to, w)
    { }

    bool operator<(const Edge& other) const {
      return data_ < other.data_;
    }

    bool operator==(const Edge& other) const {
      return data_ == other.data_;
    }

    std::tuple<Label, Vertex> DropWeight() const {
      return std::make_tuple(label(), to());
    }
  };

  std::vector<std::set<Edge>> edges_;

  static Vertex root() {
    return 1;
  }

  Weight modulus_ = 0;

  Weight modulus() const {
    return modulus_;
  }

  Graph() {
    edges_.emplace_back();
    edges_.emplace_back();
  }

  void AddEdge(Vertex v1, Vertex v2, Label l, Weight w) {
    edges_.at(v1).emplace(l, v2, WeightMod(w));
    edges_.at(v2).emplace(Inverse(l), v1, WeightMod(-w));
  }

  Vertex PushWord(Word w, Vertex s = 1, Weight weight = 0) {
    while (!w.Empty()) {
      Vertex next_id = static_cast<Vertex>(edges_.size());
      edges_.emplace_back();
      AddEdge(s, next_id, w.GetFront(), (w.size() == 1 ? weight : 0));
      s = next_id;
      w.PopFront();
    }

    return s;
  }

  void PushCycle(Word w, Vertex s = 1, Weight weight = 0) {
    if (w.Empty()) {
      modulus_ = Gcd(weight, modulus_);
    }

    auto last = w.GetBack();
    w.PopBack();

    Vertex n = PushWord(w, s, 0);
    AddEdge(n, s, last, weight);
  }

  std::tuple<Vertex, Word::size_type, Weight> ReadWord(Word w, Vertex s = 1) {
    Word::size_type read_length = 0;
    Weight current_weight = 0;
    while (!w.Empty()) {
      auto edge = edges_.at(s).lower_bound(Edge(w.GetFront(), 0, 0));
      if (edge == edges_.at(s).end() || edge->label() != w.GetFront()) {
        return std::make_tuple(s, read_length, current_weight);
      }
      auto next = std::next(edge);
      if (next != edges_.at(s).end() && next->label() == w.GetFront()) {
        throw std::runtime_error("Graph is not foled - can't traverse unfolded graphs");
      }

      s = edge->to();
      ++read_length;
      current_weight += edge->w();
      w.PopFront();
    }
    return std::make_tuple(s, read_length, WeightMod(current_weight));
  }

  void ShiftWeight(Vertex v, Weight shift) {
    std::set<Edge> reweightened_edges_for_v;
    //first, remove old edges
    for (auto&& edge : this->edges_[v]) {
      if (edge.to() != v) {
        auto opposite_count = edges_.at(edge.to()).erase(Edge(Inverse(edge.label()), v, WeightMod(-edge.w())));
        assert(opposite_count == 1);
      }
    }
    //now add new edges
    for (auto&& edge : this->edges_[v]) {
      if (edge.to() != v) {
        reweightened_edges_for_v.emplace(edge.label(), edge.to(), WeightMod(edge.w() + shift));
        edges_.at(edge.to()).emplace(Inverse(edge.label()), v, WeightMod(-edge.w() - shift));
      } else {
        reweightened_edges_for_v.emplace(edge);
      }
    }
    edges_[v] = std::move(reweightened_edges_for_v);
  }

  bool JoinVertices(Vertex v1, Vertex v2) {
    static int count = 0;
    ++count;
    if (v1 > v2) {
      std::swap(v1, v2);
    }

    if (v1 == v2) {
      return false;
    }

    auto& v2_edges = edges_[v2];

    //first, process v2 -> v2
    auto v2_edge = v2_edges.begin();
    while (v2_edge != v2_edges.end()) {
      if (v2_edge->to() == v2) {
        AddEdge(v1, v1, v2_edge->label(), v2_edge->w());
        v2_edge = v2_edges.erase(v2_edge);
      } else {
        ++v2_edge;
      }
    }

    for (auto&& edge : v2_edges) {
      Vertex next = edge.to();
      assert(next != v2);
      auto to_change_lower = edges_[next].lower_bound(
        Edge(Inverse(edge.label()), v2, std::numeric_limits<Weight>::min())
      );

      auto to_change = to_change_lower;
      while (to_change != edges_[next].end() &&
          to_change->label() == Inverse(edge.label()) &&
          to_change->to() == v2
        ) {
        AddEdge(next, v1, to_change->label(), to_change->w());
        ++to_change;
      }

      edges_[next].erase(to_change_lower, to_change);
    }
    v2_edges.clear();

    return true;
  }

  void ApplyModulus(Weight new_modulus_div) {
    modulus_ = Gcd(modulus_, std::abs(new_modulus_div));
    
    if (modulus_ == 0) {
      return;
    }

    //make all weights in range [0...modulus_)

    for(auto& v_edges : edges_) {
      std::set<Edge> new_edges;
      for (auto&& e : v_edges) {
        new_edges.emplace(e.label(), e.to(), WeightMod(e.w()));
      }
      v_edges = std::move(new_edges);
    }
  }

  void FoldEdges(Vertex o, Edge e1, Edge e2) {
    static int count = 0;
    ++count;

    assert(o);

    assert(e1.label() == e2.label());
    assert(edges_.at(o).count(e1));
    assert(edges_.at(o).count(e2));

    auto weight_diff = WeightMod(e1.w() - e2.w());
    if (weight_diff) {
      if (e1.to() == e2.to()) {
        ApplyModulus(weight_diff);
      } else {
        ShiftWeight(e1.to(), weight_diff);
      }
    }

    //this->PrintAsDot(&std::ofstream("fold-" + std::to_string(count) + "-" + std::to_string(o) + "-" + std::to_string(e1.label()) + "-1.dot"));

    JoinVertices(e1.to(), e2.to());

    //this->PrintAsDot(&std::ofstream("fold-" + std::to_string(count) + "-" + std::to_string(o) + "-" + std::to_string(e1.label()) + "-2.dot"));
  }


  bool SingleFold() {
    for (auto v1 = 1u; v1 < edges_.size(); ++v1) {
      if (edges_[v1].size() < 2) {
        continue;
      }

      auto c_edge = edges_[v1].begin();
      auto n_edge = std::next(c_edge);

      while (n_edge != edges_[v1].end()) {
        if (c_edge->label() == n_edge->label()) {
          FoldEdges(v1, *c_edge, *n_edge);
          return true;
        }         
        ++c_edge, ++n_edge;
      }
    }
    return false;
  }

  void Fold() {
    while(SingleFold()) { }
  }

  Weight WeightMod(Weight w) const {
    return modulus_ == 0 ? w : (((w % modulus_) + modulus_) % modulus_);
  }

  void Harvest(size_t k, Vertex v1, Vertex v2, Weight weight, std::vector<Word>* result) const {
    std::deque<std::tuple<Vertex, Word, Weight>> current_path = {std::make_tuple(v1, Word{ }, 0)};
  
    while (!current_path.empty()) {
      Vertex v;
      Word w;
      Weight c;
      std::tie(v, w, c) = current_path.front();
      current_path.pop_front();
      if (v == v2 && (WeightMod(c - weight) == 0 || WeightMod(c + weight) == 0)) {
        if (v1 != v2 || (w.Empty() || w.GetFront() != Inverse(w.GetBack()))) {
          result->push_back(w);
        }
      }

      if (w.size() >= k) {
        continue;
      }

      for(auto&& edge : edges_[v]) {
        if (!w.Empty() && w.GetBack() == Inverse(edge.label())) {
          continue;
        }
        Word next_word = w;
        next_word.PushBack(edge.label());
        current_path.emplace_back(edge.to(), next_word, WeightMod(c + edge.w()));
      }
    }
  }

  std::vector<Word> Harvest(size_t k, Vertex v1, Vertex v2, Weight weight = 0) const {
    std::vector<Word> result;
    Harvest(k, v1, v2, weight, &result);
    std::sort(result.begin(), result.end());
    auto unique_end = std::unique(result.begin(), result.end());
    result.erase(unique_end, result.end());
    return result;
  }

  std::vector<Word> Harvest(size_t k, Weight w) const {
    std::vector<Word> result;
    for(auto v = 1u; v < edges_.size(); ++v) {
      Harvest(k, v, v, w, &result);
    }
    std::sort(result.begin(), result.end());
    auto unique_end = std::unique(result.begin(), result.end());
    result.erase(unique_end, result.end());
    return result;
  }

  void PrintAsDot(std::ostream* out) const {
  (*out) << "digraph {\n  edge [arrowtail=dot,arrowhead=vee,arrowsize=0.5];\n  1 [color=blue];\n";
  int i = -1;
  for (auto&& vertex : edges_) {
    ++i;
    if (i == 0) {
      continue;
    }
    if (vertex.empty()) {
      continue;
    }
    for (auto&& edge : vertex) {
      (*out) << i << " -> " << edge.to()
        << " [label=\"" 
        << (edge.label() % 2 ? "-" : "") << static_cast<char>('x' + static_cast<char>(edge.label() / 2)) << ";"
        << edge.w() << "\"];\n";
    }
  }

  (*out) << "}";
}

};

TEST(NaiveGraphFolding, Harvest1) {
  Graph g;
  g.PushCycle(Word("xxyx"));
  g.PushCycle(Word("xyyX"));
  g.Fold();

  auto words = g.Harvest(10, 1, 4);

  std::vector<Word> correct = {
    { 1 },
    { 0, 0, 2 },
    { 0, 2, 2, 0, 2 },
    { 0, 2, 2, 1, 1 },
    { 0, 3, 3, 0, 2 },
    { 0, 3, 3, 1, 1 },
    { 1, 3, 1, 1, 1 },
    { 0, 0, 2, 0, 0, 0, 2 },
    { 0, 2, 2, 2, 2, 0, 2 },
    { 0, 2, 2, 2, 2, 1, 1 },
    { 0, 3, 3, 3, 3, 0, 2 },
    { 0, 3, 3, 3, 3, 1, 1 },
    { 1, 3, 1, 2, 2, 0, 2 },
    { 1, 3, 1, 2, 2, 1, 1 },
    { 1, 3, 1, 3, 3, 0, 2 },
    { 1, 3, 1, 3, 3, 1, 1 },
    { 0, 0, 2, 0, 0, 2, 2, 0, 2 },
    { 0, 0, 2, 0, 0, 2, 2, 1, 1 },
    { 0, 0, 2, 0, 0, 3, 3, 0, 2 },
    { 0, 0, 2, 0, 0, 3, 3, 1, 1 },
    { 0, 2, 2, 0, 2, 0, 0, 0, 2 },
    { 0, 2, 2, 1, 1, 3, 1, 1, 1 },
    { 0, 2, 2, 2, 2, 2, 2, 0, 2 },
    { 0, 2, 2, 2, 2, 2, 2, 1, 1 },
    { 0, 3, 3, 0, 2, 0, 0, 0, 2 },
    { 0, 3, 3, 1, 1, 3, 1, 1, 1 },
    { 0, 3, 3, 3, 3, 3, 3, 0, 2 },
    { 0, 3, 3, 3, 3, 3, 3, 1, 1 },
    { 1, 3, 1, 1, 1, 3, 1, 1, 1 },
    { 1, 3, 1, 2, 2, 2, 2, 0, 2 },
    { 1, 3, 1, 2, 2, 2, 2, 1, 1 },
    { 1, 3, 1, 3, 3, 3, 3, 0, 2 },
    { 1, 3, 1, 3, 3, 3, 3, 1, 1 }
  };

  EXPECT_EQ(correct, words);
}

TEST(NaiveGraphFolding, HarvestWeight4) {
  FoldedGraph2::Word first("xxY");
  FoldedGraph2::Word second("y");

  Graph g;
  g.PushCycle(first, 1, 1);
  g.PushCycle(second, 1, 0);
  g.Fold();

  EXPECT_EQ(0, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 1), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second));

  using Word = FoldedGraph2::Word;
  EXPECT_EQ(std::vector<FoldedGraph2::Word>({
      Word("xx"),
      Word("XX"),
      Word("xxy"),
      Word("xxY"),
      Word("XXy"),
      Word("XXY"),
      Word("yxx"),
      Word("yXX"),
      Word("Yxx"),
      Word("YXX"),
      Word("xxyy"),
      Word("xxYY"),
      Word("XXyy"),
      Word("XXYY"),
      Word("yxxy"),
      Word("yXXy"),
      Word("yyxx"),
      Word("yyXX"),
      Word("YxxY"),
      Word("YXXY"),
      Word("YYxx"),
      Word("YYXX"),
    }),
    g.Harvest(4, g.root(), g.root(), 1)
  );

}

TEST(NaiveGraphFolding, StressTest) {
#ifndef NDEBUG
  static const unsigned int kRepeat = 10000;
#else
  static const unsigned int kRepeat = 100000;
#endif
  static const unsigned int kWords = 3;
  std::mt19937_64 engine(17);
  RandomWord rw(2, 8);
  std::discrete_distribution<Weight> random_weight({0.6, 0.4, 0.1});

  for (auto i = 0u; i < kRepeat; ++i) {
#ifdef DEBUG_PRINT
    std::cout << "== " << i << " ==========================" << std::endl;
#endif
    std::vector<std::pair<Word, Weight>> words;
    Graph g;
    for (auto j = 0u; j < kWords; ++j) {
      words.emplace_back(rw(engine), random_weight(engine));
#ifdef DEBUG_PRINT
      std::cout << ::testing::PrintToString(words.back()) << "\n";
#endif
      g.PushCycle(words.back().first, g.root(), words.back().second);
    }

    g.Fold();

    if (g.modulus() != 1) {
      for (auto&& w : words) {
        auto res = g.ReadWord(w.first);
        ASSERT_EQ(std::make_tuple(1, w.first.size(), g.WeightMod(w.second)), res)
          << "Pushed " << ::testing::PrintToString(words) << ", fail on " << w.first
          << ", graph mod is " << g.modulus();
      }
    }
  }
}


} //namespace naive_graph_folding

struct Cycle {
  std::pair<Word, Weight> data_;

  Cycle(const char* word, Weight weight)
    : data_(Word(word), weight)
  { }

  const Word& word() const {
    return data_.first;
  }

  Weight weight() const {
    return data_.second;
  }
};

std::ostream& operator<<(std::ostream& out, const Cycle& c) {
  return out << "{ " << c.word() << " , " << c.weight() <<" }";
}

typedef std::pair<std::vector<Cycle>, Weight> PushReadCyclesParam;
class GraphsPushReadCycles : public ::testing::TestWithParam<PushReadCyclesParam> { };

naive_graph_folding::Graph GetNaive(const std::vector<Cycle>& cycles) {
  naive_graph_folding::Graph g;
  for(auto&& cycle : cycles) {
    g.PushCycle(cycle.word(), g.root(), cycle.weight());
  }
  //g.PrintAsDot(&std::ofstream("g-" + std::to_string(0) + ".dot")); 
  g.Fold();
  return g;
}

TEST_P(GraphsPushReadCycles, NaivePushRead) {
  static int test_num = 0;
  ++test_num;
  auto g = GetNaive(GetParam().first);

  //g.PrintAsDot(&std::ofstream("g-" + std::to_string(test_num) + ".dot")); 

  EXPECT_EQ(GetParam().second, g.modulus());
  for(auto&& cycle : GetParam().first) {
    if (g.modulus() != 1) {
      EXPECT_EQ(
        std::make_tuple(g.root(), cycle.word().size(), cycle.weight()), 
        g.ReadWord(cycle.word())
      );
    } else {
      EXPECT_EQ(
        std::make_tuple(g.root(), cycle.word().size(), 0), 
        g.ReadWord(cycle.word())
      );
    }
  }
}

FoldedGraph2 GetFolded(const std::vector<Cycle>& cycles) {
  FoldedGraph2 g;
  for(auto&& cycle : cycles) {
    g.PushCycle(cycle.word(), g.root(), cycle.weight());
  }

  return g;
}

TEST_P(GraphsPushReadCycles, FoldedGraph2PushRead) {
  auto g = GetFolded(GetParam().first);

  EXPECT_EQ(GetParam().second, g.modulus());
  for(auto&& cycle : GetParam().first) {
    if (g.modulus() != 1) {
      EXPECT_EQ(
        std::make_tuple(g.root(), cycle.word().size(), cycle.weight()), 
        g.ReadWord(cycle.word())
      );
    } else {
      EXPECT_EQ(
        std::make_tuple(g.root(), cycle.word().size(), 0), 
        g.ReadWord(cycle.word())
      );
    }
  }
}

TEST_P(GraphsPushReadCycles, CompareRootHarvest) {
  auto g_folded = GetFolded(GetParam().first);
  auto g_naive = GetNaive(GetParam().first);

  auto max_length = 0u;
  for (auto& cycle : GetParam().first) {
    max_length = std::max(max_length, static_cast<decltype(max_length)>(cycle.word().size()));
  }

  auto harvest_length = max_length + 2;

  if (harvest_length > 16) {
    harvest_length = 16;
  }

  auto folded_harvest = g_folded.Harvest(max_length, 1, 1, 1);
  auto naive_harvest = g_naive.Harvest(max_length, 1, 1, 1);

  EXPECT_EQ(naive_harvest, folded_harvest);
}

PushReadCyclesParam push_read_cycles_params[] = {
    {{{"xyxYXY", 1}, {"X", 0}}, 0},
    {{{"xyxYXY", 1}, {"x", 0}}, 0},
    {{{"xyxYXY", 1}, {"y", 0}}, 0},
    {{{"xxY", 1}, {"y", 0}}, 0},
    {{{"xyxYXY", 1}, {"yy", 0}}, 0},
    {{{"x", 1}, {"X", 0}}, 1},
    {{{"Yx", 1}, {"XyxxYxY", 0}, {"yX", 0}}, 1},
    {{{"yXY", 1}, {"YYXYXyyyx", 0}, {"XXY", 0}}, 1},
    {{{"yXY", 1}, {"YYXYXyyyx", 0}, {"XXY", 0}}, 1},
    {{{"yX", 1}, {"xYxY", 1}, {"XYXyX", 1}}, 3},
    {{{"YY", 1}, {"yy", 0}, {"YXY", 1}}, 1},
    {{{"xxyXYY", 1}, {"XY", 0}, {"xyxyx", 0}}, 1},
    {{{"yXyy", 0}, {"XYY", 0}, {"Xyy", 0}}, 0},
    {{{"yyxY", 0}, {"yyXy", 0}, {"yxY", 0}}, 0},
    {{{"yXyxY", 0}, {"xxY", 0}, {"yXXX", 0}}, 0},
    {{{"xyx", 0}, {"xyyyxYXY", 0}, {"XYx", 0}}, 0},
    {{{"yxYxYxy", 0}, {"YYYYxYx", 0}, {"Xyyyy", 0}}, 0},
    {{{"YxY", 0}, {"YYxy", 1}, {"YY", 0}}, 2},
    {{{"yyy", 0}, {"YYYY", 1}, {"yyy", 1}}, 1},
    {{{"yx", 1}, {"xyx", 0}}, 0},

    //{{{"", 1}, {"", 0}}, 0},
    //{{{"", 1}, {"", 0}}, 0},
};

INSTANTIATE_TEST_CASE_P(Examples, GraphsPushReadCycles, ::testing::ValuesIn(push_read_cycles_params));



TEST(FoldedGraph2, Trivial) {
  FoldedGraph2 g;
  EXPECT_EQ(0, g.vertex(g.root()).endpoint(0));
  EXPECT_EQ(0, g.vertex(g.root()).endpoint(1));
  EXPECT_EQ(0, g.vertex(g.root()).endpoint(2));
  EXPECT_EQ(0, g.vertex(g.root()).endpoint(3));
}

FoldedGraph2::Word AlphasToLabels(std::vector<int> alphas) {
  FoldedGraph2::Word result;
  for (auto alpha : alphas) {
    switch(alpha) {
      case 1:
        result.PushBack(0);
        break;
      case -1:
        result.PushBack(1);
        break;
      case 2:
        result.PushBack(2);
        break;
      case -2:
        result.PushBack(3);
        break;
      default:
        assert(false);
    }
  }
  return result;
}

TEST(FoldedGraph2, ReadWord) {
  FoldedGraph2 g;
  EXPECT_EQ(std::make_tuple(1, 0, 0), g.ReadWord(AlphasToLabels({1, 2})));

  EXPECT_EQ(3, g.PushWord(AlphasToLabels({1, 2})));
  EXPECT_EQ(std::make_tuple(3, 2, 0), g.ReadWord(AlphasToLabels({1, 2})));
  EXPECT_EQ(std::make_tuple(2, 1, 0), g.ReadWord(AlphasToLabels({1, 1})));

  EXPECT_EQ(std::make_tuple(1, 2, 0), g.ReadInverse(AlphasToLabels({1, 2}), 3));
  EXPECT_EQ(std::make_tuple(2, 1, 0), g.ReadInverse(AlphasToLabels({2, 2}), 3));
}

TEST(FoldedGraph2, PushWord) {
  FoldedGraph2 g;
  EXPECT_EQ(3, g.PushWord(AlphasToLabels({1, 2})));

  EXPECT_EQ(2, g.vertex(1).endpoint(0));
  EXPECT_EQ(0, g.vertex(1).endpoint(2));
  EXPECT_EQ(0, g.vertex(1).endpoint(1));
  EXPECT_EQ(0, g.vertex(1).endpoint(3));

  EXPECT_EQ(0, g.vertex(2).endpoint(0));
  EXPECT_EQ(3, g.vertex(2).endpoint(2));
  EXPECT_EQ(1, g.vertex(2).endpoint(1));
  EXPECT_EQ(0, g.vertex(2).endpoint(3));

  EXPECT_EQ(0, g.vertex(3).endpoint(0));
  EXPECT_EQ(0, g.vertex(3).endpoint(2));
  EXPECT_EQ(0, g.vertex(3).endpoint(1));
  EXPECT_EQ(2, g.vertex(3).endpoint(3));

  EXPECT_EQ(3, g.PushWord(AlphasToLabels({1, 2})));
  EXPECT_EQ(4, g.PushWord(AlphasToLabels({1, 1})));

  EXPECT_EQ(2, g.vertex(1).endpoint(0));
  EXPECT_EQ(0, g.vertex(1).endpoint(2));
  EXPECT_EQ(0, g.vertex(1).endpoint(1));
  EXPECT_EQ(0, g.vertex(1).endpoint(3));

  EXPECT_EQ(4, g.vertex(2).endpoint(0));
  EXPECT_EQ(3, g.vertex(2).endpoint(2));
  EXPECT_EQ(1, g.vertex(2).endpoint(1));
  EXPECT_EQ(0, g.vertex(2).endpoint(3));

  EXPECT_EQ(0, g.vertex(3).endpoint(0));
  EXPECT_EQ(0, g.vertex(3).endpoint(2));
  EXPECT_EQ(0, g.vertex(3).endpoint(1));
  EXPECT_EQ(2, g.vertex(3).endpoint(3));

  EXPECT_EQ(0, g.vertex(4).endpoint(0));
  EXPECT_EQ(0, g.vertex(4).endpoint(2));
  EXPECT_EQ(2, g.vertex(4).endpoint(1));
  EXPECT_EQ(0, g.vertex(4).endpoint(3));

  EXPECT_EQ(5, g.PushWord(AlphasToLabels({1, 2}), 2));

  EXPECT_EQ(2, g.vertex(1).endpoint(0));
  EXPECT_EQ(0, g.vertex(1).endpoint(2));
  EXPECT_EQ(0, g.vertex(1).endpoint(1));
  EXPECT_EQ(0, g.vertex(1).endpoint(3));

  EXPECT_EQ(4, g.vertex(2).endpoint(0));
  EXPECT_EQ(3, g.vertex(2).endpoint(2));
  EXPECT_EQ(1, g.vertex(2).endpoint(1));
  EXPECT_EQ(0, g.vertex(2).endpoint(3));

  EXPECT_EQ(0, g.vertex(3).endpoint(0));
  EXPECT_EQ(0, g.vertex(3).endpoint(2));
  EXPECT_EQ(0, g.vertex(3).endpoint(1));
  EXPECT_EQ(2, g.vertex(3).endpoint(3));

  EXPECT_EQ(0, g.vertex(4).endpoint(0));
  EXPECT_EQ(5, g.vertex(4).endpoint(2));
  EXPECT_EQ(2, g.vertex(4).endpoint(1));
  EXPECT_EQ(0, g.vertex(4).endpoint(3));

  EXPECT_EQ(0, g.vertex(5).endpoint(0));
  EXPECT_EQ(0, g.vertex(5).endpoint(2));
  EXPECT_EQ(0, g.vertex(5).endpoint(1));
  EXPECT_EQ(4, g.vertex(5).endpoint(3));

}

TEST(FoldedGraph2, JoinVertices1) {
  FoldedGraph2 g;
  EXPECT_EQ(2, g.PushWord({0}));
  EXPECT_EQ(3, g.PushWord({2}));

  g.JoinVertices(2, 3);
  EXPECT_TRUE(g.Equal(2, 3));

  EXPECT_EQ(0, g.vertex(2).endpoint(0));
  EXPECT_EQ(0, g.vertex(2).endpoint(2));
  EXPECT_EQ(1, g.vertex(2).endpoint(1));
  EXPECT_EQ(1, g.vertex(2).endpoint(3));

  EXPECT_EQ(std::make_tuple(2, 1, 0), g.ReadWord({0}));
  EXPECT_EQ(std::make_tuple(2, 1, 0), g.ReadWord({2}));
}

TEST(FoldedGraph2, JoinVertices2) {
  FoldedGraph2 g;
  EXPECT_EQ(2, g.PushWord({0}));

  g.JoinVertices(1, 2);
  EXPECT_TRUE(g.Equal(2, 1));

  EXPECT_EQ(std::make_tuple(1, 1, 0), g.ReadWord({0}));
  EXPECT_EQ(std::make_tuple(1, 0, 0), g.ReadWord({2}));
  EXPECT_EQ(std::make_tuple(1, 1, 0), g.ReadWord({0, 2}));
  EXPECT_EQ(std::make_tuple(1, 2, 0), g.ReadWord({0, 0, 2, 0}));
}

TEST(FoldedGraph2, PushCycle1) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({0, 2}));

  EXPECT_EQ(std::make_tuple(1, 2, 0), g.ReadWord({0, 2}));
  EXPECT_EQ(std::make_tuple(1, 4, 0), g.ReadWord({0, 2, 0, 2}));
}

TEST(FoldedGraph2, PushCycle2) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({0, 2, 0}));
  EXPECT_EQ(true, g.PushCycle({0, 0}));

  EXPECT_EQ(std::make_tuple(1, 2, 0), g.ReadWord({0, 0}));
  EXPECT_EQ(std::make_tuple(1, 3, 0), g.ReadWord({0, 2, 0}));
  EXPECT_EQ(std::make_tuple(1, 4, 0), g.ReadWord({0, 2, 2, 0}));
  EXPECT_EQ(std::make_tuple(1, 5, 0), g.ReadWord({0, 2, 2, 2, 0}));
}

TEST(FoldedGraph2, PushCycle3) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({0, 2, 0}));
  EXPECT_EQ(std::make_tuple(1, 3, 0), g.ReadWord({0, 2, 0}));

  EXPECT_EQ(true, g.PushCycle({0, 2, 0, 2, 0}));
  EXPECT_EQ(std::make_tuple(1, 2, 0), g.ReadWord({2, 0}));
}

TEST(FoldedGraph2, PushCycle4) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({0, 2, 0}));
  EXPECT_EQ(false, g.PushCycle({0, 2, 0}));

  EXPECT_EQ(std::make_tuple(1, 3, 0), g.ReadWord({0, 2, 0}));

  EXPECT_EQ(true, g.PushCycle({0, 2}));
  EXPECT_TRUE(g.Equal(1, 2));
  EXPECT_TRUE(g.Equal(2, 3));
}

TEST(FoldedGraph2, PushCycle5) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({1, 0, 0}));
}

TEST(FoldedGraph2, PushCycle6) {
  FoldedGraph2 g;
  EXPECT_EQ(true, g.PushCycle({0, 0, 0}));
  EXPECT_EQ(true, g.PushCycle({1, 1}));
}


TEST(FoldedGraph2, Harvest1) {
  FoldedGraph2 g;
  g.PushCycle(AlphasToLabels({1, 1, 2, 1}));
  g.PushCycle(AlphasToLabels({1, 2, 2, -1}));

  auto words = g.Harvest(10, 1, 4, 0);

  std::vector<FoldedGraph2::Word> correct = {
    { 1 },
    { 0, 0, 2 },
    { 0, 2, 2, 0, 2 },
    { 0, 3, 3, 0, 2 },
    { 1, 3, 1, 1, 1 },
    { 0, 0, 2, 0, 0, 0, 2 },
    { 0, 2, 2, 2, 2, 0, 2 },
    { 0, 3, 3, 3, 3, 0, 2 },
    { 1, 3, 1, 2, 2, 0, 2 },
    { 1, 3, 1, 2, 2, 1, 1 },
    { 1, 3, 1, 3, 3, 0, 2 },
    { 1, 3, 1, 3, 3, 1, 1 },
    { 0, 0, 2, 0, 0, 2, 2, 0, 2 },
    { 0, 0, 2, 0, 0, 3, 3, 0, 2 },
    { 0, 2, 2, 0, 2, 0, 0, 0, 2 },
    { 0, 2, 2, 2, 2, 2, 2, 0, 2 },
    { 0, 3, 3, 0, 2, 0, 0, 0, 2 },
    { 0, 3, 3, 3, 3, 3, 3, 0, 2 },
    { 1, 3, 1, 1, 1, 3, 1, 1, 1 },
    { 1, 3, 1, 2, 2, 2, 2, 0, 2 },
    { 1, 3, 1, 2, 2, 2, 2, 1, 1 },
    { 1, 3, 1, 3, 3, 3, 3, 0, 2 },
    { 1, 3, 1, 3, 3, 3, 3, 1, 1 }
  };

  EXPECT_EQ(correct, words);
}

TEST(FoldedGraph2, PushWithWeight1) {
  FoldedGraph2::Word first("xyxYXY");
  FoldedGraph2::Word second("X");

  FoldedGraph2 g;
  g.PushCycle(first, g.root(), 1);
  g.PushCycle(second, g.root(), 0);

  EXPECT_EQ(0, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 1), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second));
}

TEST(FoldedGraph2, PushWithWeight2) {
  FoldedGraph2::Word first("xyxYXY");
  FoldedGraph2::Word second("x");

  FoldedGraph2 g;
  g.PushCycle(first, g.root(), 1);
  g.PushCycle(second, g.root(), 0);

  EXPECT_EQ(0, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 1), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second));
}

TEST(FoldedGraph2, PushWithWeight3) {
  FoldedGraph2::Word first("xyxYXY");
  FoldedGraph2::Word second("y");

  FoldedGraph2 g;
  g.PushCycle(first, g.root(), 1);
  g.PushCycle(second, g.root(), 0);

  EXPECT_EQ(0, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 1), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second));
}

TEST(FoldedGraph2, HarvestWeight4) {
  FoldedGraph2::Word first("xxY");
  FoldedGraph2::Word second("y");

  FoldedGraph2 g;
  g.PushCycle(first, g.root(), 1);
  g.PushCycle(second, g.root(), 0);

  std::ofstream out("graph.dot");
  g.PrintAsDot(&out);

  EXPECT_EQ(0, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 1), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second));

  using Word = FoldedGraph2::Word;
  EXPECT_EQ(std::vector<FoldedGraph2::Word>({
      Word("xx"),
      Word("XX"),
      Word("xxy"),
      Word("xxY"),
      Word("XXy"),
      Word("XXY"),
      Word("yxx"),
      Word("yXX"),
      Word("Yxx"),
      Word("YXX"),
      Word("xxyy"),
      Word("xxYY"),
      Word("XXyy"),
      Word("XXYY"),
      Word("yxxy"),
      Word("yXXy"),
      Word("yyxx"),
      Word("yyXX"),
      Word("YxxY"),
      Word("YXXY"),
      Word("YYxx"),
      Word("YYXX"),
    }),
    g.Harvest(4, g.root(), g.root(), 1)
  );

}

TEST(FoldedGraph2, HarvestWeight5) {
  FoldedGraph2::Word first("xxY");
  FoldedGraph2::Word second("y");

  FoldedGraph2 g;
  g.PushCycle(first, g.root(), 1);
  g.PushCycle(second, 3, 0);

  EXPECT_EQ(0, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 1), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second, 3));

  using Word = FoldedGraph2::Word;

  EXPECT_EQ(0, g.modulus());
  EXPECT_EQ(std::vector<Word>({
      Word("xx"),
      Word("XX"),
      Word("xxy"),
      Word("xxY"),
      Word("XXy"),
      Word("XXY"),
      Word("yxx"),
      Word("yXX"),
      Word("Yxx"),
      Word("YXX"),
      Word("xxyy"),
      Word("xxYY"),
      Word("XXyy"),
      Word("XXYY"),
      Word("yxxy"),
      Word("yXXy"),
      Word("yyxx"),
      Word("yyXX"),
      Word("YxxY"),
      Word("YXXY"),
      Word("YYxx"),
      Word("YYXX"),
    }),
    g.Harvest(4, g.root(), g.root(), 1)
  );
}

TEST(FoldedGraph2, HarvestWeight6) {
  FoldedGraph2::Word first("xyxYXY");
  FoldedGraph2::Word second("yy");

  FoldedGraph2 g;
  g.PushCycle(first, g.root(), 1);
  g.PushCycle(second, g.root(), 0);

  EXPECT_EQ(0, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 1), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second));

  using Word = FoldedGraph2::Word;
  EXPECT_EQ(std::vector<Word>({
      Word("xyxYXy"),
      Word("xyxYXY"),
      Word("yxyXYX"),
      Word("YxyXYX"),
    }),
    g.Harvest(6, g.root(), g.root(), 1)
  );
}

TEST(FoldedGraph2, PushWithWeight7) {
  FoldedGraph2::Word first("x");
  FoldedGraph2::Word second("X");

  FoldedGraph2 g;
  g.PushCycle(first, g.root(), 1);
  g.PushCycle(second, g.root(), 0);

  EXPECT_EQ(1, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 0), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second));
}


TEST(FoldedGraph2, PushWithWeight8) {
  FoldedGraph2::Word first("Yx");
  FoldedGraph2::Word second("XyxxYxY");
  FoldedGraph2::Word third("yX");

  FoldedGraph2 g;
  g.PushCycle(first, g.root(), 1);
  g.PushCycle(second, g.root(), 0);
  g.PushCycle(third, g.root(), 0);

  EXPECT_EQ(1, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 0), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second));
  EXPECT_EQ(std::make_tuple(g.root(), third.size(), 0), g.ReadWord(third));
}

TEST(FoldedGraph2, PushWithWeight9) {
  FoldedGraph2::Word first("yXY");
  FoldedGraph2::Word second("YYXYXyyyx");
  FoldedGraph2::Word third("XXY");

  FoldedGraph2 g;
  g.PushCycle(first, g.root(), 1);
  g.PushCycle(second, g.root(), 0);
  g.PushCycle(third, g.root(), 0);

  EXPECT_EQ(1, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 0), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second));
  EXPECT_EQ(std::make_tuple(g.root(), third.size(), 0), g.ReadWord(third));
}

TEST(FoldedGraph2, PushWithWeight10) {
  FoldedGraph2::Word first("xxyXYY");
  FoldedGraph2::Word second("XY");
  FoldedGraph2::Word third("xyxyx");

  FoldedGraph2 g;
  g.PushCycle(first, g.root(), 1);
  g.PushCycle(second, g.root(), 0);
  g.PushCycle(third, g.root(), 0);

  EXPECT_EQ(1, g.modulus());
  EXPECT_EQ(std::make_tuple(g.root(), first.size(), 0), g.ReadWord(first));
  EXPECT_EQ(std::make_tuple(g.root(), second.size(), 0), g.ReadWord(second));
  EXPECT_EQ(std::make_tuple(g.root(), third.size(), 0), g.ReadWord(third));
}



TEST(FoldedGraph2, PushCycleStressRandom) {
#ifndef NDEBUG
  static const unsigned int kRepeat = 1000;
#else
  static const unsigned int kRepeat = 1000000;
#endif
  static const unsigned int kWords = 3;
  std::mt19937_64 engine;
  std::uniform_int_distribution<> random_letter(0, 3);
  std::uniform_int_distribution<size_t> random_length(2, 10);

  for (auto i = 0u; i < kRepeat; ++i) {
#ifdef DEBUG_PRINT
    std::cout << "== " << i << " ==========================" << std::endl;
#endif
    FoldedGraph2 g;
    for (auto j = 0u; j < kWords; ++j) {
      FoldedGraph2::Word w;
      size_t length = random_length(engine);
      while(w.size() < length) {
        w.PushBack(random_letter(engine));
      }

#ifdef DEBUG_PRINT
      std::cout << w << std::endl;
#endif

      g.PushCycle(w);


      ASSERT_EQ(std::make_tuple(1, length, 0), g.ReadWord(w));
    }
  }
}

TEST(FoldedGraph2, PushCycleWithWeightStressRandom) {
  static const auto kDuration = std::chrono::seconds(10);
  static const unsigned int kWords = 4;
  std::mt19937_64 engine(17);
  RandomWord rw(2, 10);
  std::discrete_distribution<Weight> random_weight({0.6, 0.4, 0.1});

  auto begin = std::chrono::steady_clock::now();

  auto repeat = 0ull;
  while (std::chrono::steady_clock::now() - begin < kDuration) {
    ++repeat;
#ifdef DEBUG_PRINT
    std::cout << "== " << i << " ==========================" << std::endl;
#endif
    std::vector<std::pair<Word, Weight>> words;
    FoldedGraph2 g;
    for (auto j = 0u; j < kWords; ++j) {
      words.emplace_back(rw(engine), random_weight(engine));
#ifdef DEBUG_PRINT
      std::cout << ::testing::PrintToString(words.back()) << "\n";
#endif
      g.PushCycle(words.back().first, g.root(), words.back().second);
    }

    if (g.modulus() != 1) {
      for (auto&& w : words) {
        auto res = g.ReadWord(w.first);
        ASSERT_EQ(std::make_tuple(1, w.first.size(), g.WeightMod(w.second)), res)
          << "Pushed " << ::testing::PrintToString(words) << ", fail on " << w.first
          << ", graph mod is " << g.modulus();
      }
    }
  }
  std::cout << std::string(13, ' ') << repeat << " repeats" << std::endl;
  ASSERT_GT(repeat, 10000);
}

TEST(FoldedGraph2, StressRootHarvestCompareWithNaive) {
  static const auto kDuration = std::chrono::seconds(10);
  static const unsigned int kWords = 2;
  std::mt19937_64 engine(17);
  RandomWord rw(2, 4);
  std::discrete_distribution<Weight> random_weight({0.6, 0.4, 0.1});

  auto begin = std::chrono::steady_clock::now();

  std::chrono::high_resolution_clock::duration harvest_folded_duration{}, harvest_naive_duration{};
  std::chrono::high_resolution_clock::time_point proc_begin;
  auto repeat = 0ull;
  while (std::chrono::steady_clock::now() - begin < kDuration) {
    ++repeat;
#ifdef DEBUG_PRINT
    std::cout << "== " << i << " ==========================" << std::endl;
#endif
    std::vector<std::pair<Word, Weight>> words;
    FoldedGraph2 g;
    naive_graph_folding::Graph g_naive;
    auto max_length = 0u;
    for (auto j = 0u; j < kWords; ++j) {
      words.emplace_back(rw(engine), random_weight(engine));
#ifdef DEBUG_PRINT
      std::cout << ::testing::PrintToString(words.back()) << "\n";
#endif
      proc_begin = std::chrono::high_resolution_clock::now();
      g.PushCycle(words.back().first, g.root(), words.back().second);
      harvest_folded_duration += (std::chrono::high_resolution_clock::now() - proc_begin);
  
      proc_begin = std::chrono::high_resolution_clock::now();
      g_naive.PushCycle(words.back().first, g.root(), words.back().second);
      harvest_naive_duration += (std::chrono::high_resolution_clock::now() - proc_begin);
      if (words.back().first.size() > max_length) {
        max_length = words.back().first.size();
      }
    }

    proc_begin = std::chrono::high_resolution_clock::now();
    g_naive.Fold();
    harvest_naive_duration += (std::chrono::high_resolution_clock::now() - proc_begin);

    proc_begin = std::chrono::high_resolution_clock::now();
    auto harvest_folded = g.Harvest(max_length + 2, g.root(), g.root(), 1);
    harvest_folded_duration += (std::chrono::high_resolution_clock::now() - proc_begin);

    proc_begin = std::chrono::high_resolution_clock::now();
    auto harvest_naive = g_naive.Harvest(max_length + 2, g.root(), g.root(), 1);
    harvest_naive_duration += (std::chrono::high_resolution_clock::now() - proc_begin);

    ASSERT_EQ(harvest_naive, harvest_folded)
      << "Pushed " << ::testing::PrintToString(words);

  }
  std::cout << std::string(13, ' ') << repeat << " repeats" << std::endl;
  std::cout << std::string(13, ' ') 
    << std::chrono::duration_cast<std::chrono::milliseconds>(harvest_folded_duration).count()
    << " vs "
    << std::chrono::duration_cast<std::chrono::milliseconds>(harvest_naive_duration).count()
    << std::endl;
  ASSERT_GT(repeat, 10000);

}

TEST(FoldedGraph2, StressFullHarvestCompareWithNaive) {
  static const auto kDuration = std::chrono::seconds(10);
  static const unsigned int kWords = 2;
  std::mt19937_64 engine(17);
  RandomWord rw(2, 4);
  std::discrete_distribution<Weight> random_weight({0.6, 0.4, 0.1});

  auto begin = std::chrono::steady_clock::now();

  std::chrono::high_resolution_clock::duration harvest_folded_duration{}, harvest_naive_duration{};
  std::chrono::high_resolution_clock::time_point proc_begin;
  auto repeat = 0ull;
  auto non_trivial = 0ull;
  while (std::chrono::steady_clock::now() - begin < kDuration) {
    ++repeat;
#ifdef DEBUG_PRINT
    std::cout << "== " << i << " ==========================" << std::endl;
#endif
    std::vector<std::pair<Word, Weight>> words;
    FoldedGraph2 g;
    naive_graph_folding::Graph g_naive;
    auto max_length = 0u;
    for (auto j = 0u; j < kWords; ++j) {
      words.emplace_back(rw(engine), random_weight(engine));
#ifdef DEBUG_PRINT
      std::cout << ::testing::PrintToString(words.back()) << "\n";
#endif
      proc_begin = std::chrono::high_resolution_clock::now();
      g.PushCycle(words.back().first, g.root(), words.back().second);
      harvest_folded_duration += (std::chrono::high_resolution_clock::now() - proc_begin);

      proc_begin = std::chrono::high_resolution_clock::now();
      g_naive.PushCycle(words.back().first, g.root(), words.back().second);
      harvest_naive_duration += (std::chrono::high_resolution_clock::now() - proc_begin);
      if (words.back().first.size() > max_length) {
        max_length = words.back().first.size();
      }
    }

    if (g.modulus() != 1) {
      ++non_trivial;
    }

    proc_begin = std::chrono::high_resolution_clock::now();
    g_naive.Fold();
    harvest_naive_duration += (std::chrono::high_resolution_clock::now() - proc_begin);

    proc_begin = std::chrono::high_resolution_clock::now();
    auto harvest_folded = g.Harvest(max_length + 2, 1);
    harvest_folded_duration += (std::chrono::high_resolution_clock::now() - proc_begin);

    proc_begin = std::chrono::high_resolution_clock::now();
    auto harvest_naive = g_naive.Harvest(max_length + 2, 1);
    harvest_naive_duration += (std::chrono::high_resolution_clock::now() - proc_begin);

    ASSERT_EQ(ReduceAndMinCycle(harvest_naive), ReduceAndMinCycle(harvest_folded))
      << "Pushed " << ::testing::PrintToString(words);

  }
  std::cout << std::string(13, ' ') << repeat << " repeats" << std::endl;
  std::cout << std::string(13, ' ') << non_trivial << " times modulus != 1" << std::endl;
  std::cout << std::string(13, ' ') 
    << std::chrono::duration_cast<std::chrono::milliseconds>(harvest_folded_duration).count()
    << " vs "
    << std::chrono::duration_cast<std::chrono::milliseconds>(harvest_naive_duration).count()
    << std::endl;
  ASSERT_GT(repeat, 5000);

}


} //namespace

} //namespace crag
