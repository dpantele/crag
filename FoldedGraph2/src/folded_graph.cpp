#include <algorithm>
#include <deque>
#include <fstream>
#include <limits>
#include <deque>
#include <string>
#include <type_traits>

#include "folded_graph2.h"

namespace crag {


using Vertex = FoldedGraph2::Vertex;
using VertexEdges = FoldedGraph2::VertexEdges;
using Label = FoldedGraph2::Label;
using Word = FoldedGraph2::Word;
using Weight = FoldedGraph2::Weight;

Weight Gcd(Weight a, Weight b) {
  while (b != 0) {
    auto c = a % b;
    a = b;
    b = c;
  }

  return a;
}

Vertex VertexEdges::endpoint(Label l) const {
  assert(l < 2 * kAlphabetSize);
  return edges_[l];
}

Vertex& VertexEdges::endpoint(Label l) {
  assert(l < 2 * kAlphabetSize);
  return edges_[l];
}

Weight VertexEdges::weight(Label l) const  {
  assert(l < 2 * kAlphabetSize);
  return weights_[l];
}

std::tuple<Vertex, Word::size_type, Weight> FoldedGraph2::ReadWord(Word w, Word::size_type length_limit, Vertex s) const {
  Weight cur_weight = 0;
  if (w.Empty()) {
    return std::make_tuple(s, 0u, cur_weight);
  }

  assert(s);

  s = GetLastCombinedWith(s);

  if (length_limit > w.size()) {
    length_limit = w.size();
  }

  for (auto i = 0u; i < length_limit; ++i) {
    auto next = GetLastCombinedWith(vertex(s).endpoint(w.GetFront()));
    cur_weight += vertex(s).weight(w.GetFront());
    if (next) {
      assert(WeightMod(vertex(s).weight(w.GetFront()) + vertex(next).weight(Inverse(w.GetFront()))) == 0);
      s = next;
    } else {
      return std::make_tuple(s, i, WeightMod(cur_weight));
    }
    w.PopFront();
  }

  return std::make_tuple(s, length_limit, WeightMod(cur_weight));
}

std::tuple<Vertex, Word::size_type, Weight> FoldedGraph2::ReadWord(Word w, Vertex s) const {
  return ReadWord(w, w.size(), s);
}

std::tuple<Vertex, Word::size_type, Weight> FoldedGraph2::ReadWord(Word w) const {
  return ReadWord(w, w.size(), root());
}

std::tuple<Vertex, Word::size_type, Weight> FoldedGraph2::ReadInverse(Word w, Word::size_type length_limit, Vertex s) const {
  if (w.Empty()) {
    return std::make_tuple(s, 0u, 0);
  }

  Weight cur_weight = 0;

  assert(s);

  s = GetLastCombinedWith(s);

  if (length_limit > w.size()) {
    length_limit = w.size();
  }

  for (auto i = 0u; i < length_limit; ++i) {
    auto next = GetLastCombinedWith(vertex(s).endpoint(Inverse(w.GetBack())));
    cur_weight += -vertex(s).weight(Inverse(w.GetBack()));
    if (next) {
      assert(vertex(s).weight(Inverse(w.GetBack())) == -vertex(next).weight(w.GetBack()));
      s = next;
    } else {
      return std::make_tuple(s, i, WeightMod(cur_weight));
    }
    w.PopBack();
  }

  return std::make_tuple(s, length_limit, WeightMod(cur_weight));
}

std::tuple<Vertex, Word::size_type, Weight> FoldedGraph2::ReadInverse(Word w, Vertex s) const {
  return ReadInverse(w, w.size(), s);
}

std::tuple<Vertex, Word::size_type, Weight> FoldedGraph2::ReadInverse(Word w) const {
  return ReadInverse(w, w.size(), root());
}


/**
* Follows the references of combined references and update them to 1-step
*/
Vertex FoldedGraph2::GetLastCombinedWith(Vertex v) const {
  if (v == kNullVertex) {
    return kNullVertex;
  }

  Vertex& combined_with = edges_[v].combined_with_;

  if (combined_with == kNullVertex) {
    return v;
  } else {
    return combined_with = GetLastCombinedWith(combined_with);
  }
}

const VertexEdges& FoldedGraph2::vertex(Vertex v) const {
  assert(v);
  assert(v < edges_.size());
  return edges_[GetLastCombinedWith(v)];
}

//VertexEdges& FoldedGraph2::vertex_mod(Vertex v) {
//  assert(v);
//  assert(v < edges_.size());
//  return edges_[GetLastCombinedWith(v)];
//}

Vertex FoldedGraph2::AddEdge(Label l, Vertex from, Vertex to) {
  if (to == kNullVertex) {
    to = vertex(from).endpoint(l);
  }

  if (to == kNullVertex) {
    assert(edges_.size() <= std::numeric_limits<int>::max());
    to = static_cast<int>(edges_.size());
    edges_.emplace_back();
  }

  assert(
      vertex(from).endpoint(l) == kNullVertex ||
      Equal(vertex(from).endpoint(l), GetLastCombinedWith(to))
  );

  edges_[GetLastCombinedWith(from)].endpoint(l) = GetLastCombinedWith(to);

  assert(
      vertex(to).endpoint(Inverse(l)) == kNullVertex ||
      Equal(vertex(to).endpoint(Inverse(l)), GetLastCombinedWith(from))
  );

  edges_[GetLastCombinedWith(to)].endpoint(Inverse(l)) = GetLastCombinedWith(from);

  return to;
}

Vertex FoldedGraph2::PushWord(Word w, Vertex s, Weight weight) {
  s = GetLastCombinedWith(s);
  if (w.Empty()) {
    return s;
  }

  Word::size_type existing_length{};
  Weight current_weight;

  std::tie(s, existing_length, current_weight) = ReadWord(w, s);
  w.PopFront(existing_length);
  while (!w.Empty()) {
    auto next = AddEdge(w.GetFront(), s);
    if (!WeightMod(current_weight - weight)) {
      assert(vertex(s).weights_[w.GetFront()] == 0);
      edges_[s].weights_[w.GetFront()] = weight - current_weight;
      edges_[next].weights_[w.GetFront()] = -(weight - current_weight);
      current_weight = weight;
    }

    s = next;
    w.PopFront();
  }

  modulus_ = Gcd(std::abs(current_weight - weight), modulus_);

  assert(s > 1);
  return s;
}

Vertex FoldedGraph2::FindInconsistentWeights() const {
  for(auto&& vertex : edges_) {
    if (vertex.combined_with_) {
      continue;
    }

    for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      auto& neigbour = edges_[GetLastCombinedWith(vertex.endpoint(label))];
      if (WeightMod(vertex.weight(label) + neigbour.weight(Inverse(label)))) {
        return GetLastCombinedWith(neigbour.endpoint(Inverse(label)));
      }
    }
  }
  return kNullVertex;
}

bool FoldedGraph2::Equal(Vertex v1, Vertex v2) {
  if (v1 == v2) {
    return true;
  }

  if (!v1 || !v2) {
    return false;
  }

  return GetLastCombinedWith(v1) == GetLastCombinedWith(v2) && GetLastCombinedWith(v1) != kNullVertex;
}

void FoldedGraph2::JoinVertices(Vertex v1, Vertex v2) {
  if (Equal(v1, v2)) {
    return;
  }

  std::deque<std::tuple<Vertex, Vertex, Vertex, Label>> vertex_to_join = {std::make_tuple(v1, v2, kNullVertex, 0)};
  auto step = 0u;
  while (!vertex_to_join.empty()) {
    ++step;
    std::ofstream out("join-" + std::to_string(step) + ".dot");
    this->PrintAsDot(&out);
    out.close();

    Vertex first, second, prev;
    Label through;

    std::tie(first, second, prev, through) = vertex_to_join.front();
    vertex_to_join.pop_front();

    first = GetLastCombinedWith(first);
    second = GetLastCombinedWith(second);
    prev = GetLastCombinedWith(prev);

    if (first == second) {
      for(auto label = 0u; label < kAlphabetSize * 2; ++label) {
        auto n = vertex(first).endpoint(label);
        if (n) {
          assert(vertex(n).endpoint(Inverse(label)));
          modulus_ = Gcd(std::abs(vertex(first).weight(label) + vertex(n).weight(Inverse(label))), modulus_);
        }
      }
      continue;
    }
    
    if (second < first) {
      std::swap(first, second);
    }

    auto& first_edges = edges_[first].edges_;
    auto& second_edges = edges_[second].edges_;
    auto& first_weights = edges_[first].weights_;
    const auto& second_weights = edges_[second].weights_;

    for (auto label = 0u; label < kAlphabetSize * 2; ++label) {
      if (first_edges[label] && second_edges[label]) {
        vertex_to_join.emplace_back(first_edges[label], second_edges[label], first, label);

        auto weight_diff = WeightMod(first_weights[label] - second_weights[label]); //eps
        //if (label & 1) {
          // x --(-l, w)--> a
          // y --(-l, v)--> b
          // want to make w - eps = w - w + v = v
          // hence weight(-l) should be decreased by eps
          // hence weight(l) should be increased by eps
          this->ShiftWeight(first_edges[label], weight_diff);
        //} else {
          // x --( l, w)--> a
          // y --( l, v)--> b
          // want to make w - eps = w - w + v = v
          // hence weight(l) should be decreased by eps
          // hence weight(-l) should be increased by eps
        //  this->ShiftWeight(first_edges[label], weight_diff);
        //} 
      }
    }

    assert(vertex(second).combined_with_ == kNullVertex);
    vertex(second).combined_with_ = first;

    for (auto label = 0u; label < kAlphabetSize * 2; ++label) {
      if (!first_edges[label]) {
        first_edges[label] = second_edges[label] ? GetLastCombinedWith(second_edges[label]) : second_edges[label];
        assert(!first_weights[label]);
        first_weights[label] = second_weights[label];
      } else {
        first_edges[label] = GetLastCombinedWith(first_edges[label]);
      }

      if (second_edges[label]) {
        assert(first_edges[label]);
        modulus_ = Gcd(std::abs(first_weights[label] - second_weights[label]), modulus_);
      }
    }
  }

  ++step;
  std::ofstream out("join-" + std::to_string(step) + ".dot");
  this->PrintAsDot(&out);
  out.close();

  assert(FindInconsistentWeights() == kNullVertex);
}

bool FoldedGraph2::PushCycle(Word w, Vertex s, Weight weight) {
  if (w.Empty()) {
    return false;
  }

  s = GetLastCombinedWith(s);

  Word::size_type existing_prefix_length{}, existing_suffix_length{};
  Vertex prefix_end{}, suffix_begin{};
  Weight current_weight = 0;
  Weight suffix_weight = 0;

  std::tie(suffix_begin, existing_suffix_length, suffix_weight) = ReadInverse(w, s);
  auto without_suffix = w;
  without_suffix.PopBack(existing_suffix_length);

  //now we need to write a word without_suffix from s to suffix_begin
  std::tie(prefix_end, existing_prefix_length, current_weight) = ReadWord(without_suffix, s);
  current_weight += suffix_weight;
  auto new_word = without_suffix;
  new_word.PopFront(existing_prefix_length);

  if (!new_word.Empty()) { //if we can't read the whole word
    //we will add new vertices from prefix_end to suffix_begin
    //first add a 'stem' if they coinside

    if (prefix_end == suffix_begin) {
      while(new_word.GetFront() == Inverse(new_word.GetBack())) {
        prefix_end = AddEdge(new_word.GetFront(), prefix_end);
        suffix_begin = prefix_end;
        new_word.PopFront();
        new_word.PopBack();
      }
    }
    assert(!new_word.Empty()); //it can't be empty, because the word itself should be reduced

    while (new_word.size() > 1) {
      prefix_end = AddEdge(new_word.GetFront(), prefix_end);
      new_word.PopFront();
    }

    assert(vertex(prefix_end).endpoint(new_word.GetFront()) == kNullVertex);
    AddEdge(new_word.GetFront(), prefix_end, suffix_begin);
    assert(edges_[prefix_end].weights_[new_word.GetFront()] == 0);
    edges_[prefix_end].weights_[new_word.GetFront()] = weight - current_weight;
    edges_[suffix_begin].weights_[Inverse(new_word.GetFront())] = current_weight - weight;
    return true;
  }

  if (prefix_end == suffix_begin) {
    modulus_ = Gcd(std::abs(current_weight - weight), modulus_);
    return false;
  }

  std::ofstream out("compl.dot");
  this->PrintAsDot(&out);
  out.close();

  auto weight_diff = WeightMod(weight - current_weight); // = eps
  //current weight should be increased by eps

  //if (weight_diff) {
  //  if (without_suffix.GetBack() & 1) {
      //a --(-l)-> b  
      //b --( l)-> a
      //weight of l should be decreased by w
      //weight of prefix_end should be shifted by +eps
      this->ShiftWeight(prefix_end, -weight_diff);
  //  } else {
      //a --( l)-> b
      //b --(-l)-> a
      //weight of l should be decreased by w
      //weight of prefix_end should be shifted by +eps
  //    this->ShiftWeight(prefix_end, -weight_diff);
  //  }
  //}

  JoinVertices(prefix_end, suffix_begin);
  return true;
}

void FoldedGraph2::CompleteWith(Word r) {
  auto initial_vertex_count = edges_.size();
  for (size_t shift = 0; shift < r.size(); ++shift, r.CyclicLeftShift()) {
    for (auto vertex = root(); vertex < initial_vertex_count; ++vertex) {
      if (edges_[vertex].combined_with_) {
        continue;
      }

      PushCycle(r, vertex);
    }
  }
}

std::vector<unsigned int> FoldedGraph2::ComputeDistances(Vertex v) const {
  v = GetLastCombinedWith(v);
  std::vector<unsigned int> distance(edges_.size(), std::numeric_limits<unsigned int>::max());
  distance[v] = 0;
  
  std::deque<Vertex> q = {v};

  while (!q.empty()) {
    auto next = GetLastCombinedWith(q.front());
    q.pop_front();

    for(auto neighbor : vertex(next).edges_) {
      if (neighbor == kNullVertex) {
        continue;
      }

      neighbor = GetLastCombinedWith(neighbor);

      if (distance[next] + 1 < distance[neighbor]) {
        distance[neighbor] = distance[next] + 1;
        q.push_back(neighbor);
      }
    }
  }

  return distance;
}

template<typename Iter>
bool IsSortedAndUnique(Iter current, Iter end) {
  if (current == end) {
    return true;
  }
  auto next = std::next(current);
  while (next != end) {
    if (!(*current < *next)) {
      return false;
    }
    ++next;
    ++current;
  }
  return true;
}

std::vector<Word> FoldedGraph2::Harvest(size_t k, Vertex v1, Vertex v2, Weight weight) const {
  std::vector<Word> result;
  Harvest(k, v1, v2, weight, &result);
  return result;
}

void FoldedGraph2::Harvest(size_t k, Vertex v1, Vertex v2, Weight weight, std::vector<Word>* result) const {
  auto initial_result_length = result->size();
  v1 = GetLastCombinedWith(v1);
  v2 = GetLastCombinedWith(v2);
  auto v2_distances = this->ComputeDistances(v2);
  std::deque<std::tuple<Vertex, Word, Weight>> current_path = {std::make_tuple(v1, Word{ }, 0)};
  
  while (!current_path.empty()) {
    Vertex v;
    Word w;
    Weight c;
    std::tie(v, w, c) = current_path.front();
    current_path.pop_front();
    if (v == v2 && WeightMod(c - weight) == 0) {
      result->push_back(w);
    }

    auto& edges = vertex(v);

    for(auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (!w.Empty() && Inverse(label) == w.GetBack()) {
        continue;
      }
      Vertex n = GetLastCombinedWith(edges.endpoint(label));
      if (n == kNullVertex) {
        continue;
      }
      if (v2_distances[n] + w.size() < k) {
        Word next_word = w;
        next_word.PushBack(label);
        current_path.emplace_back(n, next_word, c + edges.weight(label));
      }
    }
  }

  assert(IsSortedAndUnique(result->begin() + initial_result_length, result->end()));
}

std::vector<Word> FoldedGraph2::Harvest(size_t k, Weight w) const {
  std::vector<Word> result;
  for(auto v = root(); v < edges_.size(); ++v) {
    if(vertex(v).combined_with_) {
      continue;
    }

    auto current_result_size = result.size();

    Harvest(k, v, v, w, &result);
    std::inplace_merge(result.begin(), result.begin() + current_result_size, result.end());
  }
  auto unique_end = std::unique(result.begin(), result.end());
  result.erase(unique_end, result.end());
  return result;
}

void FoldedGraph2::PrintAsDot(std::ostream* out) const {
  //\n  node [shape = point];
  (*out) << "digraph {\n  edge [arrowtail=dot,arrowhead=vee,arrowsize=0.5];\n  1 [color=blue];\n";
  int i = -1;
  for (auto&& vertex : edges_) {
    ++i;
    if (i == 0) {
      continue;
    }
    if (vertex.combined_with_) {
      continue;
    }
    for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (vertex.edges_[label]) {
        (*out) << i << " -> " << GetLastCombinedWith(vertex.edges_[label]) 
          << " [label=\"" 
          << (label % 2 ? "-" : "") << static_cast<char>('x' + static_cast<char>(label / 2)) << ";"
          << vertex.weights_[label] << "\"];\n";
      }
    }
  }

  (*out) << "}";
}

void FoldedGraph2::ShiftWeight(Vertex v, Weight shift) {
  if (shift == 0) {
    return;
  }

  auto& vertex = edges_[GetLastCombinedWith(v)];

  for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
    if (vertex.edges_[label] == kNullVertex) {
      continue;
    }

    if (Equal(vertex.edges_[label], v)) {
      assert(Equal(edges_[GetLastCombinedWith(vertex.edges_[label])].endpoint(Inverse(label)), v));
      continue;
    }

    auto& inverse_edge_weight = edges_[GetLastCombinedWith(vertex.edges_[label])].weights_[Inverse(label)];
    vertex.weights_[label] += shift;
    inverse_edge_weight -= shift;
  }
}



} //namespace crag
