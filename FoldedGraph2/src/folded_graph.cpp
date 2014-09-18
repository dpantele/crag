#include <algorithm>
#include <deque>
#include <limits>
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
      assert(vertex(s).weight(w.GetFront()) == -vertex(next).weight(Inverse(w.GetFront())));
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
      if (!WeightMod(vertex.weight(label) + neigbour.weight(label))) {
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
  while (!vertex_to_join.empty()) {
    Vertex first, second, prev;
    Label through;

    std::tie(first, second, prev, through) = vertex_to_join.front();
    vertex_to_join.pop_front();

    first = GetLastCombinedWith(first);
    second = GetLastCombinedWith(second);
    prev = GetLastCombinedWith(prev);

    if (first == second) {
      for(auto label = 0u; label < kAlphabetSize * 2; ++label) {
        modulus_ = Gcd(std::abs(vertex(prev).weight(label) + vertex(first).weight(Inverse(label))), modulus_);
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

    assert(vertex(second).combined_with_ == kNullVertex);
    vertex(second).combined_with_ = first;

    if(prev != kNullVertex) {
      auto weight_diff = WeightMod(first_weights[Inverse(through)] + vertex(prev).weight(through));

      if (weight_diff) {
        for(size_t label = 0; label < kAlphabetSize * 2; ++label) {
          if (label ^ 1) {
            first_weights[label] -= weight_diff;
          } else {
            first_weights[label] += weight_diff;
          }
        }
      }
    }

    for (auto label = 0u; label < kAlphabetSize * 2; ++label) {
      if (first_edges[label] && second_edges[label]) {
        vertex_to_join.emplace_back(first_edges[label], second_edges[label], first, label);
      }
      if (!first_edges[label]) {
        first_edges[label] = second_edges[label] ? GetLastCombinedWith(second_edges[label]) : second_edges[label];
        assert(!first_weights[label]);
        first_weights[label] = second_weights[label];
      } else {
        first_edges[label] = GetLastCombinedWith(first_edges[label]);
      }
      auto weight_diff = WeightMod(first_weights[label] - second_weights[label]);
      assert(weight_diff == 0 || prev == kNullVertex || label != through);
      first_weights[label] -= weight_diff;
    }
  }

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
    //first ass a 'stem' if they coinside

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

  if (prefix_end > suffix_begin) {
    std::swap(prefix_end, suffix_begin);
  }

  auto weight_diff = WeightMod(current_weight - weight);

  if (weight_diff) {
    auto& weights = edges_[prefix_end].weights_;
    for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if(label ^ 1) {
        weights[label] -= weight_diff;
      } else {
        weights[label] += weight_diff;
      }
    }
  }

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

std::vector<Word> FoldedGraph2::Harvest(size_t k, Vertex v1, Vertex v2) const {
  v1 = GetLastCombinedWith(v1);
  v2 = GetLastCombinedWith(v2);
  auto v1_distances = this->ComputeDistances(v1);
  auto result = Harvest(k, v1, v2, v1_distances);
  std::sort(result.begin(), result.end());
  auto last = std::unique(std::make_move_iterator(result.begin()), std::make_move_iterator(result.end()));
  result.erase(last.base(), result.end());

  return result;
}

std::vector<Word> FoldedGraph2::Harvest(size_t k, Vertex v1, Vertex v2, const std::vector<unsigned int>& v1_distances) const {
  if (k == 0) {
    return {};
  }

  v1 = GetLastCombinedWith(v1);
  v2 = GetLastCombinedWith(v2);

  if (k < v1_distances[v2]) {
    return {};
  }

  std::vector<Word> result;

  for(Label label = 0; label < kAlphabetSize * 2; ++label) {
    auto endpoint = vertex(v2).endpoint(label);

    if (endpoint == kNullVertex) {
      continue;
    }

    if (endpoint == v1) {
      result.emplace_back(Word({Inverse(label)}));
    } 
    
    auto shorter_words = Harvest(k - 1, v1, endpoint, v1_distances);
    for (const auto& word : shorter_words) {
      if (word.GetBack() != label) { //don't add words with cancellations
        result.push_back(word);
        result.back().PushBack(Inverse(label));
      }
    }
  }

  if (result.size() == 1) {
    return result;
  }

  //now sort and remove equal
  return result;
}

} //namespace crag
