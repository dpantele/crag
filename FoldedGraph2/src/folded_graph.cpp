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

Vertex VertexEdges::endpoint(Label l) const {
  return edges_[l];
}

Vertex& VertexEdges::endpoint(Label l) {
  return edges_[l];
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadWord(Word w, Word::size_type length_limit, Vertex s) const {
  if (w.Empty()) {
    return std::make_tuple(s, 0u);
  }

  assert(s);

  s = GetLastCombinedWith(s);

  if (length_limit > w.size()) {
    length_limit = w.size();
  }

  for (auto i = 0u; i < length_limit; ++i) {
    auto next = GetLastCombinedWith(vertex(s).endpoint(w.GetFront()));
    if (next) {
      s = next;
    } else {
      return std::make_tuple(s, i);
    }
    w.PopFront();
  }

  return std::make_tuple(s, length_limit);
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadWord(Word w, Vertex s) const {
  return ReadWord(w, w.size(), s);
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadWord(Word w) const {
  return ReadWord(w, w.size(), root());
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadInverse(Word w, Word::size_type length_limit, Vertex s) const {
  if (w.Empty()) {
    return std::make_tuple(s, 0u);
  }

  assert(s);

  s = GetLastCombinedWith(s);

  if (length_limit > w.size()) {
    length_limit = w.size();
  }

  for (auto i = 0u; i < length_limit; ++i) {
    auto next = GetLastCombinedWith(vertex(s).endpoint(Inverse(w.GetBack())));
    if (next) {
      s = next;
    } else {
      return std::make_tuple(s, i);
    }
    w.PopBack();
  }

  return std::make_tuple(s, length_limit);
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadInverse(Word w, Vertex s) const {
  return ReadInverse(w, w.size(), s);
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadInverse(Word w) const {
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

Vertex FoldedGraph2::PushWord(Word w, Vertex s) {
  s = GetLastCombinedWith(s);
  if (w.Empty()) {
    return s;
  }

  Word::size_type existing_length{};

  std::tie(s, existing_length) = ReadWord(w, s);
  w.PopFront(existing_length);
  while (!w.Empty()) {
    s = AddEdge(w.GetFront(), s);
    w.PopFront();
  }

  assert(s > 1);
  return s;
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

  std::deque<std::pair<Vertex, Vertex>> vertex_to_join = {std::make_pair(v1, v2)};
  while (!vertex_to_join.empty()) {
    auto current = vertex_to_join.front();
    vertex_to_join.pop_front();

    current.first = GetLastCombinedWith(current.first);
    current.second = GetLastCombinedWith(current.second);

    if (current.first == current.second) {
      continue;
    }
    
    if (current.second < current.first) {
      std::swap(current.first, current.second);
    }

    auto& first_edges = edges_[current.first].edges_;
    auto& second_edges = edges_[current.second].edges_;

    assert(vertex(current.second).combined_with_ == kNullVertex);
    vertex(current.second).combined_with_ = GetLastCombinedWith(current.first);

    for (size_t label = 0; label < kAlphabetSize * 2; ++label) {
      if (first_edges[label] && second_edges[label] && !Equal(first_edges[label], second_edges[label])) {
        vertex_to_join.emplace_back(first_edges[label], second_edges[label]);
      }
      if (!first_edges[label]) {
        first_edges[label] = second_edges[label] ? GetLastCombinedWith(second_edges[label]) : second_edges[label];
      } else {
        first_edges[label] = GetLastCombinedWith(first_edges[label]);
      }
    }
  }
}

bool FoldedGraph2::PushCycle(Word w, Vertex s) {
  if (w.Empty()) {
    return false;
  }

  s = GetLastCombinedWith(s);

  Word::size_type existing_prefix_length{}, existing_suffix_length{};
  Vertex prefix_end{}, suffix_begin{};

  std::tie(prefix_end, existing_prefix_length) = ReadWord(w, s);

  auto after_prefix_word = w;
  after_prefix_word.PopFront(existing_prefix_length);

  if (!after_prefix_word.Empty()) { //if we can't read the whole word
    std::tie(suffix_begin, existing_suffix_length) = ReadInverse(w, s);

    while (existing_suffix_length + 1 < after_prefix_word.size()) {
      prefix_end = AddEdge(after_prefix_word.GetFront(), prefix_end);
      after_prefix_word.PopFront();
    }

    if (existing_suffix_length + 1 == after_prefix_word.size()) {
      //if we added 1, -1, 1, we probably now have to join (1, -1) and (1)
      auto current_endpoint = vertex(prefix_end).endpoint(after_prefix_word.GetFront());
      auto other_side = vertex(suffix_begin).endpoint(Inverse(after_prefix_word.GetFront()));
      if (current_endpoint == suffix_begin) {
        assert(other_side == prefix_end);
        return true;
      }
      
      if (current_endpoint != kNullVertex) {
        JoinVertices(current_endpoint, suffix_begin);
      } 

      if (other_side != kNullVertex) {
        JoinVertices(prefix_end, other_side);
      }
      
      AddEdge(after_prefix_word.GetFront(), prefix_end, suffix_begin);
      return true;
    }
  } else {
    if (prefix_end == s) {
      return false;
    } else {
      std::tie(suffix_begin, existing_suffix_length) = ReadInverse(w, s);
    }
  }

  /*
          suffix
        *--------*
  0  1  2  3  4  5
  ^________^
    prefix

  trace 0 1 
  
  */

  assert(w.size() >= existing_suffix_length);
  auto before_common_length = w.size() - existing_suffix_length;
  assert(before_common_length <= existing_prefix_length);

  std::tie(prefix_end, before_common_length) = ReadWord(w, before_common_length, s);
  assert(prefix_end != suffix_begin);

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
