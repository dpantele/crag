#include <deque>
#include <limits>

#include "folded_graph2.h"

namespace crag {

using Vertex = FoldedGraph2::Vertex;
using VertexEdges = FoldedGraph2::VertexEdges;
using Label = FoldedGraph2::Label;
using Word = FoldedGraph2::Word;

//! Transform label to the edges array index
inline size_t GetEdgesIndex(Label l) {
  assert(abs(l) > 0 && abs(l) <= FoldedGraph2::kAlphabetSize);
  if (l < 0) {
    // -1 -> 1, -2 -> 3, -3 -> 5, -4 -> 7
    return l * (-2) - 1;
  } else {
    // 1 -> 0, 2 -> 2, 3 -> 4, 4 -> 6 
    return ((l - 1) * 2);
  }
}

Vertex VertexEdges::endpoint(Label l) const {
  return edges_[GetEdgesIndex(l)];
}

Vertex& VertexEdges::endpoint(Label l) {
  return edges_[GetEdgesIndex(l)];
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadWord(const Word& w, size_t length_limit, Vertex s) const {
  if (w.empty()) {
    return std::make_tuple(s, 0u);
  }

  assert(s);

  if (length_limit > w.size()) {
    length_limit = w.size();
  }

  const auto end = w.begin() + length_limit;

  for (auto letter = w.begin(); letter != end; ++letter) {
    auto next = GetLastCombinedWith(vertex(s).endpoint(*letter));
    if (next) {
      s = next;
    } else {
      return std::make_tuple(s, letter - w.begin());
    }
  }

  return std::make_tuple(s, length_limit);
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadWord(const Word& w, Vertex s) const {
  return ReadWord(w, w.size(), s);
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadWord(const Word& w) const {
  return ReadWord(w, w.size(), root());
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadInverse(const Word& w, size_t length_limit, Vertex s) const {
  if (w.empty()) {
    return std::make_tuple(s, 0u);
  }

  assert(s);

  if (length_limit > w.size()) {
    length_limit = w.size();
  }

  const auto end = w.rbegin() + length_limit;

  for (auto letter = w.rbegin(); letter != end; ++letter) {
    auto next = vertex(s).endpoint(-*letter);
    if (next) {
      s = next;
    } else {
      return std::make_tuple(s, letter - w.rbegin());
    }
  }

  return std::make_tuple(s, length_limit);
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadInverse(const Word& w, Vertex s) const {
  return ReadInverse(w, w.size(), s);
}

std::tuple<Vertex, Word::size_type> FoldedGraph2::ReadInverse(const Word& w) const {
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

  assert(vertex(from).endpoint(l) == kNullVertex || Equal(vertex(from).endpoint(l), GetLastCombinedWith(to)));
  edges_[GetLastCombinedWith(from)].endpoint(l) = GetLastCombinedWith(to);

  assert(vertex(to).endpoint(-l) == kNullVertex || Equal(vertex(to).endpoint(-l), GetLastCombinedWith(from)));
  edges_[GetLastCombinedWith(to)].endpoint(-l) = GetLastCombinedWith(from);

  return to;
}

Vertex FoldedGraph2::PushWord(const Word& w, Vertex s) {
  if (w.empty()) {
    return s;
  }

  Word::size_type existing_length{};

  std::tie(s, existing_length) = ReadWord(w, s);
  for(auto letter = w.begin() + existing_length; letter != w.end(); ++letter) {
    s = AddEdge(*letter, s);
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

bool FoldedGraph2::PushCycle(const Word& w, Vertex s) {
  if (w.empty()) {
    return false;
  }

  Word::size_type existing_prefix_length{}, existing_suffix_length{};
  Vertex prefix_end{}, suffix_begin{};

  std::tie(prefix_end, existing_prefix_length) = ReadWord(w, s);
  

  auto after_prefix_letter = w.begin() + existing_prefix_length;
  Word::const_iterator begin_suffix_letter; 

  if (after_prefix_letter != w.end()) { //if we can't read the whole word
    std::tie(suffix_begin, existing_suffix_length) = ReadInverse(w, s);
    begin_suffix_letter = w.end() - existing_suffix_length;

    while (after_prefix_letter + 1 < begin_suffix_letter) {
      prefix_end = AddEdge(*after_prefix_letter, prefix_end);
      ++after_prefix_letter;
    }
    assert(after_prefix_letter != w.end());

    if (after_prefix_letter + 1 == begin_suffix_letter) {
      //if we added 1, -1, 1, we probably now have to join (1, -1) and (1)
      auto current_endpoint = vertex(prefix_end).endpoint(*after_prefix_letter);
      auto other_side = vertex(suffix_begin).endpoint(-*after_prefix_letter);
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
      
      AddEdge(*after_prefix_letter, prefix_end, suffix_begin);
      return true;
    }
  } else {
    if (prefix_end == s) {
      return false;
    } else {
      std::tie(suffix_begin, existing_suffix_length) = ReadInverse(w, s);
      begin_suffix_letter = w.end() - existing_suffix_length;
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


  auto before_common_length = begin_suffix_letter - w.begin();
  assert(before_common_length >= 0 && static_cast<unsigned int>(before_common_length) <= existing_prefix_length);

  std::tie(prefix_end, before_common_length) = ReadWord(w, before_common_length, s);
  assert(prefix_end != suffix_begin);

  JoinVertices(prefix_end, suffix_begin);
  return true;
}

Word CyclicShift(Word r) {
  if (r.size() <= 1) {
    return r;
  }

  auto front = r.front();
  auto first = r.begin();
  auto second = ++r.begin();
  for (; second != r.end(); ++first, ++second) {
    *first = *second;
  }
  r.back() = front;
  return r;
}

void FoldedGraph2::CompleteWith(Word r) {
  for (size_t shift = 0; shift < r.size(); ++shift, r = CyclicShift(std::move(r))) {
    for (auto vertex = root(); vertex <= edges_.size(); ++vertex) {
      if (edges_[vertex].combined_with_) {
        continue;
      }

      PushCycle(r, vertex);
    }
  }
}




} //namespace crag
