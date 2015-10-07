#include <algorithm>
#include <deque>
#include <fstream>
#include <limits>
#include <map>
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
      assert(WeightMod(vertex(s).weight(Inverse(w.GetBack())) + vertex(next).weight(w.GetBack())) == 0);
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

void FoldedGraph2::Combine(Vertex v1, Vertex v2) {
  v1 = GetLastCombinedWith(v1);
  v2 = GetLastCombinedWith(v2);
  if (v1 == v2) {
    return;
  }

  if (v2 < v1) {
    std::swap(v1, v2);
  }

  assert(edges_[v2].combined_with_ == kNullVertex);
  edges_[v2].combined_with_ = v1;

  Vertex v1_combined_last = v1;
  Vertex next;
  while ((next = edges_[v1_combined_last].current_combined_next_)) {
    v1_combined_last = next;
  }

  assert(edges_[v1_combined_last].current_combined_next_ == kNullVertex);

  edges_[v1_combined_last].current_combined_next_ = v2;
}

void FoldedGraph2::ForgetCombined() {
  for(auto& v : edges_) {
    v.current_combined_next_ = kNullVertex;

    if (v.combined_with_) {
      //auto combined_with = GetLastCombinedWith(v.combined_with_);
      //for (auto label = 0u; label < kAlphabetSize * 2; ++label) {
        //if (v.edges_[label] != kNullVertex) {
        //  assert(false);
        //  if (vertex(combined_with).edges_[label] != kNullVertex) {
        //    assert(Equal(vertex(combined_with).edges_[label], v.edges_[label]));
        //    assert(WeightMod(edges_[combined_with].weights_[label] - v.weights_[label]) == 0);
        //  } else {
        //    auto end = GetLastCombinedWith(v.edges_[label]);
        //    edges_[combined_with].edges_[label] = end;
        //    edges_[combined_with].weights_[label] = v.weights_[label];
        //  }

        //  v.edges_[label] = kNullVertex;
        //  v.weights_[label] = 0;
        //}

      //  assert(v.edges_[label] == kNullVertex);
      //  assert(v.weights_[label] == 0);
      //}
    } else {
      for (auto& edge : v.edges_) {
        edge = GetLastCombinedWith(edge);
      }
    }
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
  std::set<std::tuple<Vertex, Label>> prohibited_endpoints;
  Vertex vertex_id = 0;
  for(auto&& vertex : edges_) {
    for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (vertex.edges_[label] == kNullVertex) {
        prohibited_endpoints.emplace(vertex_id, Inverse(label));
      }
    }
    ++vertex_id;
  }

  vertex_id = 0;
  for(auto&& vertex : edges_) {
    for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (!vertex.edges_[label]) {
        if (vertex.weights_[label]) {
          return vertex_id;
        }
        continue;
      }

      if (prohibited_endpoints.count(std::make_tuple(vertex.edges_[label], label))) {
        return vertex_id;
      }

      auto& neigbour = edges_[vertex.endpoint(label)];
      if (neigbour.endpoint(Inverse(label)) != vertex_id) {
        return vertex_id;
      }

      if (WeightMod(vertex.weight(label) + neigbour.weight(Inverse(label)))) {
        return vertex_id;
      }
    }
    ++vertex_id;
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

  assert(FindInconsistentWeights() == kNullVertex);

  std::deque<std::tuple<Vertex, Vertex, Label>> edges_to_join;
  typedef std::map<std::tuple<Vertex, Label>, Vertex> FoldedEdges;
  FoldedEdges folded_edges;

  auto JoinEndpoints = [this, &edges_to_join, &folded_edges](Vertex v1, Vertex v2) {
    v1 = GetLastCombinedWith(v1);
    v2 = GetLastCombinedWith(v2);
    if (v1 > v2) {
      std::swap(v2, v1);
    }
    if (v1 == v2) {
      return;
    }

    Combine(v1, v2);
  
    auto& v1_edges = edges_[v1];
    auto& v2_edges = edges_[v2];
  
    for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (v1_edges.edges_[label] && v2_edges.edges_[label]) {
        edges_to_join.emplace_back(v1, v2, label);
      } else if (v2_edges.edges_[label]) {
        v1_edges.edges_[label] = v2_edges.edges_[label];
        assert(v1_edges.weights_[label] == 0);
        v1_edges.weights_[label] = v2_edges.weights_[label];
        v2_edges.edges_[label] = kNullVertex;
        v2_edges.weights_[label] = 0;
        folded_edges.emplace(std::make_tuple(v2, label), v1);

        assert(edges_[v1_edges.edges_[label]].edges_[Inverse(label)] == v2);
        edges_[v1_edges.edges_[label]].edges_[Inverse(label)] = v1;
      }
    }
  };

  JoinEndpoints(v1, v2);

  auto GetUnfolded = [&folded_edges, this](Vertex v, Label l) -> Vertex {
    auto edge = std::make_tuple(v, l);
    auto folded = folded_edges.find(edge);
    while (folded != folded_edges.end()) {
      assert(edges_[std::get<0>(edge)].edges_[std::get<1>(edge)] == kNullVertex);
      std::get<0>(edge) = folded->second;
      folded = folded_edges.find(edge);
    }
    assert(edges_[std::get<0>(edge)].edges_[std::get<1>(edge)] != kNullVertex);
    return std::get<0>(edge);
  };

  //std::deque<std::tuple<Vertex, Vertex, Vertex, Label>> vertex_to_join = {std::make_tuple(v1, v2, kNullVertex, 0)};
  auto step = 0u;
  while (!edges_to_join.empty()) {
    ++step;

    Vertex v1, v2;
    Label through;

    std::tie(v1, v2, through) = edges_to_join.front();
//#define DEBUG_PRINT
#ifdef DEBUG_PRINT
    std::ofstream out("join-" + std::to_string(step) + "-1-" + std::to_string(v1) + "-" + std::to_string(v2) + "-" + std::to_string(through) + ".dot");
    this->PrintAsDot(&out);
    out.close();
#endif

    edges_to_join.pop_front();

    assert(v1 < v2);
    assert(Equal(v1, v2));

    v1 = GetUnfolded(v1, through);
    v2 = GetUnfolded(v2, through);

    if (v1 == v2) {
      continue;
    }

    if (v1 > v2) {
      std::swap(v1, v2);
    }

    auto v1_end = edges_[v1].edges_[through];
    auto v2_end = edges_[v2].edges_[through];

    assert(v1_end != kNullVertex);
    assert(v2_end != kNullVertex);
    
    auto weights_diff = edges_[v1].weights_[through] - edges_[v2].weights_[through];

    if (Equal(v1_end, v2_end)) {
      modulus_ = Gcd(std::abs(weights_diff), modulus_);
    } else {
      ShiftWeight(v1_end, weights_diff);
    }

#ifdef DEBUG_PRINT
    PrintAsDot(&std::ofstream("join-" + std::to_string(step) + "-2-" + std::to_string(v1) + "-" + std::to_string(v2) + "-" + std::to_string(through) + ".dot"));
#endif

    assert(WeightMod(edges_[v1].weights_[through] - edges_[v2].weights_[through]) == 0);
    assert(edges_[v2_end].edges_[Inverse(through)] == v2);
    assert(WeightMod(edges_[v2_end].weights_[Inverse(through)] - edges_[v1_end].weights_[Inverse(through)]) == 0);

    if (v1_end > v2_end) {
      std::swap(v1_end, v2_end);
    }
    edges_[v1].edges_[through] = v1_end;
    edges_[v2].edges_[through] = kNullVertex;
    edges_[v2].weights_[through] = 0;

    edges_[v1_end].edges_[Inverse(through)] = v1;
    edges_[v2_end].edges_[Inverse(through)] = kNullVertex;
    edges_[v2_end].weights_[Inverse(through)] = 0;

    folded_edges.emplace(std::make_tuple(v2, through), v1);
    folded_edges.emplace(std::make_tuple(v2_end, Inverse(through)), v1_end);

#ifdef DEBUG_PRINT
    PrintAsDot(&std::ofstream("join-" + std::to_string(step) + "-3-" + std::to_string(v1) + "-" + std::to_string(v2) + "-" + std::to_string(through) + ".dot"));
#endif

    if (!Equal(v1_end, v2_end)) {
      JoinEndpoints(v1_end, v2_end);
    }
  }

  ++step;
#ifdef DEBUG_PRINT
  std::ofstream out("join-" + std::to_string(step) + ".dot");
  this->PrintAsDot(&out);
  out.close();
#endif

  ForgetCombined();

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

#ifdef DEBUG_PRINT
  std::ofstream out("compl.dot");
  this->PrintAsDot(&out);
  out.close();
#endif

  auto weight_diff = WeightMod(weight - current_weight); // = eps
  this->ShiftWeight(prefix_end, -weight_diff);

  JoinVertices(prefix_end, suffix_begin);
  return true;
}

void FoldedGraph2::CompleteWith(Word r) {
  CompleteWith(std::move(r), edges_.size());
}

void FoldedGraph2::CompleteWith(Word r, size_t max_id) {
  for (size_t shift = 0; shift < r.size(); ++shift, r.CyclicLeftShift()) {
    for (auto vertex = root(); vertex < max_id; ++vertex) {
      if (edges_[vertex].combined_with_) {
        continue;
      }

      PushCycle(r, vertex);
    }
  }
}

struct CompleteInfo {
  bool was_completed_ = false;
  unsigned short nontrivial_path_min_length_ = 0;
};


void FoldedGraph2::BoundedCompleteWith(Word r, const Word::size_type max_path_length) {
  std::vector<Word> r_permutations;
  r_permutations.reserve(r.size());
  for (size_t shift = 0; shift < r.size(); ++shift, r.CyclicLeftShift()) {
    r_permutations.push_back(r);
  }

  std::vector<CompleteInfo> vertices_complete_info;
  vertices_complete_info.resize(edges_.size());

  auto CompleteVertex = [this, &r_permutations, &vertices_complete_info](Vertex v) {
    if (edges_[v].combined_with_) {
      return;
    }
    vertices_complete_info.resize(edges_.size());
    if (vertices_complete_info[v].was_completed_) {
      return;
    }

    vertices_complete_info[v].was_completed_ = true;

    for (auto&& r : r_permutations) {
      PushCycle(r, v);
    }
  };

  auto distances = DistanceToNontrivialEdges(max_path_length);
  
  Vertex v = 0;
  for (auto d : distances) {
    if (d != 0) {
      CompleteVertex(v);
    }
    ++v;
  }
}

std::vector<Word::size_type> FoldedGraph2::DistanceToNontrivialEdges(const Word::size_type max_path_length) const {
  std::vector<Word::size_type> vertex_path_length;
  vertex_path_length.resize(edges_.size());

  std::vector<Vertex> length_decreased;
  auto CheckLengthDecreased = [&length_decreased, &vertex_path_length](Vertex v, Word::size_type new_possible_length) {
    assert(v != 0);
    if (vertex_path_length[v] == 0 || vertex_path_length[v] > new_possible_length) {
      vertex_path_length[v] = new_possible_length;
      length_decreased.push_back(v);
    }
  };
  for (Vertex v = 1; v < edges_.size(); ++v) {
    if (edges_[v].combined_with_) {
      continue;
    }
    for (Label l = 0; l < 2 * kAlphabetSize; ++l) {
      if (edges_[v].weights_[l] != 0) {
        CheckLengthDecreased(v, 1);
        CheckLengthDecreased(GetLastCombinedWith(edges_[v].edges_[l]), 1);
      }
    }
  }

  for (auto id = 0u; id < length_decreased.size(); ++id) {
    auto v = length_decreased[id];
    auto current_length = vertex_path_length[v];
    if (current_length == max_path_length) {
      continue;
    }
    for (Label l = 0; l < 2 * kAlphabetSize; ++l) {
      if (edges_[v].edges_[l]) {
        CheckLengthDecreased(GetLastCombinedWith(edges_[v].edges_[l]), current_length + 1);
      }
    }
  }

  return vertex_path_length;
}


std::vector<Word::size_type> FoldedGraph2::ComputeDistances(Vertex v, Word::size_type max_distance) const {
  v = GetLastCombinedWith(v);
  std::vector<Word::size_type> distance(edges_.size());
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

      if (distance[next] + 1 < distance[neighbor] || (distance[neighbor] == 0 && neighbor != v)) {
        if (max_distance != 0 && distance[next] >= max_distance) {
          continue;
        }
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

#define UNUSED(x) ((void)x)

//void FoldedGraph2::Harvest(
//    Word::size_type k,
//    Vertex v2,
//    Weight weight,
//    FoldedGraph2::HarvestPath* current_path,
//    std::vector<Word>* result
//  ) const {
//
//  auto initial_result_length = result->size();
//  v2 = GetLastCombinedWith(v2);
//
//  std::vector<std::vector<CurrentPath>> vertex_paths(this->edges_.size());
//
//  vertex_paths[v2].emplace_back();
//
//  while (!current_path->empty()) {
//    Vertex v;
//    Word w;
//    bool is_primary;
//    Weight c;
//    std::tie(v, w, is_primary, c) = current_path->front();
//    current_path->pop_front();
//
//    if (v == v2 && (WeightMod(c - weight) == 0 || WeightMod(c + weight) == 0)) {
//      result->push_back(w);
//    }
//
//    auto& edges = vertex(v);
//
//    for(auto label = 0u; label < 2 * kAlphabetSize; ++label) {
//      if (!w.Empty() && Inverse(label) == w.GetBack()) {
//        continue;
//      }
//      Vertex n = GetLastCombinedWith(edges.endpoint(label));
//      if (n == kNullVertex) {
//        continue;
//      }
//      if (v2_distances[n] + w.size() < k && (n == v2 || v2_distances[n] != 0)) {
//        Word next_word = w;
//        next_word.PushBack(label);
//        current_path->emplace_back(n, next_word, c + edges.weight(label));
//      }
//    }
//  }
//
//  assert(IsSortedAndUnique(result->begin() + initial_result_length, result->end()));
//  UNUSED(initial_result_length);
//}

void FoldedGraph2::Harvest(
    Word::size_type k
    , Weight w
    , Vertex origin_v
    , Vertex terminus_v
    , Label first_edge
    , std::vector<Word> *result
) const {
  struct PathQueueElement {
    Vertex v;
    Word word;
    bool is_primary; //!< If this path starts from origin
    Weight weight;
  };
  std::deque<PathQueueElement> paths_queue;

  struct HarvestPath {
    Word word;
    Weight weight;
  };

  using VertexPaths = std::vector<HarvestPath>;
  std::vector<VertexPaths> cur_vertex_paths(edges_.size());

  auto PushPath = [&](PathQueueElement path) {
    if (!path.is_primary) {
      cur_vertex_paths[path.v].push_back(HarvestPath{path.word, path.weight});
    }
    paths_queue.push_back(std::move(path));
  };

  assert(edges_[origin_v].endpoint(first_edge) != kNullVertex);

  //add a single edge as a primary path
  PushPath(
      {
          GetLastCombinedWith(edges_[origin_v].endpoint(first_edge)),
          Word({first_edge}),
          true,
          edges_[terminus_v].weight(first_edge)
      }
  );

  //add trivial path to the final vertex
  cur_vertex_paths[terminus_v].push_back(HarvestPath{Word{}, 0});

  //also push all paths from the terminal vertex as non-primary
  if (k > 1) {
    for (auto terminus_label = 0u; terminus_label < 2 * kAlphabetSize; ++terminus_label) {
      if (edges_[terminus_v].edges_[terminus_label] == kNullVertex) {
        continue;
      }
      if (terminus_label == first_edge && origin_v == terminus_v) {
        //for cycles we don't need to add the primary path, we don't want cyclic cancellations
        continue;
      }
      PushPath(
          {
              GetLastCombinedWith(edges_[terminus_v].endpoint(terminus_label)),
              Word({terminus_label}),
              false,
              edges_[terminus_v].weight(terminus_label)
          }
      );
    }
  }

  if (paths_queue.size() != 1) {
    //otherwise there are no edges either from origin or terminus

    while (!paths_queue.empty()) {
      auto path = paths_queue.front();
      paths_queue.pop_front();
//      std::cerr << "v: " << path.v
//          << "\npath: " << path.word
//          << "\nis_primary: " << path.is_primary
//          << "\nweight: " << path.weight
//          << std::endl;

      if (path.is_primary) {
        //try to find new results
        for (const HarvestPath& vertex_path : cur_vertex_paths[path.v]) {
          //no cancellations
//          std::cerr << "t_path: " << vertex_path.word
//              << "\nt_weight: " << vertex_path.weight
//              << std::endl;

          if (!vertex_path.word.Empty() && vertex_path.word.GetBack() == path.word.GetBack()) {
            continue;
          }
          auto total_path_weight = path.weight - vertex_path.weight;
          if (WeightMod(total_path_weight - w) == 0 || WeightMod(total_path_weight + w) == 0) {
            //ok, weight is compatible
            //add new result
            assert(path.word.size() + vertex_path.word.size() <= k);
            Word new_result(vertex_path.word);
            new_result.Invert();
            new_result.PushFront(path.word);

            if (new_result.GetBack() != Inverse(new_result.GetFront())) {
              result->push_back(new_result);
//            std::cerr << "New result: " << new_result << std::endl;
            }
          }
        }
      }

      //add new paths
      //for a non-primary path we extend them if they are less than k/2
      //for primary path - one step more for odd k
      if (path.word.size() < floor(k/2.) || (path.is_primary && path.word.size() < ceil(k/2.))) {
        for (auto next_label = 0u; next_label < 2 * kAlphabetSize; ++next_label) {
          if (edges_[path.v].endpoint(next_label) == kNullVertex) {
            continue;
          }

          if (next_label == Inverse(path.word.GetBack())) {
            continue;
          }

          auto next_word = path.word;
          next_word.PushBack(next_label);
          assert(next_word.size() > path.word.size());
          PushPath(
              {
                  GetLastCombinedWith(edges_[path.v].endpoint(next_label)),
                  next_word,
                  path.is_primary,
                  path.weight + edges_[path.v].weight(next_label)
              }
          );
        }
      }
    }
  } else {
    //except the case when the only edge is from origin to terminus
    if (edges_[origin_v].edges_[first_edge] == terminus_v
        && (   WeightMod(edges_[origin_v].weights_[first_edge] - w) == 0
            || WeightMod(edges_[origin_v].weights_[first_edge] + w) == 0)
        ) {
      //then just add this edge for case origin != terminus and all cycles otherwise
      auto next_power = Word();
      unsigned int i = (origin_v == terminus_v)
                       ? k
                       : 1;
      while(i-- > 0u) {
        next_power.PushBack(first_edge);
        result->emplace_back(next_power);
      }
    }
  }
}

std::vector<Word> FoldedGraph2::Harvest(Word::size_type k, Vertex origin, Vertex terminus, Weight w) const {
  std::vector<Word> result;
  if (origin == terminus && WeightMod(w) == 0) {
    result.push_back(Word{});
  }

  for(auto label = 0u; label < 2 * kAlphabetSize; ++label) {
    if (edges_[origin].edges_[label] == kNullVertex) {
      continue;
    }

    Harvest(k, w, origin, terminus, label, &result);
  }

  std::sort(result.begin(), result.end());
  auto unique_end = std::unique(result.begin(), result.end());
  result.erase(unique_end, result.end());

  return result;
}



std::vector<Word> FoldedGraph2::Harvest(Word::size_type k, Weight w) {
  std::vector<Word> result;

  if (WeightMod(w) == 0) {
    result.push_back(Word{ });
  }

  //to perform Harvest, we need a queue of current paths
  //and a list of paths to every vertex

  for(auto v = root(); v < edges_.size(); ++v) {
    if(edges_[v].combined_with_) {
      continue;
    }

    for(auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (edges_[v].edges_[label] == kNullVertex) {
        continue;
      }
      if (edges_[v].weights_[label] == 0 && WeightMod(w) != 0) {
        continue;
      }

      Harvest(k, w, v, v, label, &result);

      edges_[edges_[v].edges_[label]].edges_[Inverse(label)] = kNullVertex;
      edges_[edges_[v].edges_[label]].weights_[Inverse(label)] = 0;
      edges_[v].edges_[label] = kNullVertex;
      edges_[v].weights_[label] = 0;
    }
  }
  std::sort(result.begin(), result.end());
  auto unique_end = std::unique(result.begin(), result.end());
  result.erase(unique_end, result.end());

  return result;
}

Vertex FoldedGraph2::RestoreHarvestVertex(const Word& harvested_word) const {
  std::vector<Word> result;

  static const Weight w = 1;

  for(auto v = root(); v < edges_.size(); ++v) {
    if(edges_[v].combined_with_) {
      continue;
    }

    for(auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (edges_[v].edges_[label] == kNullVertex) {
        continue;
      }

      if (edges_[v].weights_[label] == 0) {
        continue;
      }

      Harvest(harvested_word.size(), w, v, v, label, &result);
      for (auto&& elem : result) {
        if (elem == harvested_word) {
          assert(std::get<2>(ReadWord(harvested_word, v)) == w);
          return v;
        }
      }
      result.clear();
    }
  }

  return kNullVertex;
}

Word FoldedGraph2::GetPathFromRoot(Vertex v) const {
  std::map<Vertex, Word> paths = {{root(), Word{ }}};
  std::deque<Vertex> to_do = {root()};

  while (!to_do.empty() && paths.count(v) == 0) {
    auto n = to_do.front();
    to_do.pop_front();

    auto current_path = paths[n];
    assert(!current_path.Empty());
    Label l = 0;
    for (auto&& edge : edges_[n].edges_) {
      if (edge) {
        auto next_path = current_path;
        next_path.PushBack(l);
        auto is_new = paths.emplace(edge, next_path);
        if (is_new.second) {
          to_do.push_back(edge);
        }
      }
      ++l;
    }
  }

  return paths[v];
}

void FoldedGraph2::PrintAsUdot(std::ostream* out) const {
  (*out) << "graph {\n node [shape=point];\n";
  Vertex i = -1;
  for (auto&& vertex : edges_) {
    ++i;
    if (i == 0) {
      continue;
    }

    if (vertex.combined_with_) {
      continue;
    }

    for (auto label = 0u; label < 2 * kAlphabetSize; label+=2) {
      if (vertex.edges_[label]) {
        (*out) << i << " -- " << vertex.edges_[label];
        if (vertex.weights_[label] != 0) {
          (*out) << " [color=red]";
        }

        (*out) << "\n";
      }
    }
  }

  (*out) << "}";
}

void FoldedGraph2::PrintAsDot(std::ostream* out) const {
  (*out) << "digraph {\n  edge [arrowtail=dot,arrowhead=vee,arrowsize=0.5];\n  1 [color=blue];\n";
  Vertex i = -1;
  for (auto&& vertex : edges_) {
    ++i;
    if (i == 0) {
      continue;
    }

    if (vertex.combined_with_) {
        (*out) << i << " -> " << vertex.combined_with_
          << " [color=blue];\n";
    }

    for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (vertex.edges_[label]) {
        (*out) << i << " -> " << vertex.edges_[label]
          << " [label=\"" 
          << (label % 2 ? "-" : "") << static_cast<char>('x' + static_cast<char>(label / 2)) << ";"
          << vertex.weights_[label];

        (*out) << "\"];\n";
      }
    }
  }

  (*out) << "}";
}

void FoldedGraph2::ShiftWeight(Vertex v, Weight shift) {
  if (WeightMod(shift) == 0) {
    return;
  }

  v = GetLastCombinedWith(v);
  while (v) {
    auto& vertex = edges_[v];

    for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (vertex.edges_[label] == kNullVertex) {
        continue;
      }

      vertex.weights_[label] += shift;
      edges_[vertex.edges_[label]].weights_[Inverse(label)] -= shift;
    }

    v = edges_[v].current_combined_next_;
  }
}

void FoldedGraph2::GrowHair() {
  auto current_vertex_count = edges_.size();
  for(auto vertex_id = 1u; vertex_id < current_vertex_count; ++vertex_id) {
    if (edges_[vertex_id].combined_with_) {
      continue;
    }
    for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (edges_[vertex_id].edges_[label] == kNullVertex) {
        AddEdge(label, vertex_id, kNullVertex);
      }
    }
  }
}

void FoldedGraph2::Reweight() {
  std::deque<Vertex> to_check;
  for(Vertex v = root(); v < edges_.size(); ++v) {
    if (!edges_[v].combined_with_) {
      to_check.push_back(v);
    }
  }

  std::vector<Weight> w;
  while (!to_check.empty()) {
    auto v = to_check.front();
    to_check.pop_front();

    w.clear();
    w.reserve(2 * kAlphabetSize);
    for (Label label = 0; label < 2 * kAlphabetSize; ++label) {
      if (edges_[v].edges_[label]) {
        w.push_back(edges_[v].weights_[label]);
      }
    }

    assert(!w.empty());
    Weight top_weight = w[0];

    if (w.size() == 2) {
      if (w[0] == w[1]) {
        top_weight = w[0];
      } else {
        continue;
      }
    } else if (w.size() == 3) {
      if (w[0] == w[1] || w[0] == w[2]) {
        top_weight = w[0];
      } else if (w[1] == w[2]) {
        top_weight = w[1];
      } else {
        continue;
      }
    } else if (w.size() == 4) {
      if (w[0] == w[1] && (w[1] == w[2] || w[1] == w[3])) {
        //aaba & aaab & aaaa
        top_weight = w[0];
      } else if (w[2] == w[3] && (w[0] == w[2] || w[1] == w[2])) {
        //abaa & baaa
        top_weight = w[2];
      } else {
        continue;
      }
    }
    ShiftWeight(v, -top_weight);
    if (top_weight != 0) {
      for (auto n : edges_[v].edges_) {
        if (n != kNullVertex) {
          to_check.push_back(n);
        }
      }
    }
  }

  
}

uint64_t FoldedGraph2::CountNontrivialEdges() const {
  uint64_t result = 0;
  for(auto& v : edges_) {
    if (v.combined_with_) {
      continue;
    }
    for (auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (v.weights_[label] != 0) {
        ++result;
      }
    }
  }

  return result;
}

} //namespace crag
