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

#define UNUSED(x) ((void)x)



//! Calls functor on each equivalent range of common elements from first and second
template<typename RandomIterator, typename Compare, typename F>
void MergeRanges(
      RandomIterator first_begin
    , RandomIterator first_end
    , RandomIterator second_begin
    , RandomIterator second_end
    , Compare comp
    , F f
) {
  while (first_begin != first_end && second_begin != second_end) {
    if (comp(*first_begin, *second_begin)) {
      first_begin = std::lower_bound(first_begin, first_end, *second_begin, comp);
    } else if (comp(*second_begin, *first_begin)) {
      second_begin = std::lower_bound(second_begin, second_end, *first_begin, comp);
    } else {
      auto first_equal = std::equal_range(first_begin, first_end, *first_begin, comp);
      auto second_equal = std::equal_range(second_begin, second_end, *second_begin, comp);

      f(first_equal.first, first_equal.second, second_equal.first, second_equal.second);

      first_begin = first_equal.second;
      second_begin = second_equal.second;
    }
  }
}

void FoldedGraph2::Harvest(
    Word::size_type k
    , Weight w
    , Vertex origin_v
    , Vertex terminus_v
    , std::vector<Word> *result
) const {

  struct Path {
    Vertex t_;
    Weight w_;
    Word word_;
  };

  auto HarvestPaths = [this](size_t max_length, Vertex v0, auto path_action) {
    std::deque<Path> active_paths;
    active_paths.push_back(Path{v0, 0, Word()});

    while (!active_paths.empty()) {
      Path current_path = std::move(active_paths.front());
      active_paths.pop_front();
      path_action(current_path);

      if (current_path.word_.size() == max_length) {
        continue;
      }

      for (auto next_letter = 0u; next_letter < 2 * kAlphabetSize; ++next_letter) {
        if (edges_[current_path.t_].endpoint(next_letter) == kNullVertex) {
          continue;
        }

        if (!current_path.word_.Empty() && next_letter == Inverse(current_path.word_.GetBack())) {
          continue;
        }

        auto next_word = current_path.word_;
        next_word.PushBack(next_letter);
        assert(next_word.size() > current_path.word_.size());

        active_paths.push_back(Path{
            GetLastCombinedWith(edges_[current_path.t_].endpoint(next_letter)),
            WeightMod(current_path.w_ + edges_[current_path.t_].weight(next_letter)),
            next_word
        });
      }
    }
  };

  //these should be always of length ceil(k/2.) and start from origin_v
  std::vector<Path> prefixes;

  //these should be of length not more than floor(k/2.) and start from terminus_v
  std::vector<Path> suffixes;

  assert(static_cast<size_t>(ceil(k/2.)) == ceil(k/2.));
  assert(static_cast<size_t>(floor(k/2.)) == floor(k/2.));
  if (origin_v != terminus_v) {
    HarvestPaths(static_cast<size_t>(ceil(k/2.)), origin_v, [&](const Path& prefix) {
      //if the path hits terminus, add that to the result
      if (prefix.t_ == terminus_v && WeightMod(prefix.w_ - w) == 0) {
        result->push_back(prefix.word_);
      }
      if (prefix.word_.size() == ceil(k/2.)) {
        prefixes.push_back(prefix);
      }
    });

    HarvestPaths(static_cast<size_t>(floor(k/2.)), terminus_v, [&](const Path& suffix) {
      if (!suffix.word_.Empty()) {
        //suffix may not be empty - all such paths are added while harvesting prefixes
        suffixes.push_back(suffix);
      }
    });
  } else {
    //try to do the same this in a single HarvestPaths
    HarvestPaths(static_cast<size_t>(ceil(k/2.)), origin_v, [&](const Path& path) {
      //if the path hits terminus, add that to the result
      if (path.t_ == terminus_v
          && WeightMod(path.w_ - w) == 0
          && (path.word_.Empty() || path.word_.GetFront() != Inverse(path.word_.GetBack()))
      ) {
        result->push_back(path.word_);
      }
      if (path.word_.size() == ceil(k/2.)) {
        prefixes.push_back(path);
      }
      if (path.word_.size() <= floor(k/2.) && !path.word_.Empty()) {
        suffixes.push_back(path);
      }
    });
  }

  //we will merge prefixes and suffixes so that they have the same terminus
  //and also we will look for the path where prefix.w_ + suffix_.w_ is equal
  //to @param w
  //so we sort by a) terminus b) weight and c) words
  // (to have them sorted afterwards)
  std::sort(prefixes.begin(), prefixes.end(), [](const Path& p1, const Path& p2) {
    if(p1.t_ != p2.t_) {
      return p1.t_ < p2.t_;
    }
    if (p1.w_ != p2.w_) {
      return p1.w_ < p2.w_;
    }
    return p1.word_ < p2.word_;
  });
  std::sort(suffixes.begin(), suffixes.end(), [](const Path& p1, const Path& p2) {
    if(p1.t_ != p2.t_) {
      return p1.t_ < p2.t_;
    }
    if (p1.w_ != p2.w_) {
      return p1.w_ < p2.w_;
    }
    return p1.word_ < p2.word_;
  });


  //routine which combines prefixes and suffixes of the appropriate weights
  auto ConcatenatePrefixesSuffixes =
      [this, w]
      (auto prefixes_begin, auto prefixes_end, auto suffixes_begin, auto suffixes_end, auto PushPath) {
    while (prefixes_begin != prefixes_end && suffixes_begin != suffixes_end) {
      auto current_prefix_weight = prefixes_begin->w_;
      //suffixes will be appended inversed, so their weight must be inversed now
      auto needed_suffix_weight = this->WeightMod(-(w - current_prefix_weight));
      auto suffix_weight_range = std::equal_range(
          suffixes_begin
          , suffixes_end
          , Path{kNullVertex, needed_suffix_weight, Word{}}
          , [](const Path& p1, const Path& p2) { return p1.w_ < p2.w_; }
      );
      if (suffix_weight_range.first == suffix_weight_range.second) {
        prefixes_begin = std::upper_bound(
            prefixes_begin
            , prefixes_end
            , *prefixes_begin
            , [](const Path& p1, const Path& p2) { return p1.w_ < p2.w_; }
        );
      } else {
        auto prefix_weight_range = std::equal_range(
            prefixes_begin
            , prefixes_end
            , *prefixes_begin
            , [](const Path& p1, const Path& p2) { return p1.w_ < p2.w_; }
        );

        for (auto prefix = prefix_weight_range.first; prefix != prefix_weight_range.second; ++prefix) {
          for (auto suffix = suffix_weight_range.first; suffix != suffix_weight_range.second; ++suffix) {
            if (   suffix->word_.Empty()
                || prefix->word_.Empty()
                || suffix->word_.GetBack() != prefix->word_.GetBack()
              ) {
              Word new_word = suffix->word_;
              new_word.Invert();
              new_word.PushFront(prefix->word_);
              PushPath(std::move(new_word));
            }
          }
        }
        if (this->WeightMod(-1) == -1) {
          // in this case weights of suffixes are sorted as well, and we may not consider
          // any suffixes if weight less that the weight of suffix_weight_range.second
          suffixes_begin = suffix_weight_range.second;
        }
        prefixes_begin = prefix_weight_range.second;
      }
    }
  };

  if (origin_v != terminus_v) {
    MergeRanges(
        prefixes.begin()
        , prefixes.end()
        , suffixes.begin()
        , suffixes.end()
        , [](const Path& p1, const Path& p2) {
          return p1.t_ < p2.t_;
        }
        , [&](auto prefixes_begin, auto prefixes_end, auto suffixes_begin, auto suffixes_end) {
          ConcatenatePrefixesSuffixes(
              prefixes_begin
              , prefixes_end
              , suffixes_begin
              , suffixes_end
              , [&](Word new_word) {
                result->push_back(new_word);
              }
          );
        }
    );
  } else {
    //for cycles case, we don't allow cyclic reductions
    MergeRanges(
        prefixes.begin(), prefixes.end(), suffixes.begin(), suffixes.end(), [](const Path &p1, const Path &p2) {
          return p1.t_ < p2.t_;
        }, [&](auto prefixes_begin, auto prefixes_end, auto suffixes_begin, auto suffixes_end) {
          ConcatenatePrefixesSuffixes(
              prefixes_begin, prefixes_end, suffixes_begin, suffixes_end, [&](Word new_word) {
                if (new_word.GetFront() != Inverse(new_word.GetBack())) {
                  result->push_back(std::move(new_word));
                }
              }
          );
        }
    );
  }
}


std::vector<Word> FoldedGraph2::Harvest(Word::size_type k, Vertex origin, Vertex terminus, Weight w) const {
  std::vector<Word> result;

  Harvest(k, w, origin, terminus, &result);

  std::sort(result.begin(), result.end());
  assert(IsSortedAndUnique(result.begin(), result.end()));

  return result;
}

void FoldedGraph2::Harvest(
    Word::size_type k,
    Weight w,
    Vertex origin_v,
    Vertex terminus_v,
    Label first_edge,
    std::vector<Word> *result
) const {
  origin_v = GetLastCombinedWith(origin_v);
  terminus_v = GetLastCombinedWith(terminus_v);
  auto second_v = edges_[origin_v].endpoint(first_edge);

  if (second_v == kNullVertex) {
    return;
  }

  if (k == 0) {
    return;
  }

  second_v = GetLastCombinedWith(second_v);

  std::vector<Word> this_result;
  Harvest(static_cast<CWord::size_type>(k - 1),
          WeightMod(w - edges_[origin_v].weight(first_edge)),
          second_v,
          terminus_v,
          &this_result);

  std::sort(result->begin(), result->end());

  assert(IsSortedAndUnique(this_result.begin(), this_result.end()));
  for (auto&& word : this_result) {
    if (word.GetBack() != Inverse(first_edge)) {
      word.PushFront(first_edge);
      result->push_back(word);
    }
  }
}



std::vector<Word> FoldedGraph2::Harvest(Word::size_type k, Weight w) {
  std::vector<Word> result;

  if (WeightMod(w) == 0) {
    result.push_back(Word{ });
  }

  for(auto v = root(); v < edges_.size(); ++v) {
    if(edges_[v].combined_with_) {
      continue;
    }

    //optimization: don't harvest cycles which start from
    //non-zero-weight edges for non-zero w

    if (WeightMod(w) != 0) {
      bool has_non_trivial = false;
      for(auto label = 0u; label < 2 * kAlphabetSize; ++label) {
        if (edges_[v].edges_[label] == kNullVertex) {
          continue;
        }
        if (WeightMod(edges_[v].weights_[label]) != 0) {
          has_non_trivial = true;
          break;
        }
      }
      if (!has_non_trivial) {
        continue;
      }
    }

    Harvest(k, w, v, v, &result);

    //just get rid of all edges at this vertex
    for(auto label = 0u; label < 2 * kAlphabetSize; ++label) {
      if (edges_[v].edges_[label] == kNullVertex) {
        continue;
      }
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

    Harvest(harvested_word.size(), w, v, v, &result);
    for (auto&& elem : result) {
      if (elem == harvested_word) {
        assert(std::get<2>(ReadWord(harvested_word, v)) == w);
        return v;
      }
    }
    result.clear();
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
