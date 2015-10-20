#pragma once
#ifndef CRAG_FOLDED_GRAPH2_H_
#define CRAG_FOLDED_GRAPH2_H_

#include <array>
#include <cassert>
#include <cstdint>
#include <deque>
#include <iostream>
#include <set>
#include <tuple>
#include <vector>
#include <stdint.h>

#include "compressed_word.h"

namespace crag {

//! Oriented labeled graph representing subgroup of FG[a, b]
class FoldedGraph2 {
public:
  //! Reference constant to be used as a size of the alphabet
  static const size_t kAlphabetSize = 2;

  //! Type representing vertices
  typedef unsigned int Vertex;

  //! Weight ofo edges
  typedef int64_t Weight;

  //! Type representing symbols and hence labels
  /**
   * We use the following presentation: \f$1 = a, -1 = a^{-1}, 2 = b, 3 = b^{-1}, \ldots\f$
   */
  typedef unsigned int Label;

  //! Return the label for the inverse of the symbol
  inline static Label Inverse(Label l) {
    return l ^ 1;
  }

  //! Get the ID of the root vertex. 
  static const Vertex root() {
    return kRootVertex;
  } 

  //! List of all edges from the vertex
  struct VertexEdges {
    VertexEdges()
      : edges_()
      , weights_()
    { }

    std::array<Vertex, kAlphabetSize * 2> edges_; //!< Id of the endpoint of the edge labeled by (i / 2)^(i % 2), NullVertex if no edge
    std::array<Weight, kAlphabetSize * 2> weights_; //!< Id of the endpoint of the edge labeled by (i / 2)^(i % 2), NullVertex if no edge
    mutable Vertex combined_with_ = kNullVertex; //!< When graph is folded, some vertices are sticked together. If the vertex is, this will be the id.
    Vertex current_combined_next_ = kNullVertex; //<! We need to trace which vertices were combined with the current one during vertex join

    ///@{
    //! Return the endpoint of the edge labeled by l
    Vertex endpoint(Label l) const;
    Vertex& endpoint(Label l);
    ///@}

    Weight weight(Label l) const;

  };

  //! Return the info for the vertex @ref i
  const VertexEdges& vertex(Vertex v) const;

  typedef CWord Word;

  ///@{
  //! Traces the word @ref w starting from the root.
  /**
  * @param w The word to be traced
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type, Weight> ReadWord(Word w) const;


  //! Traces the word @ref w starting from the vertex @ref s.
  /**
  * @param w The word to be traced
  * @param s The vertex to start from. The root by default
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type, Weight> ReadWord(Word w, Vertex s) const;

  //! Traces the word @ref w (but not longer than @ref length_limit) starting from the vertex @ref s.
  /**
  * @param w The word to be traced
  * @param length_limit The maximum length of prefix to be traced
  * @param s The vertex to start from. The root by default
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type, Weight> ReadWord(Word w, Word::size_type length_limit, Vertex s) const;

  //! Traces the inverse of the word @ref w starting from the root. 
  /**
  * @param w The word to be traced
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type, Weight> ReadInverse(Word w) const;

  //! Traces the inverse of the word @ref w starting from the vertex @ref s.
  /**
  * @param w The word to be traced
  * @param s The vertex to start from.
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type, Weight> ReadInverse(Word w, Vertex s) const;

  //! Traces the inverse of the word @ref w  (but not longer than @ref length_limit) starting from the vertex @ref s.
  /**
  * @param w The word to be traced
  * @param length_limit The maximum length of prefix to be traced
  * @param s The vertex to start from. The root by default
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type, Weight> ReadInverse(Word w, Word::size_type length_limit, Vertex s) const;

  ///@}


  //! Traces the word @ref w starting from the vertex @ref s, creating new vertices if required. 
  /**
  * @param w The word to be traced
  * @param s The vertex to start from. The root by default
  * @return The endpoint
  */
  Vertex PushWord(Word w, Vertex s = kRootVertex, Weight weight = 0);

  //! Traces the cycle @ref w starting from the vertex @ref s, creating new vertices if required. 
  /**
  * @param w The word to be traced
  * @param s The vertex to start from. The root by default
  * @return false if cycle already existed
  */
  bool PushCycle(Word w, Vertex s = 1, Weight weight = 0);

  //! For every vertex s and every cyclic permutation r' of the word r use pushCycle(r',s). 
  void CompleteWith(Word r);

  //! Tries to reduce the number edges with a non-trivial weight
  void Reweight();

  typedef std::deque<std::tuple<Vertex, Word, bool, Weight>> HarvestPath;
  //! Find all cycles starting at v of length up to k which can be read starting from the first_edge
  void Harvest(Word::size_type k, Weight w, Vertex origin_v, Vertex terminus_v, Label first_edge, std::vector<Word>* result) const;

  void Harvest(Word::size_type k, Weight w, Vertex origin_v, Vertex terminus_v, std::vector<Word>* result) const;


  //! Find all cycles starting at v of length up to k
  std::vector<Word> Harvest(Word::size_type k, Vertex origin, Vertex terminus, Weight w) const;

  //! Harvest all cycles of defined weight of length less than @ref k
  std::vector<Word> Harvest(Word::size_type k, Weight w = 0);

  Word RestoreConjugator(const Word& harvested_word) const {
    return GetPathFromRoot(RestoreHarvestVertex(harvested_word));
  }

  Vertex RestoreHarvestVertex(const Word& harvested_word) const;
  Word GetPathFromRoot(Vertex v) const;


  //! Returns true if vertices are equal
  bool Equal(Vertex v1, Vertex v2);

  //! Make @ref v1 and @ref v2 equal
  void JoinVertices(Vertex v1, Vertex v2);

  FoldedGraph2()
    : edges_(2)
  { }

  Weight modulus() const {
    return modulus_;
  }

  size_t size() const {
    return edges_.size();
  }

  void PrintAsDot(std::ostream* out) const;

  Weight WeightMod(Weight w) const {
    return modulus_ == 0 ? w : (((w % modulus_) + modulus_) % modulus_);
  }

  //! Get the number of edges of non-trivial weight
  uint64_t CountNontrivialEdges() const;

 private:
  //! The id of the root
  static const Vertex kRootVertex = 1;

  //! Null vertex
  static const Vertex kNullVertex = 0;

  std::vector<VertexEdges> edges_; //!< Main graph storage. Cell #i stores info about vertex #i
  Weight modulus_ = 0; //!< Something to make incomparable weigth comparable

  //! Return the last vertex in the list of combined vertices
  Vertex GetLastCombinedWith(Vertex v) const; 

  //! Make v2.combined_with_ = v1
  void Combine(Vertex v1, Vertex v2);

  //! Called in the end of JoinVertices to forget about all vertices which are not used anymore
  void ForgetCombined();

  //! Add new edge labeled by l from @ref from to @ref to (or a new vertex, if null)
  Vertex AddEdge(Label l, Vertex from, Vertex to = kNullVertex);

  //! Compute distances from the vertices to @ref v
  std::vector<unsigned int> ComputeDistances(Vertex v) const;

  //! Check that the weight of inverse edges are inverses of each other
  Vertex FindInconsistentWeights() const;

  //! All positively-labeled edges get +shift, all negtively-labeled get -shift
  void ShiftWeight(Vertex v, Weight shift);

};

} //namespace crag

#endif //CRAG_FOLDED_GRAPH2_H_