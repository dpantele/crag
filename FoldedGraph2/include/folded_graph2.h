#include <array>
#include <cassert>
#include <set>
#include <tuple>
#include <vector>




namespace crag {

//! Oriented labeled graph representing subgroup of FG[a, b]
class FoldedGraph2 {
public:
  //! Reference constant to be used as a size of the alphabet
  static const size_t kAlphabetSize = 2;

  //! Type representing vertices
  typedef unsigned int Vertex;

  //! Type representing symbols and hence labels
  typedef int Label;

  //! Get the ID of the root vertex. 
  static const Vertex root() {
    return kRootVertex;
  } 

  //! List of all edges from the vertex
  struct VertexEdges {
    VertexEdges()
      : edges_()
    { }

    std::array<Vertex, kAlphabetSize * 2> edges_; //!< Id of the endpoint of the edge labeled by (i / 2)^(i % 2), NullVertex if no edge
    mutable Vertex combined_with_ = kNullVertex; //!< When graph is folded, some vertices are sticked together. If the vertex is, this will be the id.

    ///@{
    //! Return the endpoint of the edge labeled by l
    Vertex endpoint(Label l) const;
    Vertex& endpoint(Label l);
    ///@}
  };


  //! Return the info for the vertex @ref i
  const VertexEdges& vertex(Vertex v) const;

  typedef std::vector<Label> Word;

  ///@{
  //! Traces the word @ref w starting from the root.
  /**
  * @param w The word to be traced
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type> ReadWord(const Word& w) const;


  //! Traces the word @ref w starting from the vertex @ref s.
  /**
  * @param w The word to be traced
  * @param s The vertex to start from. The root by default
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type> ReadWord(const Word& w, Vertex s) const;

  //! Traces the word @ref w (but not longer than @ref length_limit) starting from the vertex @ref s.
  /**
  * @param w The word to be traced
  * @param length_limit The maximum length of prefix to be traced
  * @param s The vertex to start from. The root by default
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type> ReadWord(const Word& w, size_t length_limit, Vertex s) const;

  //! Traces the inverse of the word @ref w starting from the root. 
  /**
  * @param w The word to be traced
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type> ReadInverse(const Word& w) const;

  //! Traces the inverse of the word @ref w starting from the vertex @ref s.
  /**
  * @param w The word to be traced
  * @param s The vertex to start from.
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type> ReadInverse(const Word& w, Vertex s) const;

  //! Traces the inverse of the word @ref w  (but not longer than @ref length_limit) starting from the vertex @ref s.
  /**
  * @param w The word to be traced
  * @param length_limit The maximum length of prefix to be traced
  * @param s The vertex to start from. The root by default
  * @return A tuple, the first element is the last vertex on the way, the second element is the length of the word traced
  */
  std::tuple<Vertex, Word::size_type> ReadInverse(const Word& w, size_t length_limit, Vertex s) const;

  ///@}


  //! Traces the word @ref w starting from the vertex @ref s, creating new vertices if required. 
  /**
  * @param w The word to be traced
  * @param s The vertex to start from. The root by default
  * @return The endpoint
  */
  Vertex PushWord(const Word& w, Vertex s = kRootVertex);

  //! Traces the cycle @ref w starting from the vertex @ref s, creating new vertices if required. 
  /**
  * @param w The word to be traced
  * @param s The vertex to start from. The root by default
  * @return false if cycle already existed
  */
  bool PushCycle(const Word& w, Vertex s = 1);

  //! For every vertex s and every cyclic permutation r' of the word r use pushCycle(r’,s). 
  void CompleteWith(Word r);

  //! Find all word of length up to k which can can be read from v1 to v2
  std::set<Word> Harvest(size_t k, Vertex v1, Vertex v2) const;

  //! Returns true if vertices are equal
  bool Equal(Vertex v1, Vertex v2);

  //! Make @ref v1 and @ref v2 equal
  void JoinVertices(Vertex v1, Vertex v2);

  FoldedGraph2()
    : edges_(2)
  { }

 private:
  //! The id of the root
  static const Vertex kRootVertex = 1;

  //! Null vertex
  static const Vertex kNullVertex = 0;

  std::vector<VertexEdges> edges_; //!< Main graph storage. Cell #i stores info about vertex #i

  //! Return the last vertex in the list of combined vertices
  Vertex GetLastCombinedWith(Vertex v) const; 

  //! Add new edge labeled by l from @ref from to @ref to (or a new vertex, if null)
  Vertex AddEdge(Label l, Vertex from, Vertex to = kNullVertex);

  //! Compute distances from the vertices to @ref v
  std::vector<unsigned int> ComputeDistances(Vertex v) const;

  std::set<Word> Harvest(size_t k, Vertex v1, Vertex v2, const std::vector<unsigned int>& v1_distances) const;

};

} //namespace crag