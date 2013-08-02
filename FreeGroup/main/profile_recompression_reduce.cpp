/*
 * profile_matching_new.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: dpantele
 */

#include <memory>
#include <functional>
#include <utility>
#include <vector>
#include <tuple>
#include <algorithm>
#include <chrono>
#include <array>

#include "slp.h"
#include "slp_recompression.h"
#include "EndomorphismSLP.h"
#include "Permutation.h"
#include "gmp_boost_pool_allocator.h"
#include "boost/functional/hash/extensions.hpp"

using crag::slp::Vertex;
using crag::slp::VertexWord;
using crag::slp::PreorderInspector;
using crag::slp::TerminalVertexTemplate;
typedef crag::slp::TerminalVertexTemplate<int> TerminalVertex;
using crag::slp::NonterminalVertex;
using crag::slp::MatchingTable;
using crag::UniformAutomorphismSLPGenerator;
using crag::EndomorphismSLP;
using crag::slp::recompression::reduce;

//namespace crag { namespace slp {
//std::string print_tree_preorder_single(const Vertex& vertex) {
//  std::ostringstream out;
//  std::unordered_map<slp::Vertex, bool> mapping;
//  mapper::SkipMappedVerticesAcceptor<std::unordered_map<slp::Vertex, bool>> acceptor(mapping);
//  Inspector<inspector::Preorder, mapper::SkipMappedVerticesAcceptor<std::unordered_map<slp::Vertex, bool>>> inspector(vertex, acceptor);
//  while (!inspector.stopped()) {
//    PrintTo(inspector.vertex(), &out);
//    out << " << ";
//    if (mapping.count(inspector.vertex().left_child()) || mapping.count(inspector.vertex().left_child().negate())) {
//      out << "(p) ";
//    }
//    PrintTo(inspector.vertex().left_child(), &out);
//    out << " >> ";
//    if (mapping.count(inspector.vertex().right_child()) || mapping.count(inspector.vertex().right_child().negate())) {
//      out << "(p) ";
//    }
//    PrintTo(inspector.vertex().right_child(), &out);
//    out << std::endl;
//    mapping.insert(std::make_pair(inspector.vertex(), true));
//
//    ++inspector;
//  }
//
//  return out.str();
//}
//
//}}

int main() {
  gmp_pool_setup();
  int REPEAT = 1;
  constexpr size_t RANK = 6;
  constexpr size_t ENDOMORPHISMS_NUMBER = 200;
  size_t seed = 112233;
  UniformAutomorphismSLPGenerator<int> generator(RANK, seed);
  auto begin = std::chrono::system_clock::now();
  int count = REPEAT;

  auto duration = begin - begin;

  while (--count >= 0) {
    auto image = EndomorphismSLP<int>::composition(ENDOMORPHISMS_NUMBER, generator).image(1);
    auto reduce_start = std::chrono::system_clock::now();
    Vertex reduced = crag::slp::recompression::reduce(image);
    duration += std::chrono::system_clock::now() - reduce_start;

    auto end = std::chrono::system_clock::now();
    std::cout << "Step duration: "<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
  }

  auto end = std::chrono::system_clock::now();
  std::cout << "Total: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()/REPEAT << std::endl;
  std::cout << "Reduce: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()/REPEAT << std::endl;
}


