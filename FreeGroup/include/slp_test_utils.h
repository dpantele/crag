/*
 * slp_test_utils.h
 *
 *  Created on: Aug 23, 2013
 *      Author: dpantele
 */

#ifndef CRAG_FREEGROUP_SLP_TEST_UTILS_H_
#define CRAG_FREEGROUP_SLP_TEST_UTILS_H_

#include "slp.h"
#include "EndomorphismSLP.h"

namespace crag {
namespace slp {

template <size_t RANK, size_t COMPOSITION_SIZE, typename TerminalSymbol = int>
class RandomEndomorphismImageGenerator {
  public:
    static size_t current_seed_ = 0;
    size_t calls_left_;
    double inverter_probability_;

    RandomEndomorphismImageGenerator(size_t total_items_to_generate, double inverter_probability = 0.5)
      : calls_left_(total_items_to_generate)
      , inverter_probability_(inverter_probability)
    { }

    Vertex operator()(size_t image_num = 1) const {
      if (!calls_left_) {
        return Vertex();
      }

      --calls_left_;
      ++current_seed_;
      srand(current_seed_);
      UniformAutomorphismSLPGenerator<TerminalSymbol> generator(RANK, current_seed_);

      generator.set_inverters_probability(inverter_probability_);

      return EndomorphismSLP<TerminalSymbol>::composition(COMPOSITION_SIZE, generator).image(image_num);
    }
};

template <typename TerminalSymbol>
class EndomorphismFileReader {
  public:
    size_t current_image_;
    EndomorphismSLP<TerminalSymbol> endomorfism_;

    EndomorphismFileReader(const std::string& filename)
      : current_image_(0)
      , endomorfism_(EndomorphismSLP<TerminalSymbol>::load_from(&std::ifstream(filename.c_str())))
    { }

    Vertex operator()() const {
      if (current_image_ >= endomorfism_.max_non_trivial_image_symbol()) {
        return Vertex();
      }

      ++current_image_;

      return endomorfism_.image(current_image_);
    }
};

inline void print_tree_preorder_single(const Vertex& vertex, std::ostream* out) {
  std::unordered_set<slp::Vertex> visited;
  auto acceptor = [&visited] (const inspector::InspectorTask& task) {
    return visited.count(task.vertex) == 0 || visited.count(task.vertex.negate());//do not accept if vertex is mapped already
  };

  Inspector<inspector::Preorder, decltype(acceptor)> inspector(vertex, acceptor);

  while (!inspector.stopped()) {
    PrintTo(inspector.vertex(), out);
    (*out) << " << ";

    if (visited.count(inspector.vertex().left_child()) || visited.count(inspector.vertex().left_child().negate())) {
      out << "(p) ";
    }
    PrintTo(inspector.vertex().left_child(), out);
    out << " >> ";

    if (visited.count(inspector.vertex().right_child()) || visited.count(inspector.vertex().right_child().negate())) {
      out << "(p) ";
    }
    PrintTo(inspector.vertex().right_child(), out);

    (*out) << std::endl;
    visited.insert(inspector.vertex());

    ++inspector;
  }
}

inline std::string tree_preorder_single_string(const Vertex& vertex) {
  std::ostringstream out;
  print_tree_preorder_single(vertex, &out);

  return out.str();
}

struct TreePreorderSinglePrinter {
  const Vertex& vertex_;
  TreePreorderSinglePrinter(const Vertex& vertex)
    : vertex_(vertex)
  { }
};

inline std::ostream& operator<<(std::ostream& out, const TreePreorderSinglePrinter& vertex) {
  print_tree_preorder_single(vertex.vertex_, &out);

  return out;
}

struct TreeTextSerializer {
  const Vertex& vertex_;
  typedef std::unordered_set<Vertex::VertexSignedId> VisitedSet;
  VisitedSet* visited_prt_;
  TreeTextSerializer(const Vertex& vertex)
    : vertex_(vertex)
    , visited_prt_(&visited_local_)
  { }

  TreeTextSerializer(const Vertex& vertex, VisitedSet* visited)
    : vertex_(vertex)
    , visited_prt_(visited)
  { }
private:
  VisitedSet visited_local_; //if the user do not want to manage this object
};

inline std::ostream& operator<<(std::ostream& out, TreeTextSerializer& vertex_serializer) {
  auto acceptor = [&vertex_serializer] (const inspector::InspectorTask& task) {
    return vertex_serializer.visited_prt_->count(task.vertex.vertex_id()) == 0 ||
           vertex_serializer.visited_prt_->count(-task.vertex.vertex_id()) == 0;//do not accept if vertex is mapped already
  };

  Inspector<inspector::Postorder, decltype(acceptor)> inspector(vertex_serializer.vertex_, acceptor);

  if (vertex_serializer.vertex_.height() > 0) {
    out << vertex_serializer.vertex_.vertex_id() << " 1" << std::endl;
  } else {
    out << vertex_serializer.vertex_.vertex_id() << " 0" << std::endl;
  }

  while (!inspector.stopped()) {
    Vertex current_vertex = inspector.vertex();
    if (current_vertex.height() > 1) {
      if (current_vertex.vertex_id() < 0) {
        current_vertex = current_vertex.negate();
      }

      out << current_vertex.vertex_id() << ' '
          << current_vertex.left_child().vertex_id()  << ' ' << (current_vertex.left_child().height() > 1)  << ' '
          << current_vertex.right_child().vertex_id() << ' ' << (current_vertex.right_child().height() > 1) << std::endl;
    }

    vertex_serializer.visited_prt_->insert(inspector.vertex());
    ++inspector;
  }

  out << "0 0 0 0 0" << std::endl;

  return out;
}

struct TreeTextDeserializer {
  Vertex::VertexSignedId root_;
  bool root_is_nonterminal_;

  typedef std::map<Vertex::VertexSignedId, Vertex> VertexStorage;
  VertexStorage* vertices_prt_;
  TreeTextDeserializer()
    : root_(0)
    , root_is_nonterminal_(false)
    , vertices_prt_(&vertices_local_)
  { }

  TreeTextDeserializer(VertexStorage* vertices)
    : root_(0)
    , root_is_nonterminal_(false)
    , vertices_prt_(vertices)
  { }

  explicit operator bool() const {
    return root_ != 0;
  }

  Vertex get_nonterminal_vertex(Vertex::VertexSignedId vertex_id) const {
    auto vertex_iter = vertices_prt_->find(vertex_id);

    if (vertex_iter == vertices_prt_->end()) {
      vertex_iter = vertices_prt_->find(-vertex_id);
    }

    if (vertex_iter == vertices_prt_->end()) {
      throw std::out_of_range("Can't find the vertex with required id.");
    }

    if (vertex_iter->first == vertex_id) {
      return vertex_iter->second;
    } else {
      assert(vertex_iter->first == -vertex_id);
      return vertex_iter->second.negate();
    }
  }

  Vertex get_vertex(Vertex::VertexSignedId vertex_id, bool is_nonterminal) const {
    return !is_nonterminal ? TerminalVertex(vertex_id) : get_nonterminal_vertex(vertex_id);
  }

  Vertex root() const  {
    return get_nonterminal_vertex(root_);
  }

  void insert_new(Vertex::VertexSignedId archived_vertex_id, Vertex left, Vertex right) {
    assert(left && right);
    assert(archived_vertex_id > 0);

    vertices_prt_->emplace(archived_vertex_id, NonterminalVertex(left, right));
  }

private:
  VertexStorage vertices_local_; //if the user do not want to manage this object
};

inline std::istream& operator>>(std::istream& in, TreeTextDeserializer& vertex_deserializer) {
  if (!in) {
    return in;
  }

  in >> vertex_deserializer.root_ >> vertex_deserializer.root_is_nonterminal_;

  Vertex::VertexSignedId vertex_id;
  Vertex::VertexSignedId left_vertex_id;
  bool is_left_nonterminal;
  Vertex::VertexSignedId right_vertex_id;
  bool is_right_nonterminal;
  do {
    in >> vertex_id >> left_vertex_id >> is_left_nonterminal >> right_vertex_id >> is_right_nonterminal;

    if (vertex_id) {
      vertex_deserializer.insert_new(
          vertex_id,
          vertex_deserializer.get_vertex(left_vertex_id, is_left_nonterminal),
          vertex_deserializer.get_vertex(right_vertex_id, is_right_nonterminal)
      );
    }
  } while (vertex_id != 0);

  return in;
}


} //namespace slp
} //namespace crag

#endif /* CRAG_FREEGROUP_SLP_TEST_UTILS_H_ */
