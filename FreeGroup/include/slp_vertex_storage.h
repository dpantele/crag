/**
 * \file slp_vertex_storage.h
 *
 * Storage for vertices and some meta-data like hashes or flags of being visited
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_VERTEX_STORAGE_H_
#define CRAG_FREEGROUP_SLP_VERTEX_STORAGE_H_

#include <deque>

#include <boost/iterator/iterator_facade.hpp>

namespace crag {
namespace slp {

struct BasicVertex; //defined in slp_vertex.h
typedef int64_t TerminalSymbol;



namespace internal {

template <typename ItemType>
struct VertexMapStoredValue {
    bool inited_ = false; //Indicates if this value was recorded
    ItemType value_;
};

//template <typename ItemType, class VertexMapType>
//class VertexMapIterator : public boost::iterator_facade<
//    VertexMapIterator<ItemType, VertexMapType>,
//    ItemType,
//    std::random_access_iterator_tag
//> {
//  public:
//    VertexMapIterator()
//      : base_iterator_()
//    { }
//
//    VertexMapIterator(BaseIterator iterator)
//      : base_iterator_(iterator)
//    { }
//
//
//  private:
//    friend class boost::iterator_core_access;
//
//    bool equal(const VertexMapIterator<ItemType, VertexMapType>& other) const {
//      return this->base_iterator_ == other.base_iterator_;
//    }
//
//    void increment() {
//      ++this->base_iterator_;
//    }
//
//    void decrement() {
//      --this->base_iterator_;
//    }
//
//    reference dereference() const {
//      return base_iterator_->value_;
//    }
//
//    void advance(difference_type n) {
//      base_iterator_ += n;
//    }
//
//    void
//
//    typedef std::conditional<
//        std::is_const<ItemType>,
//        VertexMapType::BaseContainer::const_iterator,
//        VertexMapType::BaseContainer::iterator
//    > BaseIterator;
//
//    BaseIterator base_iterator_;
//};

} //namespace internal


// Map from Vertex into ItemType
// Can't be constructed without calling method from VertexStorage
// TODO: make this a real container
// iterator should be compatible with Vertex iterator from the vertex storage
// @tparam ItemType the type of values, must be constructible from Vertex
template <typename ItemType>
class VertexMap {
  public:
    //note that there are no TerminalVertex stored in VertexStorage
    //there could be some problems with this reference-return
    //only const-access, because we use explicit modification methods
    const ItemType& operator[](const Vertex& key) const noexcept;

    bool contains(const Vertex& key) const noexcept;
    ItemType& assign(const Vertex& key, ItemType value);

    //Semi-public typedef
    typedef std::deque<internal::VertexMapStoredValue<ItemType>> BaseContainer;
  private:

    BaseContainer data_;
    size_t first_vertex_index_; // the offset of the first vertex so that we don't have to
                                // create the mapping for the whole VertexStorage
                                // if we don't use it
};

} //namespace slp
} //namespace crag
#endif // CRAG_FREEGROUP_SLP_VERTEX_STORAGE_H_
