/*
 * slp_vertex.h
 *
 *  Created on: Feb 10, 2013
 *      Author: dpantele
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_VERTEX_H_
#define CRAG_FREEGROUP_SLP_VERTEX_H_

#include <memory>
#include <iostream>
#include <cassert>
#include <array>
#include <vector>

#include <gmpxx.h>
typedef mpz_class LongInteger;

#include <boost/pool/pool_alloc.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace crag { namespace slp {

//! Type of the basic character
/**
 * Supposed to be fully compatible with int - all algorithms assumes this is int.
 */
typedef int TerminalSymbol;

class Vertex;

namespace internal {
//! Struct to be stored in the SLP container.
struct NonterminalNode {
    NonterminalNode* left_child_ = nullptr; //nullptr if left is terminal
    NonterminalNode* right_child_ = nullptr;

    constexpr static TerminalSymbol kPositive  = TerminalSymbol{} + 1;
    constexpr static TerminalSymbol kNegative  = -kPositive;

    TerminalSymbol left_terminal_ = TerminalSymbol{}; //TerminalSymbol{} if null vertex, kPositive or kNegative (depends on the sign) if left is non-terminal
    TerminalSymbol right_terminal_ = TerminalSymbol{};

    mutable LongInteger length_ = LongInteger{}; //mutable because we want length to be lazy
    mutable unsigned int height_ = 0; //the same reason
};

class VertexStorageInterface;

#ifndef DNDEGUG
//emulation of reference with existence-checking for debugging
class VertexStorageDescriptor {
  public:
    constexpr inline VertexStorageDescriptor();
    explicit inline VertexStorageDescriptor(VertexStorageInterface* storage);
    inline bool is_valid() const noexcept  { return !weak_ptr_.expired(); }

    inline VertexStorageInterface& operator*() const noexcept;
    inline VertexStorageInterface* operator->() const noexcept;
    inline bool operator==(const VertexStorageDescriptor& other) const noexcept;
    explicit operator bool() const noexcept { return bare_ptr_; }

    std::weak_ptr<VertexStorageInterface> weak_ptr_;
    VertexStorageInterface* bare_ptr_;
};
#else
typedef VertexStorageInterface* VertexStorageDescriptor;
#endif //DNDEBUG


} //namespace internal

//! Basic interface for vertex. Also represents empty vertex, use Vertex() as empty vertex
class Vertex {
  public:
    //! Consturct empty vertex. Note that this should be very rarely used, only to represent invalid vertices.
    constexpr Vertex()
      : storage_()
      , nonterminal_data_()
      , terminal_data_()
    { }


    //! Compare two vertices
    /*
     * Two vertices are equal only if they refer to the same NonterminalNode (or both has nullptr there if terminal) and have
     * have the same terminal symbol stored (for nonterminal vertices this means the same sign)
     *
     * Note that comparing two vertices obtained from the different containers is undefined behaviour.
     */
    bool operator==(const Vertex& other) const noexcept {
      assert(this->storage_ == other.storage_ || !*this || !other);
      return
          nonterminal_data_ == other.nonterminal_data_ &&
          terminal_data_ == other.terminal_data_;
    }

    bool operator!=(const Vertex& other) const noexcept {
      return !(*this == other);
    }

    //! Return vertex representing the reversed word
    inline Vertex negate() const noexcept;

    //! Get left child
    inline Vertex left_child() const noexcept;

    //! Get right child
    inline Vertex right_child() const noexcept;

    //! Return the length of the left child
    inline const LongInteger& split_point() const noexcept;

    //! Length of the word represented by this vertex.
    /**
     * Length is calculated on the first call of this function. So it should be safe to change children before that.
     */
    const LongInteger& length() const noexcept;

    //! Maximum length of the path from this vertex to the terminal + 1. 1 if terminal. 0 if null.
    /**
     * Height is calculated on the first call of this function. So it should be safe to change children before that.
     */
    unsigned int height() const noexcept;

    //! Use this function to test if the vertex is terminal
    inline bool is_terminal() const noexcept;

    //! Function to be used by std:hash<Vertex>
    inline size_t vertex_hash() const noexcept;

    //! Some order on vertices to store them in std::map
    inline bool operator<(const Vertex& other) const noexcept;

    //! True if vertex is not-null
    inline explicit operator bool() const noexcept;

    //! Function to be use in gtest's print. But can be called without it.
    void debug_print(::std::ostream* out) const;

    internal::VertexStorageDescriptor storage() { return storage_; }
    Vertex(internal::VertexStorageDescriptor storage, internal::NonterminalNode* nonterminal_data, TerminalSymbol terminal_data)
          : storage_(storage)
          , nonterminal_data_(nonterminal_data)
          , terminal_data_(terminal_data)
        { }

  protected:
    friend class internal::VertexStorageInterface;
    internal::VertexStorageDescriptor storage_{}; //to create new null vertices without explicitly indicating the storage

    internal::NonterminalNode* nonterminal_data_ = nullptr;
    TerminalSymbol terminal_data_ = TerminalSymbol{};

    static constexpr std::hash<internal::NonterminalNode*> nonterminal_ptr_hash = std::hash<internal::NonterminalNode*>();
    static constexpr std::hash<TerminalSymbol> terminal_symbol_hash = std::hash<TerminalSymbol>();

    static const LongInteger& LongZero();
    static const LongInteger& LongOne();

};

namespace internal {

class VertexStorageInterface : public std::enable_shared_from_this<VertexStorageInterface> {
  public:
    virtual ~VertexStorageInterface() {}

    virtual Vertex make_nonterminal(Vertex left, Vertex right) = 0;
    virtual Vertex make_terminal(TerminalSymbol terminal) = 0;
  protected:
    //Access to the privates of vertices for storages
    static NonterminalNode* nonterminal_data(Vertex vertex) { return vertex.nonterminal_data_; }
    static TerminalSymbol terminal_data(Vertex vertex) { return vertex.terminal_data_; }
    static VertexStorageDescriptor storage(Vertex vertex) { return vertex.storage_; }
};

}  // namespace internal

//this is not a full implementations - it has no constructors. Use appropriate storage, like BucketVertexStorage
class VertexStorage {
  public:
    VertexStorage() = delete;

    //Create non-terminal vertex in this container
    Vertex make_nonterminal(Vertex left, Vertex right) { return storage_->make_nonterminal(std::move(left), std::move(right)); }
    Vertex make_terminal(TerminalSymbol terminal) { return storage_->make_terminal(terminal); }

    VertexStorage(std::shared_ptr<internal::VertexStorageInterface> storage)
      : storage_(std::move(storage))
    { }
  protected:
    std::shared_ptr<internal::VertexStorageInterface> storage_ = nullptr;
};



//! Method for gtest printing
inline void PrintTo(const Vertex& vertex, ::std::ostream* os) {
  vertex.debug_print(os);
}

//! Terminal vertex in a SLP. Produces word of length 1.
class TerminalVertex : public Vertex {
  public:
    TerminalVertex() = delete;

    explicit TerminalVertex(TerminalSymbol terminal_symbol, VertexStorage* storage)
      : Vertex(storage->make_terminal(terminal_symbol))
    { }

    explicit TerminalVertex(const Vertex& vertex)
      : Vertex(vertex)
    {
      assert((height() == 0 && length() == 0) || (height() == 1 && length() == 1));
    }

    const TerminalSymbol& terminal_symbol() const {
      return terminal_data_;
    }

    operator TerminalSymbol() const {
      return terminal_symbol();
    }
};

inline ::std::ostream& operator << (::std::ostream& stream, const TerminalVertex& vertex) {
  return stream << vertex.terminal_symbol();
}

//! Non-terminal vertex in a SLP, represent rule A->BC
class NonterminalVertex : public Vertex {
  public:
    NonterminalVertex() = delete;

    NonterminalVertex(Vertex left, Vertex right)
      : Vertex(left.storage()->make_nonterminal(
          left,
          right
        ))
    {
    }
};


namespace internal {

class VertexBucket;

struct NonterminalNodeInBucket : public NonterminalNode {
    VertexBucket* bucket_ = nullptr;
};

//sizeof(VertexBucket) should be ~4Kb
class VertexBucket {
  public:
    constexpr static size_t kBucketLength = 0x2000ull / sizeof(NonterminalNodeInBucket);

    VertexBucket()
    { }

    NonterminalNodeInBucket* allocate() noexcept {
      assert(allocated_size_ < kBucketLength);
      data_[allocated_size_].bucket_ = this;
      return data_ + allocated_size_++;
    }

    const NonterminalNodeInBucket* begin() const noexcept {
      return data_;
    }

    NonterminalNodeInBucket* begin() noexcept {
      return data_;
    }

    bool is_full() const noexcept {
      return allocated_size_ == kBucketLength;
    }

  private:
    NonterminalNodeInBucket data_[kBucketLength] = {};
    size_t allocated_size_ = 0;
};

//Should be shared-referenced
class BucketVertexStorageImpl : public VertexStorageInterface {
  public:
    NonterminalNodeInBucket* allocate() {
      if (vertex_buckets_.empty() || vertex_buckets_.back()->is_full()) {
        vertex_buckets_.emplace_back(new VertexBucket{});
      }
      return vertex_buckets_.back()->allocate();
    }

    Vertex make_nonterminal(Vertex left, Vertex right) override {
      assert(left && right);
      assert(storage(left) == storage(right));
      assert(storage(left) == this);

      auto node = allocate();

      node->left_child_ = nonterminal_data(left);
      node->right_child_ = nonterminal_data(right);
      node->left_terminal_ = terminal_data(left);
      node->right_terminal_ = terminal_data(right);

      return Vertex(VertexStorageDescriptor(this), node, NonterminalNode::kPositive);
    }

    virtual Vertex make_terminal(TerminalSymbol terminal) override {
      return Vertex(VertexStorageDescriptor(this), nullptr, terminal);
    }
  private:
    std::vector<std::unique_ptr<VertexBucket>> vertex_buckets_;

};

}  // namespace internal

class BucketVertexStorage : public VertexStorage {
  public:
    BucketVertexStorage()
      : VertexStorage(std::make_shared<internal::BucketVertexStorageImpl>())
    { }
};

}}//namespace crag::slp

namespace std {

//! Definition of the hash for std::pair
template<typename TFirst, typename TSecond>
struct hash<pair<TFirst, TSecond>> {
private:
  constexpr static hash<TFirst> first_hash_ = hash<TFirst>();
  constexpr static hash<TSecond> second_hash_ = hash<TSecond>();
public:
  size_t operator()(const pair<TFirst, TSecond>& obj) const {
    size_t first_hash_value = first_hash_(obj.first);
    //Taken from boost/functional/hash
    return second_hash_(obj.second) + 0x9e3779b9 + (first_hash_value << 6) + (first_hash_value >> 2);
  }
};

template<typename TFirst, typename TSecond>
constexpr hash<TFirst> hash<pair<TFirst, TSecond>>::first_hash_;
template<typename TFirst, typename TSecond>
constexpr hash<TSecond> hash<pair<TFirst, TSecond>>::second_hash_;

//! Definition of the hash for SignedVertex
template<>
struct hash<crag::slp::Vertex> {
  public:
    size_t operator()(const crag::slp::Vertex& vertex) const {
      return vertex.vertex_hash();
    }
};

} //namespace std

//Definitions of inline methods
#include "slp_vertex-inl.h"
#endif /* CRAG_FREEGROUP_SLP_VERTEX_H_ */
