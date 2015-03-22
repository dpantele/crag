/**
* @file
* @brief Definition of DisjointSubset
*/
#include <cassert>
#include <cstddef>
#include <utility>

#pragma once
#ifndef _CRAG_DISJOINTSUBSET_H_
#define _CRAG_DISJOINTSUBSET_H_

//! Element of the disjoint-set forest
template <typename Label>
class DisjointSubset {
 public:

  //! Construct a single-element subset
  DisjointSubset(Label l)
    : size_(1)
    , parent_(this)
    , label_(std::move(l))
  { }

  //! Construct an invalid element
  template<typename NullPtrT>
  DisjointSubset(NullPtrT,
    typename std::enable_if<
      std::is_same<NullPtrT, std::nullptr_t>::value
        && std::is_default_constructible<Label>::value
    >::type* = 0)
    : size_(0)
    , parent_(nullptr)
    , label_()
  { }

  DisjointSubset(const DisjointSubset&) = delete;
  DisjointSubset(DisjointSubset&&) = delete;

  //! Get the label of the root of the subset
  const Label& root() const {
    return FindRoot()->label_;
  }

  //! Get the label of the root of the subset
  Label& root() {
    return FindRoot()->label_;
  }

  //! Get the size of the subset
  size_t size() const {
    return FindRoot()->size_;
  }

  //! Check if this node is the root of the tree
  bool IsRoot() const {
    return parent_ == this;
  }

  //! Check if subset is invalid
  explicit operator bool() const {
    return parent_ != nullptr;
  }

  //! Merge this subset with another one, return the label of the new root
  Label Merge(DisjointSubset<Label>* other) {
    if (!*this) {
      return other->Merge(this);
    }

    if (!*other) {
      return this->label_;
    }

    if (!this->IsRoot()) {
      return FindRoot()->Merge(other);
    }

    if (!other->IsRoot()) {
      return Merge(other->FindRoot());
    }

    if (this->size_ < other->size_) {
      return other->Merge(this);
    }

    if (this == other) {
      return this->label_;
    }

    this->size_ += other->size_;
    other->parent_ = this->parent_;

    return this->label_;
  }


 private:
  //! Halving variant of the disjoint-set find operation
  DisjointSubset* FindRoot() {
    return const_cast<DisjointSubset*>(static_cast<const DisjointSubset*>(this)->FindRoot());
  }

  //! Halving variant of the disjoint-set find operation
  const DisjointSubset* FindRoot() const {
    if (!parent_) {
      return this;
    }
    auto current = this;

    auto parent = current->parent_;
    auto grand_parent = parent->parent_;

    while(parent != grand_parent) {
      current->parent_ = grand_parent;
      current = grand_parent;

      parent = current->parent_;
      grand_parent = parent->parent_;
    }

    return parent;
  }

  size_t size_;
  mutable DisjointSubset* parent_;
  Label label_;
};

#endif //_CRAG_DISJOINTSUBSET_H_
