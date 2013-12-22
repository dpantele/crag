/*
 * \file slp_vertex-inl.h
 * 
 * Definitions of inlined functions fomr slp_vertex.h
 *
 */

namespace crag { namespace slp {

Vertex Vertex::left_child() const noexcept {
  if (nonterminal_data_) {
    assert(storage_.is_valid());
    if (terminal_data_ < 0) {
      return Vertex({storage_, nonterminal_data_->right_child_, -nonterminal_data_->right_terminal_});
    } else {
      return Vertex({storage_, nonterminal_data_->left_child_, nonterminal_data_->left_terminal_});
    }
  } else {
    return Vertex();
  }
}

Vertex Vertex::right_child() const noexcept {
  if (nonterminal_data_) {
    assert(storage_.is_valid());
    if (terminal_data_ < 0) {
      return Vertex({storage_, nonterminal_data_->left_child_, -nonterminal_data_->left_terminal_});
    } else {
      return Vertex({storage_, nonterminal_data_->right_child_, nonterminal_data_->right_terminal_});
    }
  } else {
    return Vertex();
  }
}

inline Vertex Vertex::negate() const noexcept {
  return Vertex(storage_, nonterminal_data_, -terminal_data_);
}

const LongInteger& Vertex::split_point() const noexcept {
  if (nonterminal_data_) {
    if (terminal_data_ < 0) {
      return left_child().length();
    } else {
      return right_child().length();
    }
  } else {
    return LongZero();
  }
}

size_t Vertex::vertex_hash() const noexcept {
  if (nonterminal_data_) {
    return terminal_data_ > 0 ? nonterminal_ptr_hash(nonterminal_data_) : ~nonterminal_ptr_hash(nonterminal_data_);
  } else {
    return terminal_symbol_hash(terminal_data_);
  }
}

bool Vertex::operator<(const Vertex& other) const noexcept {
  if (this->nonterminal_data_ == other.nonterminal_data_) {
    return terminal_data_ < other.terminal_data_;
  } else {
    return nonterminal_data_ < other.nonterminal_data_;
  }
}

Vertex::operator bool() const noexcept {
  return terminal_data_;
}

bool Vertex::is_terminal() const noexcept {
  return terminal_data_ && !nonterminal_data_;
}


inline void Vertex::debug_print(::std::ostream* out) const {
  if (!terminal_data_) {
    (*out) << "Vertex()";
  } else if (!nonterminal_data_) {
    (*out) << "TerminalVertex(" << terminal_data_ << ')';
  } else {
    (*out) << "NonterminalVertex(l=" << length()
           << ", h=" << height();
  }
}

#ifndef DNDEGUG
namespace internal {

constexpr VertexStorageDescriptor::VertexStorageDescriptor()
  : weak_ptr_()
  , bare_ptr_()
{ }

VertexStorageDescriptor::VertexStorageDescriptor(VertexStorageInterface* storage)
  : weak_ptr_(storage->shared_from_this())
  , bare_ptr_(storage)
{ }

VertexStorageInterface& VertexStorageDescriptor::operator*() const noexcept {
  assert(is_valid());
  return *bare_ptr_;
}

VertexStorageInterface* VertexStorageDescriptor::operator->() const noexcept {
  assert(is_valid());
  return bare_ptr_;
}

bool VertexStorageDescriptor::operator==(const VertexStorageDescriptor& other) const noexcept {
  return bare_ptr_ == other.bare_ptr_;
}

}  // namespace internal
#endif

}} //namespace crag::slp
