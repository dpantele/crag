/*
 * slp_reduce.h
 *
 *  Created on: Mar 12, 2013
 *      Author: dpantele
 */

#pragma once
#ifndef CRAG_FREEGROUP_SLP_REDUCE_H_
#define CRAG_FREEGROUP_SLP_REDUCE_H_

#include <gmpxx.h>
typedef mpz_class LongInteger;

#include "slp_vertex.h"
#include "slp_mapper.h"
#include "slp_inspector.h"
#include "slp_common_prefix.h"

namespace std {
//(x)->_mp_d

template<>
struct hash<LongInteger> {
  private:
    constexpr static hash<mp_limb_t*> limb_hasher_ = hash<mp_limb_t*>();
  public:
    size_t operator()(const LongInteger& obj) const {
      return limb_hasher_(obj.get_mpz_t()->_mp_d);
    }
};

namespace tuple_hash_detail {

template <std::size_t I, typename Tuple>
class TupleHasher {
  private:
    constexpr static hash<typename std::tuple_element<I - 1, Tuple>::type> element_hasher_ = hash<typename std::tuple_element<I - 1, Tuple>::type>();
    constexpr static TupleHasher<I - 1, Tuple> prefix_hasher_ = TupleHasher<I - 1, Tuple>();
  public:
    size_t operator()(const Tuple& obj) const {
      size_t prefix_hash_value = prefix_hasher_(obj);
      //Taken from boost/functional/hash
      return element_hasher_(std::get<I - 1>(obj)) + 0x9e3779b9 + (prefix_hash_value << 6) + (prefix_hash_value >> 2);
    }
};

template <typename Tuple>
class TupleHasher<0, Tuple> {
  public:
    constexpr size_t operator()(const Tuple& obj) const {
      //Taken from boost/functional/hash
      return 0;
    }
};

template <std::size_t I, typename Tuple>
constexpr hash<typename std::tuple_element<I - 1, Tuple>::type> TupleHasher<I, Tuple>::element_hasher_;
template <std::size_t I, typename Tuple>
constexpr TupleHasher<I - 1, Tuple> TupleHasher<I, Tuple>::prefix_hasher_;
//template <typename Tuple>
//hash<typename std::tuple_element<0, Tuple>::type> TupleHasher<0, Tuple>::element_hasher_;

//template <class T>
//inline void hash_combine(std::size_t& seed, const T& v)
//{
//    hash<T> hasher;
//    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
//}
//
//template <std::size_t I, typename T>
//inline typename std::enable_if<(I == std::tuple_size<T>::value),
//        void>::type
//hash_combine_tuple(std::size_t&, const T&)
//{ }
//
//template <std::size_t I, typename T>
//inline typename std::enable_if<(I < std::tuple_size<T>::value),
//        void>::type
//hash_combine_tuple(std::size_t& seed, const T& v)
//{
//    hash_combine(seed, std::get<I>(v));
//    hash_combine_tuple<I + 1>(seed, v);
//}
//
//template <typename T>
//inline std::size_t hash_tuple(T const& v)
//{
//    std::size_t seed = 0;
//    hash_combine_tuple<0>(seed, v);
//    return seed;
//}

}//namespace tuple_hash_detail

template<typename... T>
struct hash<tuple<T...>> {
  private:
    constexpr static tuple_hash_detail::TupleHasher<std::tuple_size<tuple<T...>>::value, tuple<T...>> tuple_haser_ = tuple_hash_detail::TupleHasher<std::tuple_size<tuple<T...>>::value, tuple<T...>>();
  public:
    size_t operator()(const tuple<T...>& obj) const {
      return tuple_haser_(obj);
    }
};

template<typename... T>
constexpr tuple_hash_detail::TupleHasher<std::tuple_size<tuple<T...>>::value, tuple<T...>> hash<tuple<T...>>::tuple_haser_;

}//namespace std

namespace crag {
namespace slp {
Vertex get_sub_slp(const Vertex& root, const LongInteger& begin, const LongInteger& end);

inline LongInteger get_cancellation_length(const Vertex& vertex, MatchingTable* matching_table) {
  return longest_common_prefix(vertex.left_child().negate(), vertex.right_child(), matching_table);
}

inline LongInteger get_cancellation_length(const Vertex& vertex) {
  MatchingTable temp;
  return get_cancellation_length(vertex, &temp);
}

struct ReducedStorage : public internal::VertexStorageEntry {
    Vertex reduced_;
    ReducedStorage(Vertex&& vertex)
      : reduced_(std::move(vertex))
    { }
};

template <typename GetCancellationLengthFunctor>
class ReducedVertex : public internal::VertexData {
  public:
    typedef ReducedStorage StorageType;

    Vertex reduced_;

    static std::unique_ptr<StorageType> create_storage_entry(const Vertex& left_child,
                                                             const Vertex& right_child,
                                                             GetCancellationLengthFunctor get_cancellation_length,
                                                             const Vertex& self) {
      auto left_reduced = std::move(left_child.get_data<ReducedVertex>(get_cancellation_length, left_child).reduced_);
      auto right_reduced = std::move(right_child.get_data<ReducedVertex>(get_cancellation_length, right_child).reduced_);
      if (!left_reduced) {
        return std::unique_ptr<StorageType>(new StorageType(std::move(right_reduced)));
      }
      NonterminalVertex result(left_reduced, right_reduced);
      LongInteger cancellation_length = get_cancellation_length(result);
      if (cancellation_length == 0) {
        if (left_reduced == left_child && right_reduced == right_child) {
          return std::unique_ptr<StorageType>(new StorageType(self.internal_abs()));
        } else if (!left_reduced) {
          return std::unique_ptr<StorageType>(new StorageType(std::move(right_reduced)));
        } else if (!right_reduced) {
          return std::unique_ptr<StorageType>(new StorageType(std::move(left_reduced)));
        } else {
          assert((result.height() > 1 && result.length() > 1) || (result.length() == 1 && result.height() == 1) || (result.length() == 0 && result.height() == 0));
          return std::unique_ptr<StorageType>(new StorageType(std::move(result)));
        }
      } else {
        Vertex remaining_prefix = get_sub_slp(left_reduced, 0, left_reduced.length() - cancellation_length);
        Vertex remaining_suffix = get_sub_slp(right_reduced, cancellation_length, right_reduced.length());

        if (!remaining_prefix && !remaining_suffix) {
          return std::unique_ptr<StorageType>(new StorageType(Vertex()));
        } else if (!remaining_prefix) {
          assert(remaining_suffix.height() >= 1);
          assert(remaining_suffix.height() != 1 || remaining_suffix.length() == 1);
          return std::unique_ptr<StorageType>(new StorageType(std::move(remaining_suffix)));
        } else if (!remaining_suffix) {
          assert(remaining_prefix.height() >= 1);
          assert(remaining_prefix.height() != 1 || remaining_prefix.length() == 1);
          return std::unique_ptr<StorageType>(new StorageType(std::move(remaining_prefix)));
        } else {
          auto result = NonterminalVertex(remaining_prefix, remaining_suffix);
          assert(result.length() > 1);
          assert(result.height() > 1);
          return std::unique_ptr<StorageType>(new StorageType(std::move(result)));
        }
      }
    }

    ReducedVertex(StorageType* storage, bool negate, GetCancellationLengthFunctor, const Vertex&)
      : reduced_(negate ? storage->reduced_.negate() : storage->reduced_)
    {
    }

    ReducedVertex(const Vertex& vertex, GetCancellationLengthFunctor, const Vertex&)
      : reduced_(vertex)
    { }
};

template <typename GetCancellationLengthFunctor>
inline Vertex base_reduce(
    const Vertex& vertex,
    GetCancellationLengthFunctor get_cancellation_length,
    std::unordered_map<Vertex, Vertex>* reduced_vertices)
{
  return vertex.get_data<ReducedVertex<GetCancellationLengthFunctor>>(get_cancellation_length, vertex).reduced_;
}

inline Vertex reduce(
    const Vertex& vertex,
    MatchingTable* matching_table,
    std::unordered_map<Vertex, Vertex>* reduced_vertices)
{
  return base_reduce(vertex,
                     [matching_table](const Vertex& vertex) { return get_cancellation_length(vertex, matching_table); },
                     reduced_vertices
  );
}


inline Vertex reduce(const Vertex& vertex) {
  MatchingTable matching_table;
  std::unordered_map<Vertex, Vertex> reduced_vertices;
  return reduce(vertex, &matching_table, &reduced_vertices);
}

}
}
#endif /* SLP_REDUCE_H_ */
