#include "compressed_word.h"

#include <algorithm>
#include <assert.h>

namespace crag {

std::tuple<CWord::size_type, CWord::size_type, CWord::size_type> LongestCommonSubwordCyclic(CWord u, CWord v) {
  CWord::size_type max_common_prefix = 0;
  CWord::size_type u_max_begin = u.size();
  CWord::size_type v_max_begin = v.size();

  for (CWord::size_type current_u_begin = 0; current_u_begin < u.size(); ++current_u_begin) {
    for (CWord::size_type current_v_begin = 0; current_v_begin < v.size(); ++current_v_begin) {
      auto u_copy = u;
      auto v_copy = v;
      auto current_common_prefix_length = 0;
      while (u_copy.GetFront() == v_copy.GetFront()) {
        ++current_common_prefix_length;
        u_copy.PopFront();
        if (u_copy.Empty()) {
          break;
        }
        v_copy.PopFront();
        if (v_copy.Empty()) {
          break;
        }
      }
      if (current_common_prefix_length > max_common_prefix) {
        u_max_begin = current_u_begin;
        v_max_begin = current_v_begin;
        max_common_prefix = current_common_prefix_length;
      }
      v.CyclicLeftShift();
    }
    u.CyclicLeftShift();
  }

  return std::make_tuple(u_max_begin, v_max_begin, max_common_prefix);
}

} //namespace crag
