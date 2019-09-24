#ifndef DDM__ALGORITHM__FOR_EACH_H__
#define DDM__ALGORITHM__FOR_EACH_H__

#include "../../ddm/iterator/GlobIter.h"
#include "../../ddm/algorithm/LocalRange.h"

#include <algorithm>


namespace ddm {

/**
 * Invoke a function on every element in a range distributed by a pattern.
 * This function has the same signature as \c std::for_each but
 * Being a collaborative operation, each unit will invoke the given
 * function on its local elements only.
 * To support compiler optimization, this const version is provided
 *
 * \tparam      ElementType   Type of the elements in the sequence
 * \tparam      UnaryFunction Function to invoke for each element
 *                            in the specified range with signature
 *                            \c (void (const ElementType &)).
 *                            Signature does not need to have \c (const &)
 *                            but must be compatible to \c std::for_each.
 *
 * \complexity  O(d) + O(nl), with \c d dimensions in the global iterators'
 *              pattern and \c nl local elements within the global range
 *
 * \ingroup     DDMAlgorithms
 */
template <
  typename ElementType,
  class    PatternType,
  class    UnaryFunction >
void for_each(
  /// Iterator to the initial position in the sequence
  const GlobIter<ElementType, PatternType> & first,
  /// Iterator to the final position in the sequence
  const GlobIter<ElementType, PatternType> & last,
  /// Function to invoke on every index in the range
  UnaryFunction                              func)
{
  /// Global iterators to local index range:
  auto index_range  = ddm::local_index_range(first, last);
  auto lbegin_index = index_range.begin;
  auto lend_index   = index_range.end;
  auto & team       = first.pattern().team();
  if (lbegin_index != lend_index) {
    // Pattern from global begin iterator:
    auto & pattern    = first.pattern();
    // Local range to native pointers:
    auto lrange_begin = (first + pattern.global(lbegin_index)).local();
    auto lrange_end   = lrange_begin + lend_index;
    std::for_each(lrange_begin, lrange_end, func);
  }
  team.barrier();
}

/**
 * Invoke a function on every element in a range distributed by a pattern.
 * Being a collaborative operation, each unit will invoke the given
 * function on its local elements only. The index passed to the function is
 * a global index.
 *
 * \tparam      ElementType            Type of the elements in the sequence
 * \tparam      UnaryFunctionWithIndex Function to invoke for each element
 *                                     in the specified range with signature
 *                                     \c void (const ElementType &, index_t)
 *                                     Signature does not need to have
 *                                     \c (const &) but must be compatible
 *                                     to \c std::for_each.
 *
 * \complexity  O(d) + O(nl), with \c d dimensions in the global iterators'
 *              pattern and \c nl local elements within the global range
 *
 * \ingroup     DDMAlgorithms
 */
template <
  typename ElementType,
  class    PatternType,
  class    UnaryFunctionWithIndex >
void for_each_with_index(
  /// Iterator to the initial position in the sequence
  const GlobIter<ElementType, PatternType> & first,
  /// Iterator to the final position in the sequence
  const GlobIter<ElementType, PatternType> & last,
  /// Function to invoke on every index in the range
  UnaryFunctionWithIndex                     func)
{
  /// Global iterators to local index range:
  auto index_range  = ddm::local_index_range(first, last);
  auto lbegin_index = index_range.begin;
  auto lend_index   = index_range.end;
  auto & team       = first.pattern().team();
  if (lbegin_index != lend_index) {
    // Pattern from global begin iterator:
    auto & pattern    = first.pattern();
    // Iterate local index range:
    for (auto lindex = lbegin_index;
         lindex != lend_index;
         ++lindex) {
      auto gindex       = pattern.global(lindex);
      auto first_offset = first.pos();
      auto element_it   = first + (gindex - first_offset);
      func(*(element_it.local()), gindex);
    }
  }
  team.barrier();
}

} // namespace ddm

#endif // DDM__ALGORITHM__FOR_EACH_H__
