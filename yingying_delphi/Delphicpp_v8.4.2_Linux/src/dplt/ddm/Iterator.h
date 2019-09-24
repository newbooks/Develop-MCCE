#ifndef DDM__ITERATOR_H__INCLUDED
#define DDM__ITERATOR_H__INCLUDED

#include "../ddm/Types.h"
#include "../ddm/Dimensional.h"
#include "../ddm/iterator/GlobIter.h"
#include "../ddm/iterator/GlobViewIter.h"

#include <iterator>

/**
 * \defgroup  DDMIteratorConcept  Multidimensional Iterator Concept
 *
 * \ingroup DDMNDimConcepts
 * \{
 * \par Description
 *
 * Definitions for multidimensional iterator expressions.
 *
 * \see DDMDimensionalConcept
 * \see DDMViewConcept
 * \see DDMRangeConcept
 *
 * \see \c ddm::view_traits
 *
 * \par Types
 *
 * \par Expressions
 *
 * Expression               | Returns | Effect | Precondition | Postcondition
 * ------------------------ | ------- | ------ | ------------ | -------------
 *
 * \par Metafunctions
 *
 * - \c ddm::iterator_traits<I>
 *
 * \par Functions
 *
 * - \c ddm::index
 *
 * \par Functions in the Range concept
 *
 * - \c ddm::distance
 * - \c ddm::size
 *
 * \}
 */

namespace ddm {

/**
 *
 * \concept{DDMIteratorConcept}
 */
template <class IndexType>
constexpr typename std::enable_if<
  std::is_integral<IndexType>::value, IndexType >::type
index(IndexType idx) {
  return idx;
}

/**
 *
 * \concept{DDMIteratorConcept}
 */
template <class Iterator>
constexpr auto index(Iterator it) -> decltype((++it).pos()) {
  return it.pos();
}

/**
 * Resolve the number of elements between two global iterators.
 *
 * The difference of global pointers is not well-defined if their range
 * spans over more than one block.
 * The corresponding invariant is:
 *   g_last == g_first + (l_last - l_first)
 * Example:
 *
 * \code
 *   unit:            0       1       0
 *   local offset:  | 0 1 2 | 0 1 2 | 3 4 5 | ...
 *   global offset: | 0 1 2   3 4 5   6 7 8   ...
 *   range:          [- - -           - -]
 * \endcode
 *
 * When iterating in local memory range [0,5[ of unit 0, the position of
 * the global iterator to return is 8 != 5
 *
 * \tparam      ElementType  Type of the elements in the range
 * \complexity  O(1)
 *
 * \ingroup     Algorithms
 *
 * \concept{DDMIteratorConcept}
 */
template<
  typename ElementType,
  class    Pattern,
  class    GlobMem,
  class    Pointer,
  class    Reference >
typename Pattern::index_type
distance(
  /// Global pointer to the initial position in the global sequence
  const GlobIter<ElementType, Pattern, GlobMem, Pointer, Reference> &
    first,
  /// Global iterator to the final position in the global sequence
  const GlobIter<ElementType, Pattern, GlobMem, Pointer, Reference> &
    last)
{
  return last - first;
}

/**
 *
 * \ingroup     Algorithms
 *
 * \concept{DDMIteratorConcept}
 */
template <class T>
constexpr std::ptrdiff_t distance(T * const first, T * const last) {
  return std::distance(first, last);
}

/**
 * Resolve the number of elements between two global pointers.
 * The difference of global pointers is not well-defined if their range
 * spans over more than one block.
 * The corresponding invariant is:
 *
 * \code
 *   g_last == g_first + (l_last - l_first)
 * \endcode
 *
 * \code
 * Example:
 *   unit:            0       1       0
 *   local offset:  | 0 1 2 | 0 1 2 | 3 4 5 | ...
 *   global offset: | 0 1 2   3 4 5   6 7 8   ...
 *   range:          [- - -           - -]
 * \endcode
 *
 * When iterating in local memory range [0,5[ of unit 0, the position of the
 * global iterator to return is 8 != 5
 *
 * \tparam      ElementType  Type of the elements in the range
 * \complexity  O(1)
 *
 * \ingroup     Algorithms
 * 
 * \concept{DDMIteratorConcept}
 */
template<typename ElementType>
ddm::default_index_t distance(
  /// Global pointer to the initial position in the global sequence
  dart_gptr_t first,
  /// Global pointer to the final position in the global sequence
  dart_gptr_t last)
{
  GlobPtr<ElementType> & gptr_first(first);
  GlobPtr<ElementType> & gptr_last(last);
  return gptr_last - gptr_first;
}

/**
 *
 * \ingroup     Algorithms
 *
 * \concept{DDMIteratorConcept}
 */
template <
  class OffsetType >
constexpr typename std::enable_if<
  std::is_integral<OffsetType>::value,
  OffsetType >::type
distance(
  OffsetType begin,
  OffsetType end) {
  return (end - begin);
}

} // namespace ddm

#endif // DDM__ITERATOR_H__INCLUDED
