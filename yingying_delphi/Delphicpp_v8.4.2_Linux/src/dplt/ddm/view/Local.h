#ifndef DDM__VIEW__LOCAL_H__INCLUDED
#define DDM__VIEW__LOCAL_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"

#include "../../ddm/view/ViewTraits.h"


namespace ddm {

/**
 * \concept{DDMViewConcept}
 */
template <class ViewType>
constexpr auto
local(ViewType & v)
-> typename std::enable_if<
     std::is_pointer< typename ViewType::iterator >::value,
     ViewType &
   >::type {
  return v;
}

/**
 * \concept{DDMViewConcept}
 */
template <class ViewType>
constexpr auto
local(const ViewType & v)
-> typename std::enable_if<
     ddm::view_traits<ViewType>::is_view::value,
     decltype(v.local())
   >::type {
  return v.local();
}

/**
 * \concept{DDMViewConcept}
 */
template <class ContainerType>
constexpr
typename std::enable_if<
  !ddm::view_traits<ContainerType>::is_view::value,
  const typename ContainerType::local_type &
>::type
local(const ContainerType & c) {
  return c.local;
}

/**
 * Convert global iterator referencing an element the active unit's
 * memory to a corresponding native pointer referencing the element.
 *
 * Precondition:  \c g_it  is local
 *
 */
template <class GlobalIterator>
constexpr auto local(
  /// Global iterator referencing element in local memory
  const GlobalIterator & g_it)
->  decltype((g_it - g_it.pos()).local()) {
  return g_it.local();
}

} // namespace ddm

#endif // DDM__VIEW__LOCAL_H__INCLUDED
