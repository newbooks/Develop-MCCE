#ifndef DDM__VIEW__GLOBAL_H__INCLUDED
#define DDM__VIEW__GLOBAL_H__INCLUDED

#include "../../ddm/view/ViewTraits.h"


namespace ddm {

/**
 * \concept{DDMViewConcept}
 */
template <class ViewType>
constexpr
typename std::enable_if<
  ddm::view_traits<ViewType>::is_view::value &&
  ddm::view_traits<ViewType>::is_local::value,
  const typename ViewType::global_type &
>::type
global(const ViewType & v) {
  return v.global();
}

/**
 * \concept{DDMViewConcept}
 */
template <class ContainerType>
constexpr
typename std::enable_if<
  !ddm::view_traits<ContainerType>::is_view::value ||
  !ddm::view_traits<ContainerType>::is_local::value,
  ContainerType &
>::type
global(ContainerType & c) {
  return c;
}

} // namespace ddm

#endif // DDM__VIEW__GLOBAL_H__INCLUDED
