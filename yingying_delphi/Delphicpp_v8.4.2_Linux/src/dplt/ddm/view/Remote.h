#ifndef DDM__VIEW__REMOTE_H__INCLUDED
#define DDM__VIEW__REMOTE_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"

#include "../../ddm/view/ViewTraits.h"


namespace ddm {

/**
 * \concept{DDMViewConcept}
 */
template <class ViewType>
constexpr auto
remote(ddm::team_unit_t unit, const ViewType & v)
-> typename std::enable_if<
     ddm::view_traits<ViewType>::is_view::value,
     decltype(v.local())
//   const typename ViewType::local_type
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
remote(ddm::team_unit_t unit, const ContainerType & c) {
  return c.local;
}


} // namespace ddm

#endif // DDM__VIEW__REMOTE_H__INCLUDED
