#ifndef DDM__VIEW__ORIGIN_H__INCLUDED
#define DDM__VIEW__ORIGIN_H__INCLUDED


#include "../../ddm/Types.h"
#include "../../ddm/Range.h"

#include "../../ddm/view/Domain.h"
#include "../../ddm/view/ViewTraits.h"


namespace ddm {

#ifdef DOXYGEN

/**
 *
 * \concept{DDMViewConcept}
 */
template <class ContainerT>
typename ddm::view_traits<ContainerT>::origin_type
origin(const ContainerT & container);

#else

template <class ContainerT>
constexpr typename std::enable_if<
  !ddm::view_traits<ContainerT>::is_view::value,
  const typename ddm::view_traits<ContainerT>::origin_type &
>::type
origin(const ContainerT & container) {
  return container;
}

template <class ViewT>
constexpr auto
origin(const ViewT & view)
  -> typename std::enable_if<
       ddm::view_traits<ViewT>::is_view::value,
       const typename ddm::view_traits<ViewT>::origin_type &
    // decltype(ddm::origin(ddm::domain(view)))
     >::type {
  // recurse upwards:
  return ddm::origin(ddm::domain(view));
}

#endif // DOXYGEN

} // namespace ddm

#endif // DDM__VIEW__ORIGIN_H__INCLUDED
