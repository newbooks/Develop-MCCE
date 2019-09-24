#ifndef DDM__VIEW__DOMAIN_H__INCLUDED
#define DDM__VIEW__DOMAIN_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"

#include "../../ddm/view/ViewTraits.h"


namespace ddm {

// ------------------------------------------------------------------------
// ddm::domain(View)

/**
 *
 * \concept{DDMViewConcept}
 */
template <
  class    ViewT,
  typename ViewValueT =
             typename std::remove_const<
               typename std::remove_reference<ViewT>::type
             >::type
>
constexpr auto
domain(ViewT && view)
  -> typename std::enable_if<
       ddm::view_traits<ViewValueT>::is_view::value,
    // const typename ddm::view_traits<ViewValueT>::domain_type &
       decltype(view.domain())
     >::type {
  return view.domain();
}

// ------------------------------------------------------------------------
// ddm::domain(Container)

/**
 *
 * \concept{DDMViewConcept}
 */
template <class ContainerT>
constexpr typename std::enable_if<
  !ddm::view_traits<ContainerT>::is_view::value,
  const ContainerT &
>::type
domain(const ContainerT & container) {
  return container;
}

} // namespace ddm

#endif // DDM__VIEW__DOMAIN_H__INCLUDED
