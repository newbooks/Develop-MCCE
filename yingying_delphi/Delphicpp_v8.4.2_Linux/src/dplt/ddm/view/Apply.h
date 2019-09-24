#ifndef DDM__VIEW__APPLY_H__INCLUDED
#define DDM__VIEW__APPLY_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"

#include "../../ddm/view/ViewMod.h"


namespace ddm {

/**
 * Inverse operation to \c ddm::domain.
 *
 * \concept{DDMViewConcept}
 */
template <class ViewTypeA, class ViewTypeB>
constexpr auto apply(
  ViewTypeA & view_a,
  ViewTypeB & view_b) -> decltype(view_a.apply(view_b)) {
  return view_a.apply(view_b);
}

/**
 * \concept{DDMViewConcept}
 */
template <class ViewType>
constexpr auto apply(
  const ViewType & view) -> decltype(view.apply()) {
  return view.apply();
}

} // namespace ddm

#endif // DDM__VIEW__APPLY_H__INCLUDED
