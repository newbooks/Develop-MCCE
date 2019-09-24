#ifndef DDM__VIEW__CHUNKED_H__INCLUDED
#define DDM__VIEW__CHUNKED_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"

#include "../../ddm/view/Domain.h"
#include "../../ddm/view/Local.h"
#include "../../ddm/view/Origin.h"
#include "../../ddm/view/Domain.h"
#include "../../ddm/view/ViewTraits.h"
#include "../../ddm/view/SetIntersect.h"


namespace ddm {

// ------------------------------------------------------------------------
// Forward-declarations
// ------------------------------------------------------------------------

template <
  class DomainType >
class ViewBlockMod;

#if 0
template <
  class ContainerType,
  class OffsetT >
constexpr auto
block(
  OffsetT               block_idx,
  const ContainerType & container)
-> typename std::enable_if<
     !ddm::view_traits<ContainerType>::is_view::value,
     decltype(container.block(0))
   >::type {
  return container.block(block_idx);
}
#endif

/**
 * Blocks view from global view
 *
 */
template <
  class ViewType,
  class OffsetT >
constexpr auto
block(
  OffsetT    block_idx,
  const ViewType & view)
-> typename std::enable_if<
     (//  ddm::view_traits<ViewType>::is_view::value &&
         !ddm::view_traits<ViewType>::is_local::value   ),
     ViewBlockMod<ViewType>
   >::type {
  return ViewBlockMod<ViewType>(view, block_idx);
}

/**
 * Blocks view from local view
 *
 */
template <
  class ViewType,
  class OffsetT >
constexpr auto
block(
  OffsetT          block_idx,
  const ViewType & view)
-> typename std::enable_if<
     (// ddm::view_traits<ViewType>::is_view::value &&
         ddm::view_traits<ViewType>::is_local::value   ),
     decltype(ddm::block(block_idx, ddm::local(ddm::origin(view))))
   >::type {
  return ddm::local(ddm::origin(view)).block(block_idx);
}

} // namespace ddm

#endif // DDM__VIEW__CHUNKED_H__INCLUDED
