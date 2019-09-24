#ifndef DDM__VIEW__BLOCK_H__INCLUDED
#define DDM__VIEW__BLOCK_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"

#include "../../ddm/view/ViewTraits.h"
#include "../../ddm/view/ViewMod.h"


namespace ddm {

#if 0
/**
 *
 * \concept{DDMViewConcept}
 */
template <
  class ViewT,
  class BlockIndexT >
constexpr typename std::enable_if<
  ddm::view_traits<ViewT>::is_view::value,
  ddm::ViewBlockMod<ViewT>
>::type
block(BlockIndexT block_index, const ViewT & view) {
  return ViewBlockMod<ViewT>(view, block_index);
}
#endif

}

#endif // DDM__VIEW__BLOCK_H__INCLUDED
