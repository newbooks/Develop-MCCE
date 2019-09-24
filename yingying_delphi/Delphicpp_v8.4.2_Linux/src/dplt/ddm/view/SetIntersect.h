#ifndef DDM__VIEW__SET_INTERSECT_H__INCLUDED
#define DDM__VIEW__SET_INTERSECT_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"

#include "../../ddm/view/Sub.h"


namespace ddm {

/**
 * \concept{DDMViewConcept}
 */
template <
  class ViewTypeA,
  class ViewTypeB >
constexpr auto
intersect(
  const ViewTypeA & va,
  const ViewTypeB & vb)
  -> decltype(ddm::sub(0, 0, va))
{
  return ddm::sub(
           ddm::index(va).pre()[
             std::max(
               *ddm::begin(ddm::index(va)),
               *ddm::begin(ddm::index(vb))
             )
           ],
           ddm::index(va).pre()[
             std::min(
               *ddm::end(ddm::index(va)),
               *ddm::end(ddm::index(vb))
             )
           ],
           va
         );
}

} // namespace ddm

#endif // DDM__VIEW__SET_INTERSECT_H__INCLUDED
