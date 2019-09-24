#ifndef DDM__VIEW__SUB_H__INCLUDED
#define DDM__VIEW__SUB_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"

#include "../../ddm/view/ViewMod.h"
#include "../../ddm/view/NViewMod.h"


namespace ddm {

// -------------------------------------------------------------------------
// View Modifiers (not coupled with origin memory / index space):
// -------------------------------------------------------------------------

// Sub-space slice, view dimensions maintain domain dimensions

/**
 * \concept{DDMViewConcept}
 */
template <
  dim_t SubDim   = 0,
  dim_t NViewDim,
  class OffsetFirstT,
  class OffsetFinalT >
constexpr ViewSubMod<ViewOrigin<NViewDim>, SubDim>
sub(OffsetFirstT begin,
    OffsetFinalT end) {
  return ViewSubMod<ViewOrigin<NViewDim>, SubDim>(begin, end);
}

/**
 * \concept{DDMViewConcept}
 */
template <
  dim_t SubDim   = 0,
  dim_t NViewDim,
  class IndexRangeT >
constexpr ViewSubMod<ViewOrigin<NViewDim>, SubDim>
sub(const IndexRangeT & range) {
  return sub<SubDim>(ddm::begin(range),
                     ddm::end(range));
}

// Sub-space projection, view reduces domain by one dimension

#if 0
/**
 * \concept{DDMViewConcept}
 */
template <
  dim_t SubDim = 0,
  class OffsetT >
constexpr ViewSubMod<ViewOrigin, SubDim>
sub(
    OffsetT offset) {
  return ViewSubMod<ViewOrigin, SubDim>(offset);
}
#endif

// -------------------------------------------------------------------------
// View Proxies (coupled with origin memory / index space):
// -------------------------------------------------------------------------

// Sub-space slice, view dimensions maintain domain dimensions

#if 0
/**
 * \concept{DDMViewConcept}
 */
template <
  dim_t    SubDim  = 0,
  class    DomainT,
  class    OffsetFirstT,
  class    OffsetFinalT,
  typename DomainValueT =
             typename std::remove_const<
               typename std::remove_reference<DomainT>::type
             >::type
>
constexpr auto
sub(
    OffsetFirstT    begin,
    OffsetFinalT    end,
    DomainT       & domain)
  -> typename std::enable_if<
       ddm::view_traits<
         DomainValueT
       >::rank::value == 1,
       ViewSubMod<DomainT, SubDim>
     >::type {
  return ViewSubMod<DomainT, SubDim>(
           domain,
           begin,
           end);
}
#endif

#if 1
template <
  dim_t SubDim  = 0,
  class DomainT,
  class OffsetFirstT,
  class OffsetFinalT,
  typename DomainValueT =
  //         typename std::remove_const<
               typename std::remove_reference<DomainT>::type
  //         >::type
>
constexpr auto
sub(
    OffsetFirstT    begin,
    OffsetFinalT    end,
    DomainT      && domain)
  -> typename std::enable_if<
       ddm::view_traits<
         DomainValueT
       >::rank::value == 1,
       ViewSubMod<DomainValueT, SubDim>
     >::type {
  return ViewSubMod<DomainValueT, SubDim>(
           std::forward<DomainT>(domain),
           begin,
           end);
}
#endif

// =========================================================================
// Multidimensional Views
// =========================================================================

template <
  dim_t SubDim  = 0,
  class DomainT,
  class OffsetFirstT,
  class OffsetFinalT >
constexpr auto
sub(
    OffsetFirstT    begin,
    OffsetFinalT    end,
    const DomainT & domain)
  -> typename std::enable_if<
       (ddm::view_traits<DomainT>::rank::value > 1),
       NViewSubMod<DomainT, SubDim, ddm::view_traits<DomainT>::rank::value>
     >::type {
  return NViewSubMod<
           DomainT,
           SubDim,
           ddm::view_traits<DomainT>::rank::value
         >(domain,
           begin,
           end);
}

} // namespace ddm

#endif // DDM__VIEW__SUB_H__INCLUDED
