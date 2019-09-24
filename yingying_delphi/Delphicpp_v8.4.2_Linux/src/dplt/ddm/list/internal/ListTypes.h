#ifndef DDM__INTERNAL__LIST__LIST_TYPES_H__INCLUDED
#define DDM__INTERNAL__LIST__LIST_TYPES_H__INCLUDED

#include "../../dart-impl/dart_types.h"

namespace ddm {
namespace internal {

template<typename ElementType>
struct ListNode
{
private:
  typedef ListNode<ElementType> self_t;

public:
  ElementType  value;
  self_t     * lprev = nullptr;
  self_t     * lnext = nullptr;
  dart_gptr_t  gprev = DART_GPTR_NULL;
  dart_gptr_t  gnext = DART_GPTR_NULL;
};

} // namespace internal
} // namespace ddm

#endif // DDM__INTERNAL__LIST__LIST_TYPES_H__INCLUDED
