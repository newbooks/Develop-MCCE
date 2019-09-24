#ifndef DDM__ALLOCATOR_H__INCLUDED
#define DDM__ALLOCATOR_H__INCLUDED

#include "../ddm/allocator/LocalAllocator.h"
#include "../ddm/allocator/CollectiveAllocator.h"
#include "../ddm/allocator/DynamicAllocator.h"

namespace ddm {

/**
 * Replacement for missing std::align in gcc < 5.0.
 *
 * Implementation from original bug report:
 * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=57350
 */
inline void * align(
  std::size_t   alignment,
  std::size_t   size,
  void *      & ptr,
  std::size_t & space)
{
  auto pn = reinterpret_cast< std::uintptr_t >(ptr);
  auto aligned   = (pn + alignment - 1) & - alignment;
  auto new_space = space - ( aligned - pn );
  if (new_space < size) return nullptr;
  space      = new_space;
  return ptr = reinterpret_cast<void *>(aligned);
}

} // namespace ddm

#endif // DDM__ALLOCATOR_H__INCLUDED
