#ifndef DDM__INTERNAL__ALLOCATOR__GLOB_DYNAMIC_MEM_TYPES_H__INCLUDED
#define DDM__INTERNAL__ALLOCATOR__GLOB_DYNAMIC_MEM_TYPES_H__INCLUDED

namespace ddm {
namespace internal {

template<
  typename SizeType,
  typename ElementType >
struct glob_dynamic_mem_bucket_type
{
  SizeType      size;
  ElementType * lptr;
  dart_gptr_t   gptr;
  bool          attached;
};

} // namespace internal
} // namespace ddm

#endif // DDM__INTERNAL__ALLOCATOR__GLOB_DYNAMIC_MEM_TYPES_H__INCLUDED
