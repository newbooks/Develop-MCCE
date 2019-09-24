#ifndef DDM__VIEW__INDEX_RANGE_H__INCLUDED
#define DDM__VIEW__INDEX_RANGE_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"

#include <array>


namespace ddm {

// N-dimensional index range
template <
  dim_t    NDim,
  typename IndexType >
class IndexRange
{
  typedef IndexRange<NDim, IndexType> self_t;

  // One-dimensional index range in every dimension:
  std::array< IndexRange<1, IndexType>, NDim > _ranges;

public:
  template <dim_t SDim>
  constexpr self_t sub(IndexType first, IndexType last) const {
    return self_t(*this); // _ranges[SDim].sub(first, last)
  }
};

// Specialization for one-dimensional index range
template <
  typename IndexType >
class IndexRange<1, IndexType>
{
  typedef IndexRange<1, IndexType> self_t;

  IndexType _first;
  IndexType _last;

public:
  constexpr IndexRange(IndexType first, IndexType last)
    : _first(first)
    , _last(last)
  { }

  template <typename RangeT>
  explicit constexpr IndexRange(RangeT range)
    : IndexRange(ddm::begin(range),
                 ddm::end(range))
  { }
};

} // namespace ddm

#endif // DDM__VIEW__INDEX_RANGE_H__INCLUDED
