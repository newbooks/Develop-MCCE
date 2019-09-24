#ifndef DDM__ALGORITHM__NONE_OF_H__
#define DDM__ALGORITHM__NONE_OF_H__

#include "../../ddm/iterator/GlobIter.h"
#include "../../ddm/algorithm/LocalRange.h"
#include "../../ddm/algorithm/Operation.h"
#include "../dart-impl/dart_communication.h"
#include "../../ddm/algorithm/Find_if.h"

namespace ddm {

/**
 * \ingroup     DDMAlgorithms
 */
template<
    typename ElementType,
    class PatternType>
GlobIter<ElementType, PatternType> any_of(
    /// Iterator to the initial position in the sequence
    GlobIter<ElementType, PatternType>   first,
    /// Iterator to the final position in the sequence
    GlobIter<ElementType, PatternType>   last,
    /// Predicate applied to the elements in range [first, last)
    UnaryPredicate p)
{
    
    return find_if(first, last, p) == last;

} // namespace ddm

#endif // DDM__ALGORITHM__NONE_OF_H__
