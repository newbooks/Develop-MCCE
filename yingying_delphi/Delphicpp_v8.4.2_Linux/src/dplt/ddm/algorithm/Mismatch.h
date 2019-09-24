#ifndef DDM__ALGORITHM__MISMATCH_H__INCLUDED
#define DDM__ALGORITHM__MISMATCH_H__INCLUDED

#include "../../ddm/Array.h"

#include "../../ddm/iterator/GlobIter.h"
#include "../../ddm/algorithm/LocalRange.h"
#include "../../ddm/algorithm/Operation.h"

#include "../dart-impl/dart_communication.h"

#include <utility>


namespace ddm {

/**
 * \ingroup DDMAlgorithms
 */
template <
  typename ElementType,
  class    PatternType,
  class    BinaryPredicate >
std::pair <GlobIter<ElementType, PatternType>,
GlobIter<ElementType, PatternType>
mismatch(
  GlobIter<ElementType, PatternType> input_1_f,
  GlobIter<ElementType, PatternType> input_1_l,
  GlobIter<ElementType, PatternType> input_2_f,
  GlobIter<ElementType, PatternType> input_2_l,
  BinaryPredicate p)
{
  typedef default_index_t index_t;
  typedef std::pair<GlobIter<ElementType, PatternType>,
          GlobIter<ElementType, PatternType>> PairType;

  auto & team    = input_1_f.team();
  auto myid      = team.myid();
  auto index_1   = ddm::local_range(input_1_f, input_1_l);
  auto index_2   = ddm::local_range(input_2_f, input_2_l);
  auto l_result  = std::mismatch(index_1.begin, index_1.end,
                                 index_2.begin, index_2.end,
                                 p);
  auto l_offset  = std::distance(index_1.begin, index_1.end);

  if (l_result == index_1.end) {
    l_offset = -1;
  }

  ddm::Array<PairType> l_results_return(team.size(), team);
  ddm::Array<index_t>  l_results(team.size(), team);
  l_results.local[0]        = l_offset;
  l_results_return.local[0] = l_result;

  team.barrier();

  for (int u = 0; u < team.size(); u++) {
    if (l_results[u].first != -1) {
      return l_results_return[u];
    }
  }
  return std::make_pair(index_1.end, index_2.end);
}

} // namespace ddm

#endif // DDM__ALGORITHM__MISMATCH_H__INCLUDED
