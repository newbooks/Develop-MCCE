#ifndef DDM__ALGORITHM__TRANSFORM_H__
#define DDM__ALGORITHM__TRANSFORM_H__

#include "../../ddm/GlobRef.h"
#include "../../ddm/GlobAsyncRef.h"

#include "../../ddm/algorithm/LocalRange.h"
#include "../../ddm/algorithm/Operation.h"
#include "../../ddm/algorithm/Accumulate.h"

#include "../../ddm/iterator/GlobIter.h"

#include "../../ddm/internal/Config.h"
#include "../../ddm/util/Trace.h"

#include "../dart-impl/dart_communication.h"

#include <iterator>

#ifdef DDM_ENABLE_OPENMP
#include <omp.h>
#endif

namespace ddm {

namespace internal {

/**
 * Wrapper of the blocking DART accumulate operation.
 */
template< typename ValueType >
inline dart_ret_t transform_blocking_impl(
  dart_gptr_t        dest,
  ValueType        * values,
  size_t             nvalues,
  dart_operation_t   op)
{
  static_assert(ddm::dart_datatype<ValueType>::value != DART_TYPE_UNDEFINED,
      "Cannot accumulate unknown type!");

  dart_ret_t result = dart_accumulate(
                        dest,
                        reinterpret_cast<void *>(values),
                        nvalues,
                        ddm::dart_datatype<ValueType>::value,
                        op);
  dart_flush(dest);
  return result;
}

/**
 * Wrapper of the non-blocking DART accumulate operation.
 */
template< typename ValueType >
dart_ret_t transform_impl(
  dart_gptr_t        dest,
  ValueType        * values,
  size_t             nvalues,
  dart_operation_t   op)
{
  static_assert(ddm::dart_datatype<ValueType>::value != DART_TYPE_UNDEFINED,
      "Cannot accumulate unknown type!");

  dart_ret_t result = dart_accumulate(
                        dest,
                        reinterpret_cast<void *>(values),
                        nvalues,
                        ddm::dart_datatype<ValueType>::value,
                        op);
  dart_flush_local(dest);
  return result;
}

} // namespace internal

/**
 * Apply a given function to elements in a range and store the result in
 * another range, beginning at \c out_first.
 * Corresponding to \c MPI_Accumulate, the unary operation is executed
 * atomically on single elements.
 *
 * Precondition: All elements in the input range are contained in a single
 * block so that
 *
 * <tt>
 *   g_out_last == g_out_first + (l_in_last - l_in_first)
 * </tt>
 *
 * Semantics:
 *
 * <tt>
 *   unary_op(in_first[0]), unary_op(in_first[1]), ..., unary_op(in_first[n])
 * </tt>
 *
 * \returns  Output iterator to the element past the last element transformed.
 *
 * \ingroup  DDMAlgorithms
 */
template<
  typename ValueType,
  class InputIt,
  class OutputIt,
  class UnaryOperation >
OutputIt transform(
  InputIt        in_first,
  InputIt        in_last,
  OutputIt       out_first,
  UnaryOperation unary_op);

/**
 * Apply a given function to pairs of elements from two ranges and store the
 * result in another range, beginning at \c out_first.
 *
 * Corresponding to \c MPI_Accumulate, the binary operation is executed
 * atomically on single elements.
 *
 * Precondition: All elements in the input range are contained in a single
 * block so that
 *
 *   g_out_last == g_out_first + (l_in_last - l_in_first)
 *
 * Semantics:
 *
 *   binary_op(in_a[0], in_b[0]),
 *   binary_op(in_a[1], in_b[1]),
 *   ...,
 *   binary_op(in_a[n], in_b[n])
 *
 * Example:
 * \code
 *   gptr_diff_t num_transformed_elements =
 *                 ddm::distance(
 *                   ddm::transform(in.begin(), in.end(),  // A
 *                                   out.begin(),           // B
 *                                   out.begin(),           // C = op(A, B)
 *                                   ddm::plus<int>()),    // op
 *                   out.end());
 *
 * \endcode
 *
 * \returns  Output iterator to the element past the last element transformed.
 * \see      ddm::accumulate
 * \see      DDMReduceOperations
 *
 * \tparam   InputIt         Iterator on first (local) input range
 * \tparam   GlobInputIt     Iterator on second (global) input range
 * \tparam   GlobOutputIt    Iterator on global result range
 * \tparam   BinaryOperation Reduce operation type
 *
 * \ingroup  DDMAlgorithms
 */
template<
  typename ValueType,
  class InputIt,
  class GlobInputIt,
  class GlobOutputIt,
  class BinaryOperation >
GlobOutputIt transform(
  /// Iterator on begin of first local range
  InputIt         in_a_first,
  /// Iterator after last element of local range
  InputIt         in_a_last,
  /// Iterator on begin of second local range
  GlobInputIt     in_b_first,
  /// Iterator on first element of global output range 
  GlobOutputIt    out_first,
  /// Reduce operation
  BinaryOperation binary_op);

/**
 * Transform operation on ranges with identical distribution and start
 * offset.
 * In this case, no communication is needed as all output values can be
 * obtained from input values in local memory:
 *
 * \note
 * This function does not execute the transformation as atomic operation
 * on elements. Use \c ddm::transform if concurrent access to elements is
 * possible.
 *
 * <pre>
 *   input a: [ u0 | u1 | u2 | ... ]
 *              op   op   op   ...
 *   input b: [ u0 | u1 | u2 | ... ]
 *              =    =    =    ...
 *   output:  [ u0 | u1 | u2 | ... ]
 * </pre>
 */
template<
  typename ValueType,
  class InputAIt,
  class InputBIt,
  class OutputIt,
  class BinaryOperation >
OutputIt transform_local(
  InputAIt        in_a_first,
  InputAIt        in_a_last,
  InputBIt        in_b_first,
  OutputIt        out_first,
  BinaryOperation binary_op)
{
  DDM_LOG_DEBUG("ddm::transform_local()");
  DDM_ASSERT_MSG(in_a_first.pattern() == in_b_first.pattern(),
                  "ddm::transform_local: "
                  "distributions of input ranges differ");
  DDM_ASSERT_MSG(in_a_first.pattern() == out_first.pattern(),
                  "ddm::transform_local: "
                  "distributions of input- and output ranges differ");
  // Local subrange of input range a:
  auto local_range_a   = ddm::local_range(in_a_first, in_a_last);
  ValueType * lbegin_a = local_range_a.begin;
  ValueType * lend_a   = local_range_a.end;
  if (lbegin_a == lend_a) {
    // Local input range is empty, return initial output iterator to indicate
    // that no values have been transformed:
    DDM_LOG_DEBUG("ddm::transform_local", "local range empty");
    return out_first;
  }
  // Global offset of first local element:
  auto g_offset_first    = in_a_first.pattern().global(0);
  // Number of elements in global ranges:
  auto num_gvalues       = ddm::distance(in_a_first, in_b_first);
  DDM_LOG_TRACE_VAR("ddm::transform_local", num_gvalues);
  // Number of local elements:
  DDM_LOG_TRACE("ddm::transform_local", "local elements:", lend_a-lbegin_a);
  // Local subrange of input range b:
  ValueType * lbegin_b   = (in_b_first + g_offset_first).local();
  // Local pointer of initial output element:
  ValueType * lbegin_out = (out_first  + g_offset_first).local();
  // Generate output values:
#ifdef DDM_ENABLE_OPENMP
  ddm::util::UnitLocality uloc;
  auto n_threads = uloc.num_domain_threads();
  DDM_LOG_DEBUG("ddm::transform_local", "thread capacity:",  n_threads);
  if (n_threads > 1) {
    auto l_size = lend_a - lbegin_a;
    // TODO: Vectorize.
    // Documentation of Intel MIC intrinsics, see:
    // https://software.intel.com/de-de/node/523533
    // https://software.intel.com/de-de/node/523387
    #pragma omp parallel for num_threads(n_threads) schedule(static)
    for (int i = 0; i < l_size; i++) {
      lbegin_out[i] = binary_op(lbegin_a[i], lbegin_b[i]);
    }
    return out_first + num_gvalues;
  }
#endif
  // No OpenMP or insufficient number of threads for parallelization:
  for (; lbegin_a != lend_a; ++lbegin_a, ++lbegin_b, ++lbegin_out) {
    *lbegin_out = binary_op(*lbegin_a, *lbegin_b);
  }
  // Return out_end iterator past final transformed element;
  return out_first + num_gvalues;
}

/**
 * Local lhs input ranges on global output range.
 *
 */
template<
  typename ValueType,
  class InputIt,
  class GlobInputIt,
  class GlobOutputIt,
  class BinaryOperation >
GlobOutputIt transform(
  InputIt         in_a_first,
  InputIt         in_a_last,
  GlobInputIt     in_b_first,
  GlobOutputIt    out_first,
  BinaryOperation binary_op)
{
  DDM_LOG_DEBUG("ddm::transform(af, al, bf, outf, binop)");
  // Outut range different from rhs input range is not supported yet
  auto in_first = in_a_first;
  auto in_last  = in_a_last;
  std::vector<ValueType> in_range;
  if (in_b_first == out_first) {
    // Output range is rhs input range: C += A
    // Input is (in_a_first, in_a_last).
  } else {
    // Output range different from rhs input range: C = A+B
    // Input is (in_a_first, in_a_last) + (in_b_first, in_b_last):
    std::transform(
      in_a_first, in_a_last,
      in_b_first,
      std::back_inserter(in_range),
      binary_op);
    in_first = in_range.data();
    in_last  = in_first + in_range.size();
  }

  ddm::util::Trace trace("transform");

  // Resolve local range from global range:
  // Number of elements in local range:
  size_t num_local_elements     = std::distance(in_first, in_last);
  // Global iterator to dart_gptr_t:
  dart_gptr_t dest_gptr         = out_first.dart_gptr();
  // Send accumulate message:
  trace.enter_state("transform_blocking");
  ddm::internal::transform_blocking_impl(
      dest_gptr,
      in_first,
      num_local_elements,
      binary_op.dart_operation());
  trace.exit_state("transform_blocking");
  // The position past the last element transformed in global element space
  // cannot be resolved from the size of the local range if the local range
  // spans over more than one block. Otherwise, the difference of two global
  // iterators is not well-defined. The corresponding invariant is:
  //   g_out_last == g_out_first + (l_in_last - l_in_first)
  // Example:
  //   unit:            0       1       0
  //   local offset:  | 0 1 2 | 0 1 2 | 3 4 5 | ...
  //   global offset: | 0 1 2   3 4 5   6 7 8   ...
  //   range:          [- - -           - -]
  // When iterating in local memory range [0,5[ of unit 0, the position of the
  // global iterator to return is 8 != 5
  // For ranges over block borders, we would have to resolve the global
  // position past the last element transformed from the iterator's pattern
  // (see ddm::PatternIterator).
  return out_first + num_local_elements;
}

/**
 * Specialization of \c ddm::transform for global lhs input range.
 */
template<
  typename ValueType,
  class PatternType,
  class GlobInputIt,
  class GlobOutputIt,
  class BinaryOperation >
GlobOutputIt transform(
  GlobIter<ValueType, PatternType> in_a_first,
  GlobIter<ValueType, PatternType> in_a_last,
  GlobInputIt                      in_b_first,
  GlobOutputIt                     out_first,
  BinaryOperation                  binary_op = ddm::plus<ValueType>())
{
  DDM_LOG_DEBUG("ddm::transform(gaf, gal, gbf, goutf, binop)");
  auto in_first = in_a_first;
  auto in_last  = in_a_last;

  if (in_b_first == out_first) {
    // Output range is rhs input range: C += A
    // Input is (in_a_first, in_a_last).
  } else {
    DDM_THROW(
      ddm::exception::NotImplemented,
      "ddm::transform is only implemented for out = op(in,out)");
    // Output range different from rhs input range: C = A+B
    // Input is (in_a_first, in_a_last) + (in_b_first, in_b_last):

    // TODO:
    // in_range.allocate(...);
    // in_first = in_range.begin();
    // in_last  = in_range.end();
  }

  ddm::util::Trace trace("transform");

  // Pattern of input ranges a and b, and output range:
  auto pattern_in_a = in_a_first.pattern();
  auto pattern_in_b = in_b_first.pattern();
  auto pattern_out  = out_first.pattern();

#if __NON_ATOMIC__
  // Fast path: check if transform operation is local-only:
  if (pattern_in_a == pattern_in_b &&
      pattern_in_a == pattern_out) {
    // Identical pattern in all ranges
    if (in_a_first.pos() == in_b_first.pos() &&
        in_a_first.pos() == out_first.pos()) {
      trace.enter_state("local");
      // All units operate on local ranges that have identical distribution:
      auto out_last = ddm::transform_local<ValueType>(
                        in_a_first,
                        in_a_last,
                        in_b_first,
                        out_first,
                        binary_op);
      trace.exit_state("local");
      return out_last;
    }
  }
#endif
  // Resolve teams from global iterators:
  ddm::Team & team_in_a        = pattern_in_a.team();
  DDM_ASSERT_MSG(
    team_in_a == pattern_in_b.team(),
    "ddm::transform: Different teams in input ranges");
  DDM_ASSERT_MSG(
    team_in_a == pattern_out.team(),
    "ddm::transform: Different teams in input- and output ranges");
  // Resolve local range from global range:
  auto l_index_range_in_a       = local_index_range(in_a_first, in_a_last);
  DDM_LOG_TRACE_VAR("ddm::transform", l_index_range_in_a.begin);
  DDM_LOG_TRACE_VAR("ddm::transform", l_index_range_in_a.end);
  // Local range to global offset:
  auto global_offset            = pattern_in_a.global(
                                    l_index_range_in_a.begin);
  DDM_LOG_TRACE_VAR("ddm::transform", global_offset);
  // Number of elements in local range:
  size_t num_local_elements     = l_index_range_in_a.end -
                                  l_index_range_in_a.begin;
  DDM_LOG_TRACE_VAR("ddm::transform", num_local_elements);
  // Global iterator to dart_gptr_t:
  dart_gptr_t dest_gptr         = (out_first + global_offset).dart_gptr();
  // Native pointer to local sub-range:
  ValueType * l_values          = (in_a_first + global_offset).local();
  // Send accumulate message:
  trace.enter_state("transform_blocking");
  ddm::internal::transform_blocking_impl(
      dest_gptr,
      l_values,
      num_local_elements,
      binary_op.dart_operation());
  trace.exit_state("transform_blocking");

  return out_first + global_offset + num_local_elements;
}

/**
 * Specialization of \c ddm::transform as non-blocking operation.
 *
 * \tparam   InputIt         Iterator on first (typically local) input range
 * \tparam   GlobInputIt     Iterator on second (typically global) input range
 * \tparam   GlobOutputIt    Iterator on global result range
 * \tparam   BinaryOperation Reduce operation type
 */
template<
  typename ValueType,
  class InputIt,
  class GlobInputIt,
  class BinaryOperation >
GlobAsyncRef<ValueType> transform(
  InputIt                 in_a_first,
  InputIt                 in_a_last,
  GlobInputIt             in_b_first,
  GlobAsyncRef<ValueType> out_first,
  BinaryOperation         binary_op  = ddm::plus<ValueType>()) {
  DDM_THROW(
    ddm::exception::NotImplemented,
    "Async variant of ddm::transform is not implemented");
}

} // namespace ddm

#endif // DDM__ALGORITHM__TRANSFORM_H__
