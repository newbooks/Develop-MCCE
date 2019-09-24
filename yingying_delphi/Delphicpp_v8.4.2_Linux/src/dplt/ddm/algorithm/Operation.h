#ifndef DDM__ALGORITHM__OPERATION_H__
#define DDM__ALGORITHM__OPERATION_H__

#include "../dart-impl/dart_types.h"
#include <functional>

/**
 * \defgroup DDMReduceOperations DDM Reduce Operations
 *
 * Distributed reduce operation types.
 *
 */

namespace ddm {

/**
 * Base type of all reduce operations, primarily acts as a container of a
 * \c dart_operation_t.
 *
 * \ingroup  DDMReduceOperations
 */
template< typename ValueType, dart_operation_t OP >
class ReduceOperation {

public:
  typedef ValueType value_type;

public:

  constexpr dart_operation_t dart_operation() const {
    return _op;
  }

private:
  static constexpr const dart_operation_t _op = OP;

};

/**
 * Reduce operands to their minimum value.
 *
 * \see      dart_operation_t::DART_OP_MIN
 *
 * \ingroup  DDMReduceOperations
 */
template< typename ValueType >
struct min : public ReduceOperation<ValueType, DART_OP_MIN> {

public:

  ValueType operator()(
    const ValueType & lhs,
    const ValueType & rhs) const {
    return (lhs < rhs) ? lhs : rhs;
  }
};

/**
 * Reduce operands to their maximum value.
 *
 * \see      dart_operation_t::DART_OP_MAX
 *
 * \ingroup  DDMReduceOperations
 */
template< typename ValueType >
struct max : public ReduceOperation<ValueType, DART_OP_MAX> {

public:

  ValueType operator()(
    const ValueType & lhs,
    const ValueType & rhs) const {
    return (lhs > rhs) ? lhs : rhs;
  }
};

/**
 * Reduce operands to their sum.
 *
 * \see      dart_operation_t::DART_OP_SUM
 *
 * \ingroup  DDMReduceOperations
 */
template< typename ValueType >
struct plus : public ReduceOperation<ValueType, DART_OP_SUM> {

public:

  ValueType operator()(
    const ValueType & lhs,
    const ValueType & rhs) const {
    return lhs + rhs;
  }
};

/**
 * Reduce operands to their product.
 *
 * \see      dart_operation_t::DART_OP_PROD
 *
 * \ingroup  DDMReduceOperations
 */
template< typename ValueType >
struct multiply : public ReduceOperation<ValueType, DART_OP_PROD> {

public:

  ValueType operator()(
    const ValueType & lhs,
    const ValueType & rhs) const {
    return lhs * rhs;
  }
};

/**
 * Returns second operand. Used as replace reduce operation
 *
 * \see      dart_operation_t::DART_OP_REPLACE
 *
 * \ingroup  DDMReduceOperations
 */
template< typename ValueType >
struct second : public ReduceOperation<ValueType, DART_OP_REPLACE> {

public:

  ValueType operator()(
    const ValueType & lhs,
    const ValueType & rhs) const {
    return rhs;
  }
};

}  // namespace ddm

#endif // DDM__ALGORITHM__OPERATION_H__
