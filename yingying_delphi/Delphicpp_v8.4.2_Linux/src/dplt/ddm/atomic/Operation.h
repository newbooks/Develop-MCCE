#ifndef DDM__ATOMIC_OPERATION_H_
#define DDM__ATOMIC_OPERATION_H_

#include "../../ddm/atomic/GlobAtomicRef.h"

namespace ddm {

// forward decls
template<typename T>
class Atomic;

/**
 * Routines to perform atomic operations on atomics residing in the
 * global address space
 *
 * \code
 *  int N = ddm::size();
 *  ddm::Array<ddm::Atomic<int>> array(N);
 *  ddm::fill(array.begin(), array.end(), 0);
 *  // each unit adds 1 to each array position
 *  for(auto & el : array){
 *    ddm::atomic::add(el, 1);
 *  }
 *  // postcondition:
 *  // array = {N,N,N,N,...}
 * \endcode
 */
namespace atomic {

/**
 * Get the value of the shared atomic.
 */
template<typename T>
T load(const ddm::GlobRef<ddm::Atomic<T>> & ref){
  return ref.load();
}

/**
 * Set the value of the atomic reference 
 */
template<typename T>
void store(const ddm::GlobRef<ddm::Atomic<T>> & ref,
           const T & value)
{
  ref.store(value);
}

/**
 * Atomically sets the value of the atomic reference and returns
 * the old value
 */
template<typename T>
T exchange(const ddm::GlobRef<ddm::Atomic<T>> & ref,
           const T & value)
{
  return ref.exchange(value);
}

/**
 * Atomically compares the value with the value of expected and if those are
 * bitwise-equal, replaces the former with desired.
 * 
 * \return  True if value is exchanged
 */
template<typename T>
bool compare_exchange(
  const ddm::GlobRef<ddm::Atomic<T>> & ref,
  const T & expected,
  const T & desired)
{
  return ref.compare_exchange(expected, desired);
}

/**
 * Atomically executes specified operation on the referenced shared value.
 */
template<
  typename T,
  typename BinaryOp >
void op(
  const ddm::GlobRef<ddm::Atomic<T>> & ref,
  const BinaryOp  binary_op,
  /// Value to be added to global atomic variable.
  const T & value)
{
  ref.op(binary_op, value);
}

/**
 * Atomic fetch-and-op operation on the referenced shared value.
 *
 * \return  The value of the referenced shared variable before the
 *          operation.
 */
template<
  typename T,
  typename BinaryOp >
T fetch_op(
  const ddm::GlobRef<ddm::Atomic<T>> & ref,
  const BinaryOp  binary_op,
  /// Value to be added to global atomic variable.
  const T & value)
{
  return ref.fetch_op(binary_op, value);
}

/**
 * Atomic add operation on the referenced shared value.
 */
template<typename T>
typename std::enable_if<
  std::is_integral<T>::value,
  void>::type
add(
  const ddm::GlobRef<ddm::Atomic<T>> & ref,
  const T & value)
{
  ref.add(value);
}

/**
 * Atomic subtract operation on the referenced shared value.
 */
template<typename T>
typename std::enable_if<
  std::is_integral<T>::value,
  void>::type
sub(
  const ddm::GlobRef<ddm::Atomic<T>> & ref,
  const T & value)
{
  ref.sub(value);
}

/**
 * Atomic fetch-and-add operation on the referenced shared value.
 *
 * \return  The value of the referenced shared variable before the
 *          operation.
 */
template<typename T>
typename std::enable_if<
  std::is_integral<T>::value,
  T>::type
fetch_add(
  const ddm::GlobRef<ddm::Atomic<T>> & ref,
  /// Value to be added to global atomic variable.
  const T & value)
{
  return ref.fetch_add(value);
}

/**
 * Atomic fetch-and-sub operation on the referenced shared value.
 *
 * \return  The value of the referenced shared variable before the
 *          operation.
 */
template<typename T>
typename std::enable_if<
  std::is_integral<T>::value,
  T>::type
fetch_sub(
  const ddm::GlobRef<ddm::Atomic<T>> & ref,
  /// Value to be subtracted from global atomic variable.
  const T & value)
{
  return ref.fetch_sub(value);
}
  
} // namespace atomic
} // namespace ddm

#endif // DDM__ATOMIC_OPERATION_H_

