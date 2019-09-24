#ifndef DDM__ATOMIC_H__INCLUDED
#define DDM__ATOMIC_H__INCLUDED

#include "../ddm/internal/TypeInfo.h"

#include <iostream>


namespace ddm {

/**
 * Type wrapper to mark any trivial type atomic.
 *
 * If one unit writes to an atomic object while another unit reads from it,
 * the behavior is well-defined. The DDM version follows as closely as
 * possible the interface of \c std::atomic
 * However as data has to be transferred between
 * units using DART, the atomicity guarantees are set by the DART
 * implementation.
 *
 * \note \c Atomic objects have to be placed in a DDM container,
 *       and can only be accessed using \c GlobRef<ddm::Atomic<T>> .
 *       Local accesses to atomic elements are not allowed.
 *
 *       \code
 *       ddm::Array<ddm::Atomic<int>> array(100);
 *       array.local[10].load()               // not allowed
 *       ddm::atomic::load(array.local[10])  // not allowed
 *       ddm::atomic::load(array.lbegin())   // not allowed
 *       \endcode
 * \endnote
 * 
 * \code
 *   ddm::Array<ddm::Atomic<int>> array(100);
 *   // supported as Atomic<value_t>(value_t T) is available
 *   ddm::fill(array.begin(), array.end(), 0);
 *
 *   if(ddm::myid() == 0){
 *     array[10].store(5);
 *   }
 *   ddm::barrier();
 *   array[10].add(1);
 *   // postcondition:
 *   // array[10] == ddm::size() + 5
 * \endcode
 */
template<typename T>
class Atomic {
private:
  T _value;
  typedef Atomic<T> self_t;

public:
  typedef T value_type;

  constexpr Atomic()                                 = default;
  constexpr Atomic(const Atomic<T> & other)          = default;
  self_t & operator=(const self_t & other)           = default;

  /**
   * Initializes the underlying value with desired.
   * The initialization is not atomic
   */
  constexpr Atomic(T value)
  : _value(value) { }

  /**
   * Disabled assignment as this violates the atomic semantics
   *
   * TODO: Assignment semantics are not well-defined:
   *       - Constructor Atomic(T)  is default-defined
   *       - Assignment  Atomic=(T) is deleted
   */
  T operator=(T value) = delete;

  /**
   * As \c Atomic is implemented as phantom type,
   * the value has to be queried using the \c ddm::GlobRef
   */
  operator T()         = delete;

  constexpr bool operator==(const self_t & other) const {
    return _value == other._value;
  }

  constexpr bool operator!=(const self_t & other) const {
    return !(*this == other);
  }

  template<typename T_>
  friend std::ostream & operator<<(
    std::ostream     & os,
    const Atomic<T_> & at);

}; // class Atomic

/**
 * type traits for \c ddm::Atomic
 *
 * true if type is atomic
 * false otherwise
 */
template<typename T>
struct is_atomic {
  static constexpr bool value = false;
};
template<typename T>
struct is_atomic<ddm::Atomic<T>> {
  static constexpr bool value = true;
};

template<typename T>
std::ostream & operator<<(
  std::ostream    & os,
  const Atomic<T> & at)
{
  std::ostringstream ss;
  ss << ddm::internal::typestr(at) << "<phantom>";
  return operator<<(os, ss.str());
}

} // namespace ddm

#include "../ddm/atomic/GlobAtomicRef.h"
#include "../ddm/atomic/Operation.h"

#endif // DDM__ATOMIC_H__INCLUDED

