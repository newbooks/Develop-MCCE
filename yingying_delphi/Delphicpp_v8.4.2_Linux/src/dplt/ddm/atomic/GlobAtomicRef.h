#ifndef DDM__ATOMIC_GLOBREF_H_
#define DDM__ATOMIC_GLOBREF_H_

#include "../../ddm/GlobPtr.h"
#include "../../ddm/algorithm/Operation.h"


namespace ddm {

// forward decls
template<typename T>
class Atomic;

template<typename T>
class Shared;

/**
 * Specialization for atomic values. All atomic operations are 
 * \c const as the \c GlobRef does not own the atomic values.
 */
template<typename T>
class GlobRef<ddm::Atomic<T>>
{
    /*
    Used to force using small data types
  static_assert(ddm::dart_datatype<T>::value != DART_TYPE_UNDEFINED,
    "ddm::GlobRef<Atomic<T>> only valid on integral and floating point types");
    */
    /* Now changed to copyable data types */
  static_assert(std::is_trivially_copyable<T>::value,
      "Element type must be trivially copyable");
  
  template<typename U>
  friend std::ostream & operator<<(
    std::ostream & os,
    const GlobRef<U> & gref);
  
  friend class Shared<T>;
  
public:
  typedef T
    value_type;
  typedef GlobRef<const ddm::Atomic<T>>
    const_type;
  
private:
  typedef ddm::Atomic<T>      atomic_t;
  typedef GlobRef<atomic_t>      self_t;
  
public:
  /**
   * Default constructor, creates an GlobRef object referencing an element in
   * global memory.
   */
  GlobRef()
  : _gptr(DART_GPTR_NULL) {
  }
  
  /**
   * Constructor, creates an GlobRef object referencing an element in global
   * memory.
   */
  template<typename PatternT>
  explicit GlobRef(
    /// Pointer to referenced object in global memory
    GlobPtr<atomic_t, PatternT> & gptr)
  : GlobRef(gptr.dart_gptr())
  { }
  
  /**
   * Constructor, creates an GlobRef object referencing an element in global
   * memory.
   */
  template<typename PatternT>
  GlobRef(
    /// Pointer to referenced object in global memory
    const GlobPtr<atomic_t, PatternT> & gptr)
  : GlobRef(gptr.dart_gptr())
  { }
  
  /**
   * Constructor, creates an GlobRef object referencing an element in global
   * memory.
   */
  explicit GlobRef(dart_gptr_t dart_gptr)
  : _gptr(dart_gptr)
  {
    DDM_LOG_TRACE_VAR("GlobRef(dart_gptr_t)", dart_gptr);
  }
  
  /**
   * Copy constructor.
   */
  GlobRef(
    /// GlobRef instance to copy.
    const GlobRef<atomic_t> & other)
  : _gptr(other._gptr)
  { }
  
  self_t & operator=(const self_t & other) = delete;
  
  inline bool operator==(const self_t & other) const noexcept
  {
    return _gptr == other._gptr;
  }
  
  inline bool operator!=(const self_t & other) const noexcept
  {
    return !(*this == other);
  }
  
  inline bool operator==(const T & value) const = delete;
  inline bool operator!=(const T & value) const = delete;
  
  operator T() const {
    return load();
  }
  
  operator GlobPtr<T>() const {
    DDM_LOG_TRACE("GlobRef.GlobPtr()", "conversion operator");
    DDM_LOG_TRACE_VAR("GlobRef.T()", _gptr);
    return GlobPtr<atomic_t>(_gptr);
  }
  
  dart_gptr_t dart_gptr() const {
    return _gptr;
  }
  
  /**
   * Checks whether the globally referenced element is in
   * the calling unit's local memory.
   */
  bool is_local() const {
    return GlobPtr<T>(_gptr).is_local();
  }
  
  /// atomically assigns value
  GlobRef<atomic_t> operator=(const T & value) const {
    store(value);
    return *this;
  }
  
  /**
   * Set the value of the shared atomic variable.
   */
  void store(const T & value) const
  {
    DDM_LOG_DEBUG_VAR("GlobRef<Atomic>.store()", value);
    DDM_LOG_TRACE_VAR("GlobRef<Atomic>.store",   _gptr);
    dart_ret_t ret = dart_accumulate(
                       _gptr,
                       reinterpret_cast<const char * const>(&value),
                       1,
                       ddm::dart_datatype<T>::value,
                       DART_OP_REPLACE);
    dart_flush(_gptr);
    DDM_ASSERT_EQ(DART_OK, ret, "dart_accumulate failed");
    DDM_LOG_DEBUG("GlobRef<Atomic>.store >");
  }

  /// atomically fetches value
  T load() const
  {
    DDM_LOG_DEBUG("GlobRef<Atomic>.load()");
    DDM_LOG_TRACE_VAR("GlobRef<Atomic>.load", _gptr);
    value_type nothing;
    value_type result;
    dart_ret_t ret = dart_fetch_and_op(
                       _gptr,
                       reinterpret_cast<void * const>(&nothing),
                       reinterpret_cast<void * const>(&result),
                       ddm::dart_datatype<T>::value,
                       DART_OP_NO_OP);
    dart_flush_local(_gptr);
    DDM_ASSERT_EQ(DART_OK, ret, "dart_accumulate failed");
    DDM_LOG_DEBUG_VAR("GlobRef<Atomic>.get >", result);
    return result;
  }
  
  /**
   * Atomically executes specified operation on the referenced shared value.
   */
  template<typename BinaryOp>
  void op(
    BinaryOp  binary_op,
    /// Value to be added to global atomic variable.
    const T & value) const
  {
    DDM_LOG_DEBUG_VAR("GlobRef<Atomic>.op()", value);
    DDM_LOG_TRACE_VAR("GlobRef<Atomic>.op",   _gptr);
    value_type acc = value;
    DDM_LOG_TRACE("GlobRef<Atomic>.op", "dart_accumulate");
    dart_ret_t ret = dart_accumulate(
                       _gptr,
                       reinterpret_cast<char *>(&acc),
                       1,
                       ddm::dart_datatype<T>::value,
                       binary_op.dart_operation());
    dart_flush(_gptr);
    DDM_ASSERT_EQ(DART_OK, ret, "dart_accumulate failed");
    DDM_LOG_DEBUG_VAR("GlobRef<Atomic>.op >", acc);
  }
  
  /**
   * Atomic fetch-and-op operation on the referenced shared value.
   *
   * \return  The value of the referenced shared variable before the
   *          operation.
   */
  template<typename BinaryOp>
  T fetch_op(
    BinaryOp  binary_op,
    /// Value to be added to global atomic variable.
    const T & value) const
  {
    DDM_LOG_DEBUG_VAR("GlobRef<Atomic>.fetch_op()", value);
    DDM_LOG_TRACE_VAR("GlobRef<Atomic>.fetch_op",   _gptr);
    DDM_LOG_TRACE_VAR("GlobRef<Atomic>.fetch_op",   typeid(value).name());
    value_type acc;
    dart_ret_t ret = dart_fetch_and_op(
                       _gptr,
                       reinterpret_cast<const void * const>(&value),
                       reinterpret_cast<void * const>(&acc),
                       ddm::dart_datatype<T>::value,
                       binary_op.dart_operation());
    dart_flush(_gptr);
    DDM_ASSERT_EQ(DART_OK, ret, "dart_fetch_op failed");
    DDM_LOG_DEBUG_VAR("GlobRef<Atomic>.fetch_op >", acc);
    return acc;
  }
  
  /**
   * atomically exchanges value
   */
  T exchange(const T & value) const {
    return fetch_op(ddm::second<T>(), value);
  }
  
  /**
   * Atomically compares the value with the value of expected and if those are
   * bitwise-equal, replaces the former with desired.
   * 
   * \return  True if value is exchanged
   * 
   * \see \c ddm::atomic::compare_exchange
   */
  bool compare_exchange(const T & expected, const T & desired) const {
    static_assert(std::is_integral<T>::value,
      "GlobRef<Atomic>.compare_exchange only valid on integral types!");
    DDM_LOG_DEBUG_VAR("GlobRef<Atomic>.compare_exchange()", desired);
    DDM_LOG_TRACE_VAR("GlobRef<Atomic>.compare_exchange",   _gptr);
    DDM_LOG_TRACE_VAR("GlobRef<Atomic>.compare_exchange",   expected);
    DDM_LOG_TRACE_VAR("GlobRef<Atomic>.compare_exchange",   typeid(desired).name());
    value_type result;
    dart_ret_t ret = dart_compare_and_swap(
                       _gptr,
                       reinterpret_cast<const void * const>(&desired),
                       reinterpret_cast<const void * const>(&expected),
                       reinterpret_cast<void * const>(&result),
                       ddm::dart_datatype<T>::value);
    dart_flush(_gptr);
    DDM_ASSERT_EQ(DART_OK, ret, "dart_compare_and_swap failed");
    DDM_LOG_DEBUG_VAR("GlobRef<Atomic>.compare_exchange >", (expected == result));
    return (expected == result);
  }
  
  /*
   * ---------------------------------------------------------------------------
   * ------------ specializations for atomic integral types --------------------
   * ---------------------------------------------------------------------------
   *
   *  As the check for integral type is already implemented in constructor, 
   *  no check is performed here
   */
  
  /**
   * DDM specific variant which is faster than \cfetch_add
   * but does not return value
   */
  void add (const T & value) const
  {
    op(ddm::plus<T>(), value);
  }
  
  /**
   * Atomic fetch-and-add operation on the referenced shared value.
   *
   * \return  The value of the referenced shared variable before the
   *          operation.
   */
  T fetch_add (
    /// Value to be added to global atomic variable.
    const T & value) const
  {
    return fetch_op(ddm::plus<T>(), value);
  }
  
  /**
   * DDM specific variant which is faster than \cfetch_sub
   * but does not return value
   */
  void sub (const T & value) const
  {
    op(ddm::plus<T>(), -value);
  }
  
  /**
   * Atomic fetch-and-sub operation on the referenced shared value.
   *
   * \return  The value of the referenced shared variable before the
   *          operation.
   */
  T fetch_sub (
    /// Value to be subtracted from global atomic variable.
    const T & value) const
  {
    return fetch_op(ddm::plus<T>(), -value);
  }
  
  /// prefix atomically increment value by one
  T operator++ () const {
    return fetch_add(1) + 1;
  }
  
  /// postfix atomically increment value by one
  T operator++ (int) const {
    return fetch_add(1);
  }
  
  /// prefix atomically decrement value by one
  T operator-- () const {
    return fetch_sub(1) - 1;
  }
  
  /// postfix atomically decrement value by one
  T operator-- (int) const {
    return fetch_sub(1);
  }
  
  /// atomically increment value by ref
  T operator+=(const T & value) const {
    return fetch_add(value) + value;
  }
  
  /// atomically decrement value by ref
  T operator-=(const T & value) const {
    return fetch_sub(value) - value;
  }

private:
  /**
   * Assignment of \c GlobRef<Atomic> is prohibited but in some cases
   * (e.g. \c ddm::Shared) it is necessary and useful.
   * 
   * \note this operation is not atomic
   */
  void _set_dart_gptr(const dart_gptr_t & gptr){
    _gptr = gptr;
  }
  
private:
  dart_gptr_t _gptr;
  
};

} // namespace ddm

#endif // DDM__ATOMIC_GLOBREF_H_

