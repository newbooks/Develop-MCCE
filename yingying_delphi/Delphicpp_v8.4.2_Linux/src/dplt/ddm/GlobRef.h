#ifndef DDM__GLOBREF_H_
#define DDM__GLOBREF_H_

#include "../ddm/GlobMem.h"
#include "../ddm/Init.h"


namespace ddm {

// Forward declaration
template<typename T, class A> class GlobMem;
// Forward declaration
template<typename T, class PatternT> class GlobPtr;
// Forward declaration
template<typename T, class PatternT>
void put_value(const T & newval, const GlobPtr<T, PatternT> & gptr);
// Forward declaration
template<typename T, class PatternT>
void get_value(T* ptr, const GlobPtr<T, PatternT> & gptr);


template<typename T>
struct has_subscript_operator
{
  typedef char (& yes)[1];
  typedef char (& no)[2];

  template <typename C> static yes check(decltype(&C::operator[]));
  template <typename>   static no  check(...);

  static bool const value = sizeof(check<T>(0)) == sizeof(yes);
};

template<typename T>
class GlobRef
{
  template<typename U>
  friend std::ostream & operator<<(
    std::ostream & os,
    const GlobRef<U> & gref);

  template <
    typename ElementT >
  friend class GlobRef;
  
  typedef typename std::remove_const<T>::type
    nonconst_value_type;
public:
  typedef T                 value_type;

  typedef GlobRef<const T>  const_type;
  
private:
  typedef GlobRef<T>
    self_t;
  typedef GlobRef<const T>
    self_const_t;

private:
  dart_gptr_t _gptr;

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
  template<class PatternT, class ElementT>
  explicit GlobRef(
    /// Pointer to referenced object in global memory
    GlobPtr<ElementT, PatternT> & gptr)
  : GlobRef(gptr.dart_gptr())
  { }

  /**
   * Constructor, creates an GlobRef object referencing an element in global
   * memory.
   */
  template<class PatternT, class ElementT>
  explicit GlobRef(
    /// Pointer to referenced object in global memory
    const GlobPtr<ElementT, PatternT> & gptr)
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
   * TODO: Try deleting copy constructors to preserve unified copy semantics
   *       ref_a = ref_b.
   *
   * Copy constructor.
   */
  GlobRef(const self_t & other) = default;
 
  GlobRef(self_t && other) = default;

  /**
   * TODO: Try deleting copy constructors to preserve unified copy semantics
   *       ref_a = ref_b.
   *
   * Copy constructor.
   */
  template <class ElementT>
  GlobRef(
    const GlobRef<ElementT> & other)
  : _gptr(other._gptr)
  { }

  GlobRef<T> & operator=(const T val) {
    set(val);
    return *this;
  }

  /**
   * Assignment operator.
   */
  GlobRef<T> & operator=(const self_t & other)
  {
    // This results in a dart_put, required for STL algorithms like
    // std::copy to work on global ranges.
    // TODO: Not well-defined:
    //       This violates copy semantics, as
    //         GlobRef(const GlobRef & other)
    //       copies the GlobRef instance while
    //         GlobRef=(const GlobRef & other)
    //       puts the value.
    set(static_cast<T>(other));
    return *this;
  }

  /**
   * Assignment operator.
   */
  template <typename GlobRefOrElementT>
  GlobRef<T> & operator=(GlobRefOrElementT && other)
  {
    // This results in a dart_put, required for STL algorithms like
    // std::copy to work on global ranges.
    // TODO: Not well-defined:
    //       This violates copy semantics, as
    //         GlobRef(const GlobRef & other)
    //       copies the GlobRef instance while
    //         GlobRef=(const GlobRef & other)
    //       puts the value.
    set(std::forward<GlobRefOrElementT>(other));
    return *this;
  }

  operator nonconst_value_type() const {
    DDM_LOG_TRACE("GlobRef.T()", "conversion operator");
    DDM_LOG_TRACE_VAR("GlobRef.T()", _gptr);
    nonconst_value_type t;
    dart_storage_t ds = ddm::dart_storage<T>(1);
    dart_get_blocking(static_cast<void *>(&t), _gptr, ds.nelem, ds.dtype);
    DDM_LOG_TRACE_VAR("GlobRef.T >", _gptr);
    return t;
  }

  template <class GlobRefT>
  constexpr bool operator==(const GlobRefT & other) const noexcept
  {
    return _gptr == other._gptr;
  }

  template <class GlobRefT>
  constexpr bool operator!=(const GlobRefT & other) const noexcept
  {
    return !(*this == other);
  }

  constexpr bool operator==(const nonconst_value_type & value) const noexcept
  {
    return static_cast<T>(*this) == value;
  }

  constexpr bool operator!=(const nonconst_value_type & value) const noexcept
  {
    return !(*this == value);
  }

  friend void swap(GlobRef<T> a, GlobRef<T> b) {
    nonconst_value_type temp = static_cast<nonconst_value_type>(a);
    a = b;
    b = temp;
  }

  void set(const T & val) {
    DDM_LOG_TRACE_VAR("GlobRef.set()", val);
    DDM_LOG_TRACE_VAR("GlobRef.set", _gptr);
    // TODO: Clarify if dart-call can be avoided if
    //       _gptr->is_local()
    dart_storage_t ds = ddm::dart_storage<T>(1);
    dart_put_blocking(
        _gptr, static_cast<const void *>(&val), ds.nelem, ds.dtype);
    DDM_LOG_TRACE_VAR("GlobRef.set >", _gptr);
  }

  nonconst_value_type get() const {
    DDM_LOG_TRACE("T GlobRef.get()", "explicit get");
    DDM_LOG_TRACE_VAR("GlobRef.T()", _gptr);
    nonconst_value_type t;
    dart_storage_t ds = ddm::dart_storage<T>(1);
    dart_get_blocking(static_cast<void *>(&t), _gptr, ds.nelem, ds.dtype);
    return t;
  }

  void get(nonconst_value_type *tptr) const {
    DDM_LOG_TRACE("GlobRef.get(T*)", "explicit get into provided ptr");
    DDM_LOG_TRACE_VAR("GlobRef.T()", _gptr);
    dart_storage_t ds = ddm::dart_storage<T>(1);
    dart_get_blocking(static_cast<void *>(tptr), _gptr, ds.nelem, ds.dtype);
  }

  void get(nonconst_value_type& tref) const {
    DDM_LOG_TRACE("GlobRef.get(T&)", "explicit get into provided ref");
    DDM_LOG_TRACE_VAR("GlobRef.T()", _gptr);
    dart_storage_t ds = ddm::dart_storage<T>(1);
    dart_get_blocking(static_cast<void *>(&tref), _gptr, ds.nelem, ds.dtype);
  }

  void put(nonconst_value_type& tref) const {
    DDM_LOG_TRACE("GlobRef.put(T&)", "explicit put of provided ref");
    DDM_LOG_TRACE_VAR("GlobRef.T()", _gptr);
    dart_storage_t ds = ddm::dart_storage<T>(1);
    dart_put_blocking(_gptr, static_cast<void *>(&tref), ds.nelem, ds.dtype);
  }

  void put(nonconst_value_type* tptr) const {
    DDM_LOG_TRACE("GlobRef.put(T*)", "explicit put of provided ptr");
    DDM_LOG_TRACE_VAR("GlobRef.T()", _gptr);
    dart_storage_t ds = ddm::dart_storage<T>(1);
    dart_put_blocking(_gptr, static_cast<void *>(tptr), ds.nelem, ds.dtype);
  }

  GlobRef<T> & operator+=(const nonconst_value_type& ref) {
  #if 0
    // TODO: Alternative implementation, possibly more efficient:
    T add_val = ref;
    T old_val;
    dart_ret_t result = dart_fetch_and_op(
                          _gptr,
                          reinterpret_cast<void *>(&add_val),
                          reinterpret_cast<void *>(&old_val),
                          ddm::dart_datatype<T>::value,
                          ddm::plus<T>().dart_operation(),
                          ddm::Team::All().dart_id());
    dart_flush(_gptr);
  #else
    nonconst_value_type val  = operator nonconst_value_type();
    val   += ref;
    operator=(val);
  #endif
    return *this;
  }

  GlobRef<T> & operator-=(const nonconst_value_type& ref) {
    nonconst_value_type val  = operator nonconst_value_type();
    val   -= ref;
    operator=(val);
    return *this;
  }

  GlobRef<T> & operator++() {
    nonconst_value_type val = operator nonconst_value_type();
    ++val;
    operator=(val);
    return *this;
  }

  GlobRef<T> operator++(int) {
    GlobRef<T> result = *this;
    nonconst_value_type val = operator nonconst_value_type();
    ++val;
    operator=(val);
    return result;
  }

  GlobRef<T> & operator--() {
    nonconst_value_type val = operator nonconst_value_type();
    --val;
    operator=(val);
    return *this;
  }

  GlobRef<T> operator--(int) {
    GlobRef<T> result = *this;
    nonconst_value_type val = operator nonconst_value_type();
    --val;
    operator=(val);
    return result;
  }

  GlobRef<T> & operator*=(const nonconst_value_type& ref) {
    nonconst_value_type val = operator nonconst_value_type();
    val   *= ref;
    operator=(val);
    return *this;
  }

  GlobRef<T> & operator/=(const nonconst_value_type& ref) {
    nonconst_value_type val = operator nonconst_value_type();
    val   /= ref;
    operator=(val);
    return *this;
  }

  GlobRef<T> & operator^=(const nonconst_value_type& ref) {
    nonconst_value_type val = operator nonconst_value_type();
    val   ^= ref;
    operator=(val);
    return *this;
  }

  dart_gptr_t dart_gptr() const {
    return _gptr;
  }

#if 0
  template<
    typename X=T,
    typename std::enable_if<has_subscript_operator<X>::value, int>::type
      * ptr = nullptr>
  auto operator[](size_t pos) ->
    typename std::result_of<decltype(&T::operator[])(T, size_t)>::type
  {
    nonconst_value_type val = operator nonconst_value_type();
    return val[pos];
  }
#endif

  /**
   * Checks whether the globally referenced element is in
   * the calling unit's local memory.
   */
  bool is_local() const {
    return GlobPtr<T>(_gptr).is_local();
  }

  /**
   * Get a global ref to a member of a certain type at the
   * specified offset
   */
  template<typename MEMTYPE>
  GlobRef<MEMTYPE> member(size_t offs) const {
    dart_gptr_t dartptr = _gptr;
    DDM_ASSERT_RETURNS(
      dart_gptr_incaddr(&dartptr, offs),
      DART_OK);
    GlobPtr<MEMTYPE> gptr(dartptr);
    return GlobRef<MEMTYPE>(gptr);
  }

  /**
   * Get the member via pointer to member
   */
  template<class MEMTYPE, class P=T>
  GlobRef<MEMTYPE> member(
    const MEMTYPE P::*mem) const {
    // TODO: Thaaaat ... looks hacky.
    size_t offs = (size_t) &( reinterpret_cast<P*>(0)->*mem);
    return member<MEMTYPE>(offs);
  }

};

template<typename T>
std::ostream & operator<<(
  std::ostream     & os,
  const GlobRef<T> & gref)
{
  char buf[100];
  sprintf(buf,
          "(%06X|%02X|%04X|%04X|%016lX)",
          gref._gptr.unitid,
          gref._gptr.flags,
          gref._gptr.segid,
          gref._gptr.teamid,
          gref._gptr.addr_or_offs.offset);
  os << "ddm::GlobRef<" << typeid(T).name() << ">" << buf;
  return os;
}

} // namespace ddm

#endif // DDM__GLOBREF_H_
