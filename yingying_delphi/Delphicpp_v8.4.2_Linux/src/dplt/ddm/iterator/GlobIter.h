#ifndef DDM__GLOB_ITER_H__INCLUDED
#define DDM__GLOB_ITER_H__INCLUDED

#include "../../ddm/Pattern.h"
#include "../../ddm/GlobRef.h"
#include "../../ddm/GlobPtr.h"

#include <functional>
#include <sstream>

namespace ddm {

#ifndef DOXYGEN

// Forward-declaration
template<
  typename ElementType,
  class    PatternType,
  class    GlobMemType,
  class    PointerType,
  class    ReferenceType >
class GlobStencilIter;
// Forward-declaration
template<
  typename ElementType,
  class    PatternType,
  class    GlobMemType,
  class    PointerType,
  class    ReferenceType >
class GlobViewIter;

#endif // DOXYGEN

/**
 * \defgroup  DDMGlobalIteratorConcept  Global Iterator Concept
 * Concept for iterators in global index space.
 *
 * \ingroup DDMConcept
 * \{
 * \par Description
 *
 * \par Methods
 * Return Type             | Method                 | Parameters     | Description                                             |
 * ----------------------- | ---------------------- | -------------- | ------------------------------------------------------- |
 * <tt>dart_gptr_t</tt>    | <tt>dart_gptr</tt>     | &nbsp;         | DART global pointer on the iterator's current position. |
 *
 * \}
 */


/**
 * Iterator on Partitioned Global Address Space.
 *
 * \concept{DDMGlobalIteratorConcept}
 */
template<
  typename ElementType,
  class    PatternType,
  class    GlobMemType
             = GlobMem< typename std::remove_const<ElementType>::type >,
  class    PointerType   = GlobPtr<ElementType, PatternType>,
  class    ReferenceType = GlobRef<ElementType> >
class GlobIter
: public std::iterator<
           std::random_access_iterator_tag,
           ElementType,
           typename PatternType::index_type,
           PointerType,
           ReferenceType >
{
private:
  typedef GlobIter<
            ElementType,
            PatternType,
            GlobMemType,
            PointerType,
            ReferenceType>
    self_t;

  typedef typename std::remove_const<ElementType>::type
    nonconst_value_type;
public:
  typedef          ElementType                         value_type;

  typedef          ReferenceType                        reference;
  typedef typename ReferenceType::const_type      const_reference;

  typedef          PointerType                            pointer;
  typedef typename PointerType::const_type          const_pointer;

  typedef typename GlobMemType::local_pointer       local_pointer;
  typedef typename GlobMemType::local_pointer          local_type;

  typedef          PatternType                       pattern_type;
  typedef typename PatternType::index_type             index_type;

private:
  typedef GlobIter<
            const ElementType,
            PatternType,
            GlobMemType,
            const_pointer,
            const_reference >
    self_const_t;

public:
  typedef std::integral_constant<bool, false>       has_view;

public:
  // For ostream output
  template <
    typename T_,
    class    P_,
    class    GM_,
    class    Ptr_,
    class    Ref_ >
  friend std::ostream & operator<<(
           std::ostream & os,
           const GlobIter<T_, P_, GM_, Ptr_, Ref_> & it);

  // For conversion to GlobStencilIter
  template<
    typename T_,
    class    P_,
    class    GM_,
    class    Ptr_,
    class    Ref_ >
  friend class GlobStencilIter;

  // For conversion to GlobViewIter
  template<
    typename T_,
    class    P_,
    class    GM_,
    class    Ptr_,
    class    Ref_ >
  friend class GlobViewIter;

  // For comparison operators
  template<
    typename T_,
    class    P_,
    class    GM_,
    class    Ptr_,
    class    Ref_ >
  friend class GlobIter;

private:
  static const dim_t      NumDimensions = PatternType::ndim();
  static const MemArrange Arrangement   = PatternType::memory_order();

protected:
  /// Global memory used to dereference iterated values.
  GlobMemType          * _globmem;
  /// Pattern that specifies the iteration order (access pattern).
  const PatternType    * _pattern;
  /// Current position of the iterator in global canonical index space.
  index_type             _idx             = 0;
  /// Maximum position allowed for this iterator.
  index_type             _max_idx         = 0;
  /// Unit id of the active unit
  team_unit_t            _myid;
  /// Pointer to first element in local memory
  local_pointer          _lbegin          = nullptr;

public:
  /**
   * Default constructor.
   */
  GlobIter()
  : _globmem(nullptr),
    _pattern(nullptr),
    _idx(0),
    _max_idx(0),
    _myid(ddm::Team::All().myid()),
    _lbegin(nullptr)
  {
    DDM_LOG_TRACE_VAR("GlobIter()", _idx);
    DDM_LOG_TRACE_VAR("GlobIter()", _max_idx);
  }

  /**
   * Constructor, creates a global iterator on global memory following
   * the element order specified by the given pattern.
   */
  GlobIter(
    GlobMemType       * gmem,
      const PatternType & pat,
      index_type          position = 0)
  : _globmem(gmem),
    _pattern(&pat),
    _idx(position),
    _max_idx(pat.size() - 1),
    _myid(pat.team().myid()),
    _lbegin(_globmem->lbegin())
  {
    DDM_LOG_TRACE_VAR("GlobIter(gmem,pat,idx,abs)", _idx);
    DDM_LOG_TRACE_VAR("GlobIter(gmem,pat,idx,abs)", _max_idx);
  }

  /**
   * Copy constructor.
   */
  template <class GlobIterT>
  GlobIter(
    const GlobIterT & other)
  : _globmem(other._globmem)
  , _pattern(other._pattern)
  , _idx    (other._idx)
  , _max_idx(other._max_idx)
  , _myid   (other._myid)
  , _lbegin (other._lbegin)
  { }

  /**
   * Assignment operator.
   */
  template <
    typename T_,
    class    P_,
    class    GM_,
    class    Ptr_,
    class    Ref_ >
  self_t & operator=(
    const GlobIter<T_, P_, GM_, Ptr_, Ref_ > & other)
  {
    _globmem = other._globmem;
    _pattern = other._pattern;
    _idx     = other._idx;
    _max_idx = other._max_idx;
    _myid    = other._myid;
    _lbegin  = other._lbegin;
  }

  /**
   * The number of dimensions of the iterator's underlying pattern.
   */
  static dim_t ndim()
  {
    return NumDimensions;
  }

  /**
   * Type conversion operator to \c GlobPtr.
   *
   * \return  A global reference to the element at the iterator's position
   */
  explicit operator const_pointer() const {
    DDM_LOG_TRACE_VAR("GlobIter.GlobPtr()", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    index_type idx    = _idx;
    index_type offset = 0;
    // Convert iterator position (_idx) to local index and unit.
    if (_idx > _max_idx) {
      // Global iterator pointing past the range indexed by the pattern
      // which is the case for .end() iterators.
      idx     = _max_idx;
      offset += _idx - _max_idx;
    }
    // Global index to local index and unit:
    local_pos_t local_pos = _pattern->local(idx);
    DDM_LOG_TRACE_VAR("GlobIter.GlobPtr >", local_pos.unit);
    DDM_LOG_TRACE_VAR("GlobIter.GlobPtr >", local_pos.index);
    // Create global pointer from unit and local offset:
    const_pointer gptr(
      _globmem->at(team_unit_t(local_pos.unit), local_pos.index)
    );
    return gptr + offset;
  }

  /**
   * Type conversion operator to \c GlobPtr.
   *
   * \return  A global reference to the element at the iterator's position
   */
  explicit operator pointer() {
    DDM_LOG_TRACE_VAR("GlobIter.GlobPtr()", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    index_type idx    = _idx;
    index_type offset = 0;
    // Convert iterator position (_idx) to local index and unit.
    if (_idx > _max_idx) {
      // Global iterator pointing past the range indexed by the pattern
      // which is the case for .end() iterators.
      idx     = _max_idx;
      offset += _idx - _max_idx;
    }
    // Global index to local index and unit:
    local_pos_t local_pos = _pattern->local(idx);
    DDM_LOG_TRACE_VAR("GlobIter.GlobPtr >", local_pos.unit);
    DDM_LOG_TRACE_VAR("GlobIter.GlobPtr >", local_pos.index);
    // Create global pointer from unit and local offset:
    pointer gptr(
      _globmem->at(team_unit_t(local_pos.unit), local_pos.index)
    );
    return gptr + offset;
  }

  /**
   * Explicit conversion to \c dart_gptr_t.
   *
   * \return  A DART global pointer to the element at the iterator's
   *          position
   */
  dart_gptr_t dart_gptr() const
  {
    DDM_LOG_TRACE_VAR("GlobIter.dart_gptr()", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    index_type idx    = _idx;
    index_type offset = 0;
    // Convert iterator position (_idx) to local index and unit.
    if (_idx > _max_idx) {
      // Global iterator pointing past the range indexed by the pattern
      // which is the case for .end() iterators.
      idx     = _max_idx;
      offset += _idx - _max_idx;
      DDM_LOG_TRACE_VAR("GlobIter.dart_gptr", _max_idx);
      DDM_LOG_TRACE_VAR("GlobIter.dart_gptr", idx);
      DDM_LOG_TRACE_VAR("GlobIter.dart_gptr", offset);
    }
    // Global index to local index and unit:
    local_pos_t local_pos = _pattern->local(idx);
    DDM_LOG_TRACE("GlobIter.dart_gptr",
                   "unit:",        local_pos.unit,
                   "local index:", local_pos.index);
    // Global pointer to element at given position:
    const_pointer gptr(
      _globmem->at(
        team_unit_t(local_pos.unit),
        local_pos.index)
    );
    DDM_LOG_TRACE_VAR("GlobIter.dart_gptr >", gptr);
    return (gptr + offset).dart_gptr();
  }

  /**
   * Dereference operator.
   *
   * \return  A global reference to the element at the iterator's position.
   */
  reference operator*()
  {
    DDM_LOG_TRACE("GlobIter.*()", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    index_type idx = _idx;
    // Global index to local index and unit:
    local_pos_t local_pos = _pattern->local(idx);
    DDM_LOG_TRACE("GlobIter.* >",
                   "unit:", local_pos.unit, "index:", local_pos.index);
    // Global reference to element at given position:
    return reference(
             _globmem->at(local_pos.unit,
                          local_pos.index));
  }

  /**
   * Dereference operator.
   *
   * \return  A global reference to the element at the iterator's position.
   */
  const_reference operator*() const
  {
    DDM_LOG_TRACE("GlobIter.*", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    index_type idx = _idx;
    // Global index to local index and unit:
    local_pos_t local_pos = _pattern->local(idx);
    DDM_LOG_TRACE_VAR("GlobIter.*", local_pos.unit);
    DDM_LOG_TRACE_VAR("GlobIter.*", local_pos.index);
    // Global reference to element at given position:
    return const_reference(
             _globmem->at(local_pos.unit,
                          local_pos.index));
  }

  /**
   * Subscript operator, returns global reference to element at given
   * global index.
   */
  reference operator[](
    /// The global position of the element
    index_type g_index)
  {
    DDM_LOG_TRACE("GlobIter.[]", g_index);
    index_type idx = g_index;
    typedef typename pattern_type::local_index_t
      local_pos_t;
    // Global index to local index and unit:
    local_pos_t local_pos = _pattern->local(idx);
    DDM_LOG_TRACE_VAR("GlobIter.[]", local_pos.unit);
    DDM_LOG_TRACE_VAR("GlobIter.[]", local_pos.index);
    // Global reference to element at given position:
    return reference(
             _globmem->at(local_pos.unit,
                          local_pos.index));
  }

  /**
   * Subscript operator, returns global reference to element at given
   * global index.
   */
  const_reference operator[](
    /// The global position of the element
    index_type g_index) const
  {
    DDM_LOG_TRACE("GlobIter.[]", g_index);
    index_type idx = g_index;
    typedef typename pattern_type::local_index_t
      local_pos_t;
    // Global index to local index and unit:
    local_pos_t local_pos = _pattern->local(idx);
    DDM_LOG_TRACE_VAR("GlobIter.[]", local_pos.unit);
    DDM_LOG_TRACE_VAR("GlobIter.[]", local_pos.index);
    // Global reference to element at given position:
    return const_reference(
             _globmem->at(local_pos.unit,
                          local_pos.index));
  }

  /**
   * Checks whether the element referenced by this global iterator is in
   * the calling unit's local memory.
   */
  constexpr bool is_local() const
  {
    return (_myid == lpos().unit);
  }

  /**
   * Convert global iterator to native pointer.
   *
   * TODO: Evaluate alternative:
   *         auto l_idx_this = _container.pattern().local(this->pos());
   *         return (l_idx_this.unit == _myid
   *                 ? _lbegin + l_idx_this
   *                 : nullptr
   *                );
   */
  local_pointer local() const
  {
    DDM_LOG_TRACE_VAR("GlobIter.local=()", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    index_type idx    = _idx;
    index_type offset = 0;
    DDM_LOG_TRACE_VAR("GlobIter.local=", _max_idx);
    // Convert iterator position (_idx) to local index and unit.
    if (_idx > _max_idx) {
      // Global iterator pointing past the range indexed by the pattern
      // which is the case for .end() iterators.
      idx     = _max_idx;
      offset += _idx - _max_idx;
    }
    DDM_LOG_TRACE_VAR("GlobIter.local=", idx);
    DDM_LOG_TRACE_VAR("GlobIter.local=", offset);
    // Global index to local index and unit:
    local_pos_t local_pos = _pattern->local(idx);
    DDM_LOG_TRACE_VAR("GlobIter.local= >", local_pos.unit);
    DDM_LOG_TRACE_VAR("GlobIter.local= >", local_pos.index);
    if (_myid != local_pos.unit) {
      // Iterator position does not point to local element
      return nullptr;
    }
    return (_lbegin + local_pos.index + offset);
  }

  /**
   * Unit and local offset at the iterator's position.
   */
  inline typename pattern_type::local_index_t lpos() const
  {
    DDM_LOG_TRACE_VAR("GlobIter.lpos()", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    index_type idx    = _idx;
    index_type offset = 0;
    // Convert iterator position (_idx) to local index and unit.
    if (_idx > _max_idx) {
      // Global iterator pointing past the range indexed by the pattern
      // which is the case for .end() iterators.
      idx    = _max_idx;
      offset = _idx - _max_idx;
      DDM_LOG_TRACE_VAR("GlobIter.lpos", _max_idx);
      DDM_LOG_TRACE_VAR("GlobIter.lpos", idx);
      DDM_LOG_TRACE_VAR("GlobIter.lpos", offset);
    }
    // Global index to local index and unit:
    local_pos_t local_pos = _pattern->local(idx);
    local_pos.index += offset;
    DDM_LOG_TRACE("GlobIter.lpos >",
                   "unit:",        local_pos.unit,
                   "local index:", local_pos.index);
    return local_pos;
  }

  /**
   * Map iterator to global index domain.
   */
  constexpr const self_t & global() const noexcept {
    return *this;
  }

  /**
   * Map iterator to global index domain.
   */
  self_t & global() {
    return *this;
  }

  /**
   * Position of the iterator in global index space.
   */
  constexpr index_type pos() const noexcept
  {
    return _idx;
  }

  /**
   * Position of the iterator in global index range.
   */
  constexpr index_type gpos() const noexcept
  {
    return _idx;
  }

  /**
   * Whether the iterator's position is relative to a view.
   *
   * TODO:
   * should be iterator trait:
   *   ddm::iterator_traits<GlobIter<..>>::is_relative()::value
   */
  constexpr bool is_relative() const noexcept
  {
    return false;
  }

  /**
   * The instance of \c GlobMem used by this iterator to resolve addresses
   * in global memory.
   */
  constexpr const GlobMemType & globmem() const noexcept
  {
    return *_globmem;
  }

  /**
   * The instance of \c GlobMem used by this iterator to resolve addresses
   * in global memory.
   */
  inline GlobMemType & globmem()
  {
    return *_globmem;
  }

  /**
   * Prefix increment operator.
   */
  inline self_t & operator++()
  {
    ++_idx;
    return *this;
  }

  /**
   * Postfix increment operator.
   */
  inline self_t operator++(int)
  {
    self_t result = *this;
    ++_idx;
    return result;
  }

  /**
   * Prefix decrement operator.
   */
  inline self_t & operator--()
  {
    --_idx;
    return *this;
  }

  /**
   * Postfix decrement operator.
   */
  inline self_t operator--(int)
  {
    self_t result = *this;
    --_idx;
    return result;
  }

  inline self_t & operator+=(index_type n)
  {
    _idx += n;
    return *this;
  }

  inline self_t & operator-=(index_type n)
  {
    _idx -= n;
    return *this;
  }

  constexpr self_t operator+(index_type n) const noexcept
  {
    return self_t(
      _globmem,
      *_pattern,
      _idx + static_cast<index_type>(n));
  }

  constexpr self_t operator-(index_type n) const noexcept
  {
    return self_t(
      _globmem,
      *_pattern,
      _idx - static_cast<index_type>(n));
  }

  template <class GlobIterT>
  constexpr auto operator+(
    const GlobIterT & other) const noexcept
    -> typename std::enable_if<
         !std::is_integral<GlobIterT>::value,
         index_type
       >::type
  {
    return _idx + other._idx;
  }

  template <class GlobIterT>
  constexpr auto operator-(
    const GlobIterT & other) const noexcept
    -> typename std::enable_if<
         !std::is_integral<GlobIterT>::value,
         index_type
       >::type
  {
    return _idx - other._idx;
  }

  template <class GlobIterT>
  constexpr bool operator<(const GlobIterT & other) const noexcept
  {
    return (_idx < other._idx);
  }

  template <class GlobIterT>
  constexpr bool operator<=(const GlobIterT & other) const noexcept
  {
    return (_idx <= other._idx);
  }

  template <class GlobIterT>
  constexpr bool operator>(const GlobIterT & other) const noexcept
  {
    return (_idx > other._idx);
  }

  template <class GlobIterT>
  constexpr bool operator>=(const GlobIterT & other) const noexcept
  {
    return (_idx >= other._idx);
  }

  template <class GlobIterT>
  constexpr bool operator==(const GlobIterT & other) const noexcept
  {
    return _idx == other._idx;
  }

  template <class GlobIterT>
  constexpr bool operator!=(const GlobIterT & other) const noexcept
  {
    return _idx != other._idx;
  }

  constexpr const PatternType & pattern() const noexcept
  {
    return *_pattern;
  }

  constexpr ddm::Team & team() const noexcept
  {
    return _pattern->team();
  }

}; // class GlobIter


template <
  typename ElementType,
  class    Pattern,
  class    GlobMem,
  class    Pointer,
  class    Reference >
std::ostream & operator<<(
  std::ostream & os,
  const ddm::GlobIter<
          ElementType, Pattern, GlobMem, Pointer, Reference> & it)
{
  std::ostringstream ss;
  ddm::GlobPtr<const ElementType, Pattern> ptr(it.dart_gptr());
  ss << "ddm::GlobIter<" << typeid(ElementType).name() << ">("
     << "idx:"  << it._idx << ", "
     << "gptr:" << ptr << ")";
  return operator<<(os, ss.str());
}

} // namespace ddm

#include "../../ddm/iterator/GlobViewIter.h"

#endif // DDM__GLOB_ITER_H__INCLUDED
