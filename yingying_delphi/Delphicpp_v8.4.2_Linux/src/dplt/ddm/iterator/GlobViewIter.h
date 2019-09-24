#ifndef DDM__GLOB_VIEW_ITER_H__INCLUDED
#define DDM__GLOB_VIEW_ITER_H__INCLUDED

#include "../../ddm/Pattern.h"
#include "../../ddm/Allocator.h"
#include "../../ddm/GlobRef.h"
#include "../../ddm/GlobPtr.h"
#include "../../ddm/GlobMem.h"

#include "../../ddm/iterator/GlobIter.h"

#include <cstddef>
#include <functional>
#include <sstream>
#include <iostream>

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
#endif

/**
 * Global iterator on an index set specified by a view modifier.
 *
 * \concept{DDMGlobalIteratorConcept}
 */
template<
  typename ElementType,
  class    PatternType,
  class    GlobMemType   = GlobMem<ElementType>,
  class    PointerType   = GlobPtr<ElementType, PatternType>,
  class    ReferenceType = GlobRef<ElementType> >
class GlobViewIter
: public std::iterator<
           std::random_access_iterator_tag,
           ElementType,
           typename PatternType::index_type,
           PointerType,
           ReferenceType > {
private:
  typedef GlobViewIter<
            ElementType,
            PatternType,
            GlobMemType,
            PointerType,
            ReferenceType>
    self_t;

  typedef typename PatternType::viewspec_type
    ViewSpecType;
  typedef typename PatternType::index_type
    IndexType;

public:
  typedef ElementType                             value_type;

  typedef ReferenceType                            reference;
  typedef typename reference::const_type     const_reference;

  typedef PointerType                                pointer;
  typedef typename pointer::const_type         const_pointer;

  typedef PatternType                           pattern_type;
  typedef IndexType                               index_type;

private:
  typedef GlobViewIter<
            const ElementType,
            PatternType,
            GlobMemType,
            const_pointer,
            const_reference >
    self_const_t;

public:
  typedef std::integral_constant<bool, true>        has_view;

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
      const GlobViewIter<T_, P_, GM_, Ptr_, Ref_> & it);

  // For conversion to GlobStencilIter
  template<
    typename T_,
    class    P_,
    class    Ptr_,
    class    GM_,
    class    Ref_ >
  friend class GlobStencilIter;

  template<
    typename T_,
    class    P_,
    class    GM_,
    class    Ptr_,
    class    Ref_ >
  friend class GlobViewIter;

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
  GlobMemType                * _globmem;
  /// Pattern that specifies the iteration order (access pattern).
  const PatternType          * _pattern;
  /// View that specifies the iterator's index range relative to the global
  /// index range of the iterator's pattern.
  const ViewSpecType         * _viewspec;
  /// Current position of the iterator relative to the iterator's view.
  IndexType                    _idx             = 0;
  /// The iterator's view index start offset.
  IndexType                    _view_idx_offset = 0;
  /// Maximum position relative to the viewspec allowed for this iterator.
  IndexType                    _max_idx         = 0;
  /// Unit id of the active unit
  team_unit_t                  _myid;
  /// Pointer to first element in local memory
  ElementType                * _lbegin          = nullptr;

public:
  /**
   * Default constructor.
   */
  GlobViewIter()
  : _globmem(nullptr),
    _pattern(nullptr),
    _viewspec(nullptr),
    _idx(0),
    _view_idx_offset(0),
    _max_idx(0),
    _myid(ddm::Team::GlobalUnitID()),
    _lbegin(nullptr)
  {
    DDM_LOG_TRACE_VAR("GlobViewIter()", _idx);
    DDM_LOG_TRACE_VAR("GlobViewIter()", _max_idx);
  }

  /**
   * Constructor, creates a global iterator on global memory following
   * the element order specified by the given pattern and view spec.
   */
  template <class GlobMemType_>
  GlobViewIter(
    GlobMemType_         * gmem,
      const PatternType    & pat,
    const ViewSpecType   & viewspec,
      IndexType              position          = 0,
    IndexType              view_index_offset = 0)
  : _globmem(reinterpret_cast<GlobMemType *>(gmem)),
    _pattern(&pat),
    _viewspec(&viewspec),
    _idx(position),
    _view_idx_offset(view_index_offset),
    _max_idx(viewspec.size() - 1),
    _myid(ddm::Team::GlobalUnitID()),
    _lbegin(_globmem->lbegin())
  {
    DDM_LOG_TRACE_VAR("GlobViewIter(gmem,pat,vs,idx,abs)", _idx);
    DDM_LOG_TRACE_VAR("GlobViewIter(gmem,pat,vs,idx,abs)", _max_idx);
    DDM_LOG_TRACE_VAR("GlobViewIter(gmem,pat,vs,idx,abs)", *_viewspec);
    DDM_LOG_TRACE_VAR("GlobViewIter(gmem,pat,vs,idx,abs)", _view_idx_offset);
  }

  /**
   * Constructor, creates a global iterator on global memory following
   * the element order specified by the given pattern and view spec.
   */
  template <class GlobMemType_>
  GlobViewIter(
    GlobMemType_      * gmem,
      const PatternType & pat,
      IndexType           position          = 0,
    IndexType           view_index_offset = 0)
  : _globmem(reinterpret_cast<GlobMemType *>(gmem)),
    _pattern(&pat),
    _viewspec(nullptr),
    _idx(position),
    _view_idx_offset(view_index_offset),
    _max_idx(pat.size() - 1),
    _myid(ddm::Team::GlobalUnitID()),
    _lbegin(_globmem->lbegin())
  {
    DDM_LOG_TRACE_VAR("GlobViewIter(gmem,pat,idx,abs)", _idx);
    DDM_LOG_TRACE_VAR("GlobViewIter(gmem,pat,idx,abs)", _max_idx);
    DDM_LOG_TRACE_VAR("GlobViewIter(gmem,pat,idx,abs)", _view_idx_offset);
  }

  /**
   * Constructor, creates a global view iterator from a global iterator.
   */
  template<class PtrT, class RefT>
  GlobViewIter(
    const GlobIter<ElementType, PatternType, PtrT, RefT> & other,
    const ViewSpecType                                   & viewspec,
    IndexType                                              view_idx_offset = 0)
  : _globmem(other._globmem),
    _pattern(other._pattern),
    _viewspec(&viewspec),
    _idx(other._idx),
    _view_idx_offset(view_idx_offset),
    _max_idx(other._max_idx),
    _myid(other._myid),
    _lbegin(other._lbegin)
  {
    DDM_LOG_TRACE_VAR("GlobViewIter(GlobIter)", _idx);
    DDM_LOG_TRACE_VAR("GlobViewIter(GlobIter)", _max_idx);
    DDM_LOG_TRACE_VAR("GlobViewIter(GlobIter)", *_viewspec);
    DDM_LOG_TRACE_VAR("GlobViewIter(GlobIter)", _view_idx_offset);
  }

  /**
   * Copy constructor.
   */
  template <class GlobViewIterT>
  GlobViewIter(
    const GlobViewIterT & other)
  : _globmem        (other._globmem)
  , _pattern        (other._pattern)
  , _viewspec       (other._viewspec)
  , _idx            (other._idx)
  , _view_idx_offset(other._view_idx_offset)
  , _max_idx        (other._max_idx)
  , _myid           (other._myid)
  , _lbegin         (other._lbegin)
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
    const GlobViewIter<T_, P_, GM_, Ptr_, Ref_ > & other)
  {
    _globmem         = other._globmem;
    _pattern         = other._pattern;
    _viewspec        = other._viewspec;
    _idx             = other._idx;
    _view_idx_offset = other._view_idx_offset;
    _max_idx         = other._max_idx;
    _myid            = other._myid;
    _lbegin          = other._lbegin;
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
  operator pointer() const
  {
    DDM_LOG_TRACE_VAR("GlobViewIter.GlobPtr()", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    IndexType idx    = _idx;
    IndexType offset = 0;
    DDM_LOG_TRACE_VAR("GlobViewIter.GlobPtr", _max_idx);
    // Convert iterator position (_idx) to local index and unit.
    if (_idx > _max_idx) {
      // Global iterator pointing past the range indexed by the pattern
      // which is the case for .end() iterators.
      idx     = _max_idx;
      offset += _idx - _max_idx;
    }
    DDM_LOG_TRACE_VAR("GlobViewIter.GlobPtr", idx);
    DDM_LOG_TRACE_VAR("GlobViewIter.GlobPtr", offset);
    // Global index to local index and unit:
    local_pos_t local_pos;
    if (_viewspec == nullptr) {
      // No viewspec mapping required:
      local_pos        = _pattern->local(idx);
    } else {
      // Viewspec projection required:
      auto glob_coords = coords(idx);
      local_pos        = _pattern->local_index(glob_coords);
    }
    DDM_LOG_TRACE_VAR("GlobViewIter.GlobPtr >", local_pos.unit);
    DDM_LOG_TRACE_VAR("GlobViewIter.GlobPtr >", local_pos.index + offset);
    // Create global pointer from unit and local offset:
    PointerType gptr(
      _globmem->at(local_pos.unit, local_pos.index)
    );
    gptr += offset;
    return gptr;
  }

  /**
   * Explicit conversion to \c dart_gptr_t.
   *
   * \return  A DART global pointer to the element at the iterator's
   *          position
   */
  dart_gptr_t dart_gptr() const
  {
    DDM_LOG_TRACE_VAR("GlobViewIter.dart_gptr()", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    IndexType idx    = _idx;
    IndexType offset = 0;
    // Convert iterator position (_idx) to local index and unit.
    if (_idx > _max_idx) {
      // Global iterator pointing past the range indexed by the pattern
      // which is the case for .end() iterators.
      idx     = _max_idx;
      offset += _idx - _max_idx;
      DDM_LOG_TRACE_VAR("GlobViewIter.dart_gptr", _max_idx);
      DDM_LOG_TRACE_VAR("GlobViewIter.dart_gptr", idx);
      DDM_LOG_TRACE_VAR("GlobViewIter.dart_gptr", offset);
    }
    // Global index to local index and unit:
    local_pos_t local_pos;
    if (_viewspec == nullptr) {
      // No viewspec mapping required:
      local_pos        = _pattern->local(idx);
    } else {
      // Viewspec projection required:
      DDM_LOG_TRACE_VAR("GlobViewIter.dart_gptr", _viewspec);
      auto glob_coords = coords(idx);
      local_pos        = _pattern->local_index(glob_coords);
    }
    DDM_LOG_TRACE("GlobViewIter.dart_gptr",
                   "unit:",        local_pos.unit,
                   "local index:", local_pos.index);
    // Global pointer to element at given position:
    ddm::GlobPtr<ElementType, PatternType> gptr(
      _globmem->at(
        local_pos.unit,
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
  reference operator*() const
  {
    DDM_LOG_TRACE("GlobViewIter.*", _idx, _view_idx_offset);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    IndexType idx = _idx;
    // Global index to local index and unit:
    local_pos_t local_pos;
    if (_viewspec == nullptr) {
      // No viewspec mapping required:
      local_pos        = _pattern->local(idx);
    } else {
      // Viewspec projection required:
      auto glob_coords = coords(idx);
      local_pos        = _pattern->local_index(glob_coords);
    }
    DDM_LOG_TRACE_VAR("GlobViewIter.*", local_pos.unit);
    DDM_LOG_TRACE_VAR("GlobViewIter.*", local_pos.index);
    // Global reference to element at given position:
    return ReferenceType(
             _globmem->at(local_pos.unit,
                          local_pos.index));
  }

  /**
   * Subscript operator, returns global reference to element at given
   * global index.
   */
  reference operator[](
    /// The global position of the element
    index_type g_index) const
  {
    DDM_LOG_TRACE("GlobViewIter.[]", g_index, _view_idx_offset);
    IndexType idx = g_index;
    typedef typename pattern_type::local_index_t
      local_pos_t;
    // Global index to local index and unit:
    local_pos_t local_pos;
    if (_viewspec == nullptr) {
      // No viewspec mapping required:
      local_pos        = _pattern->local(idx);
    } else {
      // Viewspec projection required:
      auto glob_coords = coords(idx);
      local_pos        = _pattern->local_index(glob_coords);
    }
    DDM_LOG_TRACE_VAR("GlobViewIter.[]", local_pos.unit);
    DDM_LOG_TRACE_VAR("GlobViewIter.[]", local_pos.index);
    // Global reference to element at given position:
    return ReferenceType(
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
   */
  ElementType * local() const
  {
    DDM_LOG_TRACE_VAR("GlobViewIter.local=()", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    IndexType idx    = _idx;
    IndexType offset = 0;
    DDM_LOG_TRACE_VAR("GlobViewIter.local=", _max_idx);
    // Convert iterator position (_idx) to local index and unit.
    if (_idx > _max_idx) {
      // Global iterator pointing past the range indexed by the pattern
      // which is the case for .end() iterators.
      idx     = _max_idx;
      offset += _idx - _max_idx;
    }
    DDM_LOG_TRACE_VAR("GlobViewIter.local=", idx);
    DDM_LOG_TRACE_VAR("GlobViewIter.local=", offset);

    // Global index to local index and unit:
    local_pos_t local_pos;
    if (_viewspec == nullptr) {
      // No viewspec mapping required:
      local_pos        = _pattern->local(idx);
    } else {
      DDM_LOG_TRACE_VAR("GlobViewIter.local=", *_viewspec);
      // Viewspec projection required:
      auto glob_coords = coords(idx);
      local_pos        = _pattern->local_index(glob_coords);
    }
    DDM_LOG_TRACE_VAR("GlobViewIter.local= >", local_pos.unit);
    DDM_LOG_TRACE_VAR("GlobViewIter.local= >", local_pos.index);
    if (_myid != local_pos.unit) {
      // Iterator position does not point to local element
      return nullptr;
    }
    return (_lbegin + local_pos.index + offset);
  }

  /**
   * Map iterator to global index domain by projecting the iterator's view.
   */
  inline GlobIter<ElementType, PatternType, GlobMemType> global() const
  {
    auto g_idx = gpos();
    return ddm::GlobIter<ElementType, PatternType, GlobMemType>(
             _globmem,
             *_pattern,
             g_idx
           );
  }

  /**
   * Position of the iterator in its view's iteration space and the view's
   * offset in global index space.
   */
  inline index_type pos() const
  {
    DDM_LOG_TRACE("GlobViewIter.pos()",
                   "idx:", _idx, "vs_offset:", _view_idx_offset);
    return _idx + _view_idx_offset;
  }

  /**
   * Position of the iterator in its view's iteration space, disregarding
   * the view's offset in global index space.
   */
  inline index_type rpos() const
  {
    return _idx;
  }

  /**
   * Position of the iterator in global index range.
   * Projects iterator position from its view spec to global index domain.
   */
  inline index_type gpos() const
  {
    DDM_LOG_TRACE_VAR("GlobViewIter.gpos()", _idx);
    if (_viewspec == nullptr) {
      // No viewspec mapping required:
      DDM_LOG_TRACE_VAR("GlobViewIter.gpos >", _idx);
      return _idx;
    } else {
      IndexType idx    = _idx;
      IndexType offset = 0;
      DDM_LOG_TRACE_VAR("GlobViewIter.gpos", *_viewspec);
      DDM_LOG_TRACE_VAR("GlobViewIter.gpos", _max_idx);
      // Convert iterator position (_idx) to local index and unit.
      if (_idx > _max_idx) {
        // Global iterator pointing past the range indexed by the pattern
        // which is the case for .end() iterators.
        idx    = _max_idx;
        offset = _idx - _max_idx;
      }
      // Viewspec projection required:
      auto g_coords = coords(idx);
      DDM_LOG_TRACE_VAR("GlobViewIter.gpos", _idx);
      DDM_LOG_TRACE_VAR("GlobViewIter.gpos", g_coords);
      auto g_idx    = _pattern->memory_layout().at(g_coords);
      DDM_LOG_TRACE_VAR("GlobViewIter.gpos", g_idx);
      g_idx += offset;
      DDM_LOG_TRACE_VAR("GlobViewIter.gpos >", g_idx);
      return g_idx;
    }
  }

  /**
   * Unit and local offset at the iterator's position.
   * Projects iterator position from its view spec to global index domain.
   */
  inline typename pattern_type::local_index_t lpos() const
  {
    DDM_LOG_TRACE_VAR("GlobViewIter.lpos()", _idx);
    typedef typename pattern_type::local_index_t
      local_pos_t;
    IndexType idx    = _idx;
    IndexType offset = 0;
    // Convert iterator position (_idx) to local index and unit.
    if (_idx > _max_idx) {
      // Global iterator pointing past the range indexed by the pattern
      // which is the case for .end() iterators.
      idx    = _max_idx;
      offset = _idx - _max_idx;
      DDM_LOG_TRACE_VAR("GlobViewIter.lpos", _max_idx);
      DDM_LOG_TRACE_VAR("GlobViewIter.lpos", idx);
      DDM_LOG_TRACE_VAR("GlobViewIter.lpos", offset);
    }
    // Global index to local index and unit:
    local_pos_t local_pos;
    if (_viewspec == nullptr) {
      // No viewspec mapping required:
      local_pos        = _pattern->local(idx);
    } else {
      // Viewspec projection required:
      DDM_LOG_TRACE_VAR("GlobViewIter.lpos", *_viewspec);
      auto glob_coords = coords(idx);
      local_pos        = _pattern->local_index(glob_coords);
    }
    local_pos.index += offset;
    DDM_LOG_TRACE("GlobViewIter.lpos >",
                   "unit:",        local_pos.unit,
                   "local index:", local_pos.index);
    return local_pos;
  }

  /**
   * Whether the iterator's position is relative to a view.
   */
  inline bool is_relative() const noexcept
  {
    return _viewspec != nullptr;
  }

  /**
   * The view that specifies this iterator's index range.
   */
  inline ViewSpecType viewspec() const
  {
    if (_viewspec != nullptr) {
      return *_viewspec;
    }
    return ViewSpecType(_pattern->memory_layout().extents());
  }

  /**
   * The instance of \c GlobMem used by this iterator to resolve addresses
   * in global memory.
   */
  constexpr const GlobMemType & globmem() const
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
  self_t & operator++()
  {
    ++_idx;
    return *this;
  }

  /**
   * Postfix increment operator.
   */
  self_t operator++(int)
  {
    self_t result = *this;
    ++_idx;
    return result;
  }

  /**
   * Prefix decrement operator.
   */
  self_t & operator--()
  {
    --_idx;
    return *this;
  }

  /**
   * Postfix decrement operator.
   */
  self_t operator--(int)
  {
    self_t result = *this;
    --_idx;
    return result;
  }

  self_t & operator+=(index_type n)
  {
    _idx += n;
    return *this;
  }

  self_t & operator-=(index_type n)
  {
    _idx -= n;
    return *this;
  }

  self_t operator+(index_type n) const
  {
    if (_viewspec == nullptr) {
      self_t res(
        _globmem,
        *_pattern,
        _idx + static_cast<IndexType>(n),
        _view_idx_offset);
      return res;
    }
    self_t res(
      _globmem,
      *_pattern,
      *_viewspec,
      _idx + static_cast<IndexType>(n),
      _view_idx_offset);
    return res;
  }

  self_t operator-(index_type n) const
  {
    if (_viewspec == nullptr) {
      self_t res(
        _globmem,
        *_pattern,
        _idx - static_cast<IndexType>(n),
        _view_idx_offset);
      return res;
    }
    self_t res(
      _globmem,
      *_pattern,
      *_viewspec,
      _idx - static_cast<IndexType>(n),
      _view_idx_offset);
    return res;
  }

  index_type operator+(
    const self_t & other) const
  {
    return _idx + other._idx;
  }

  index_type operator-(
    const self_t & other) const
  {
    return _idx - other._idx;
  }

  inline bool operator<(const self_t & other) const
  {
    // NOTE:
    // This function call is significantly slower than the explicit
    // implementation in operator== and operator!=.
    return compare(other,
                   std::less<index_type>(),
                   std::less<pointer>());
  }

  inline bool operator<=(const self_t & other) const
  {
    // NOTE:
    // This function call is significantly slower than the explicit
    // implementation in operator== and operator!=.
    return compare(other,
                   std::less_equal<index_type>(),
                   std::less_equal<pointer>());
  }

  inline bool operator>(const self_t & other) const
  {
    // NOTE:
    // This function call is significantly slower than the explicit
    // implementation in operator== and operator!=.
    return compare(other,
                   std::greater<index_type>(),
                   std::greater<pointer>());
  }

  inline bool operator>=(const self_t & other) const
  {
    // NOTE:
    // This function call is significantly slower than the explicit
    // implementation in operator== and operator!=.
    return compare(other,
                   std::greater_equal<index_type>(),
                   std::greater_equal<pointer>());
  }

  inline bool operator==(const self_t & other) const
  {
    // NOTE: See comments in method compare().
    if (_viewspec == other._viewspec) {
      // Same viewspec pointer
      return _idx == other._idx;
    }
    if ((_viewspec != nullptr && other._viewspec != nullptr) &&
        (*_viewspec) == *(other._viewspec)) {
      // Viewspec instances are equal
      return _idx == other._idx;
    }
    auto lhs_local = lpos();
    auto rhs_local = other.lpos();
    return (lhs_local.unit  == rhs_local.unit &&
            lhs_local.index == rhs_local.index);
  }

  inline bool operator!=(const self_t & other) const
  {
    // NOTE: See comments in method compare().
    if (_viewspec == other._viewspec) {
      // Same viewspec pointer
      return _idx != other._idx;
    }
    if ((_viewspec != nullptr && other._viewspec != nullptr) &&
        (*_viewspec) == *(other._viewspec)) {
      // Viewspec instances are equal
      return _idx != other._idx;
    }
    auto lhs_local = lpos();
    auto rhs_local = other.lpos();
    return (lhs_local.unit  != rhs_local.unit ||
            lhs_local.index != rhs_local.index);
  }

  inline const PatternType & pattern() const
  {
    return *_pattern;
  }

  inline ddm::Team & team() const
  {
    return _pattern->team();
  }

private:
  /**
   * Compare position of this global iterator to the position of another
   * global iterator with respect to viewspec projection.
   */
  template<
    class GlobIndexCmpFun,
    class GlobPtrCmpFun >
  bool compare(
    const self_t & other,
    const GlobIndexCmpFun & gidx_cmp,
    const GlobPtrCmpFun   & gptr_cmp) const
  {
#if __REMARK__
    // Usually this is a best practice check, but it's an infrequent case
    // so we rather avoid this comparison:
    if (this == &other) {
      return true;
    }
#endif
    // NOTE:
    // Do not check _idx first, as it would never match for comparison with
    // an end iterator.
    if (_viewspec == other._viewspec) {
      // Same viewspec pointer
      return gidx_cmp(_idx, other._idx);
    }
    if ((_viewspec != nullptr && other._viewspec != nullptr) &&
        (*_viewspec) == *(other._viewspec)) {
      // Viewspec instances are equal
      return gidx_cmp(_idx, other._idx);
    }
    // View projection at lhs and/or rhs set.
    // Convert both to GlobPtr (i.e. apply view projection) and compare.
    //
    // NOTE:
    // This conversion is quite expensive but will never be necessary
    // if both iterators have been created from the same range.
    // Example:
    //   a.block(1).begin() == a.block(1).end().
    // does not require viewspace projection while
    //   a.block(1).begin() == a.end()
    // does. The latter case should be avoided for this reason.
    const pointer lhs_dart_gptr(dart_gptr());
    const pointer rhs_dart_gptr(other.dart_gptr());
    return gptr_cmp(lhs_dart_gptr, rhs_dart_gptr);
  }

  /**
   * Convert a global offset within the global iterator's range to
   * corresponding global coordinates with respect to viewspec projection.
   *
   * NOTE:
   * This method could be specialized for NumDimensions = 1 for performance
   * tuning.
   */
  std::array<IndexType, NumDimensions> coords(
    IndexType glob_index) const
  {
    DDM_LOG_TRACE_VAR("GlobViewIter.coords()", glob_index);
    // Global cartesian coords of current iterator position:
    std::array<IndexType, NumDimensions> glob_coords;
    if (_viewspec != nullptr) {
      DDM_LOG_TRACE_VAR("GlobViewIter.coords", *_viewspec);
      // Create cartesian index space from extents of view projection:
      CartesianIndexSpace<NumDimensions, Arrangement, IndexType> index_space(
        _viewspec->extents());
      // Initialize global coords with view coords (coords of iterator
      // position in view index space):
      glob_coords = index_space.coords(glob_index);
      // Apply offset of view projection to view coords:
      for (dim_t d = 0; d < NumDimensions; ++d) {
        auto dim_offset  = _viewspec->offset(d);
        glob_coords[d]  += dim_offset;
      }
    } else {
      glob_coords = _pattern->memory_layout().coords(glob_index);
    }
    DDM_LOG_TRACE_VAR("GlobViewIter.coords >", glob_coords);
    return glob_coords;
  }

}; // class GlobViewIter

/**
 * Resolve the number of elements between two global iterators.
 *
 * The difference of global pointers is not well-defined if their range
 * spans over more than one block.
 * The corresponding invariant is:
 *   g_last == g_first + (l_last - l_first)
 * Example:
 *
 * \code
 *   unit:            0       1       0
 *   local offset:  | 0 1 2 | 0 1 2 | 3 4 5 | ...
 *   global offset: | 0 1 2   3 4 5   6 7 8   ...
 *   range:          [- - -           - -]
 * \endcode
 *
 * When iterating in local memory range [0,5[ of unit 0, the position of the
 * global iterator to return is 8 != 5
 *
 * \tparam      ElementType  Type of the elements in the range
 * \complexity  O(1)
 *
 * \ingroup     Algorithms
 */
template<
  typename ElementType,
  class    Pattern,
  class    GlobMem,
  class    Pointer,
  class    Reference >
auto distance(
  const GlobViewIter<ElementType, Pattern, GlobMem, Pointer, Reference> &
    first,
  /// Global iterator to the final position in the global sequence
  const GlobViewIter<ElementType, Pattern, GlobMem, Pointer, Reference> &
    last)
-> typename Pattern::index_type
{
  return last - first;
}

template <
  typename ElementType,
  class    Pattern,
  class    GlobMem,
  class    Pointer,
  class    Reference >
std::ostream & operator<<(
  std::ostream & os,
  const ddm::GlobViewIter<
          ElementType, Pattern, GlobMem, Pointer, Reference> & it)
{
  std::ostringstream ss;
  ddm::GlobPtr<ElementType, Pattern> ptr(it);
  ss << "ddm::GlobViewIter<" << typeid(ElementType).name() << ">("
     << "idx:"  << it._idx << ", "
     << "gptr:" << ptr << ")";
  return operator<<(os, ss.str());
}

} // namespace ddm

#endif // DDM__GLOB_VIEW_ITER_H__INCLUDED
