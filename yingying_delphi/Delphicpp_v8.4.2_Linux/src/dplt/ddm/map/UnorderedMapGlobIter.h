#ifndef DDM__MAP__UNORDERED_MAP_GLOB_ITER_H__INCLUDED
#define DDM__MAP__UNORDERED_MAP_GLOB_ITER_H__INCLUDED

#include "../dart-impl/dart.h"

#include "../../ddm/Types.h"
#include "../../ddm/GlobPtr.h"
#include "../../ddm/GlobSharedRef.h"
#include "../../ddm/Allocator.h"
#include "../../ddm/Team.h"
#include "../../ddm/Onesided.h"

#include "../../ddm/map/UnorderedMapLocalIter.h"

#include "../../ddm/internal/Logging.h"

#include <type_traits>
#include <list>
#include <vector>
#include <iterator>
#include <sstream>
#include <iostream>

namespace ddm {

// Forward-declaration
template<
  typename Key,
  typename Mapped,
  typename Hash,
  typename Pred,
  typename Alloc >
class UnorderedMap;

template<
  typename Key,
  typename Mapped,
  typename Hash,
  typename Pred,
  typename Alloc >
class UnorderedMapGlobIter
: public std::iterator<
           std::random_access_iterator_tag,
           std::pair<const Key, Mapped>,
           ddm::default_index_t,
           ddm::GlobPtr< std::pair<const Key, Mapped> >,
           ddm::GlobSharedRef< std::pair<const Key, Mapped> > >
{
  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  friend class UnorderedMapGlobIter;

  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  friend std::ostream & ddm::operator<<(
    std::ostream & os,
    const ddm::UnorderedMapGlobIter<K_, M_, H_, P_, A_> & it);

private:
  typedef UnorderedMapGlobIter<Key, Mapped, Hash, Pred, Alloc>
    self_t;

  typedef UnorderedMap<Key, Mapped, Hash, Pred, Alloc>
    map_t;

  typedef UnorderedMapLocalIter<Key, Mapped, Hash, Pred, Alloc>
    local_iterator;
  typedef UnorderedMapLocalIter<Key, Mapped, Hash, Pred, Alloc>
    const_local_iterator;

public:
  typedef typename map_t::value_type                              value_type;
#if 0
  typedef typename map_t::index_type                              index_type;
  typedef typename map_t::size_type                                size_type;
#else
  typedef ddm::default_index_t                                   index_type;
  typedef ddm::default_size_t                                     size_type;
#endif

  typedef ddm::GlobPtr<      value_type>                            pointer;
  typedef ddm::GlobPtr<const value_type>                      const_pointer;
  typedef ddm::GlobSharedRef<      value_type>                    reference;
  typedef ddm::GlobSharedRef<const value_type>              const_reference;

  typedef       value_type *                                     raw_pointer;
  typedef const value_type *                               const_raw_pointer;

  typedef typename
    std::conditional<
      std::is_const<value_type>::value,
      typename map_t::const_local_iterator,
      typename map_t::local_iterator
    >::type
    local_pointer;

  typedef struct {
    team_unit_t unit;
    index_type  index;
  } local_index;

public:
  /**
   * Default constructor.
   */
  UnorderedMapGlobIter()
  : UnorderedMapGlobIter(nullptr)
  { }

  /**
   * Constructor, creates iterator at specified global position.
   */
  UnorderedMapGlobIter(
    map_t       * map,
    index_type    position)
  : _map(map),
    _idx(0),
    _myid(map->team().myid()),
    _idx_unit_id(0),
    _idx_local_idx(0)
  {
    DDM_LOG_TRACE_VAR("UnorderedMapGlobIter(map,pos)", _idx);
    increment(position);
    DDM_LOG_TRACE("UnorderedMapGlobIter(map,pos) >");
  }

  /**
   * Constructor, creates iterator at local position relative to the
   * specified unit's local iteration space.
   */
  UnorderedMapGlobIter(
    map_t         * map,
    team_unit_t     unit,
    index_type      local_index)
  : _map(map),
    _idx(0),
    _myid(map->team().myid()),
    _idx_unit_id(unit),
    _idx_local_idx(local_index)
  {
    DDM_LOG_TRACE("UnorderedMapGlobIter(map,unit,lidx)()");
    DDM_LOG_TRACE_VAR("UnorderedMapGlobIter(map,unit,lidx)", unit);
    DDM_LOG_TRACE_VAR("UnorderedMapGlobIter(map,unit,lidx)", local_index);
    // Unit and local offset to global position:
    size_type unit_l_cumul_size_prev = 0;
    if (unit > 0) {
      unit_l_cumul_size_prev = _map->_local_cumul_sizes[unit-1];
    }
    _idx = unit_l_cumul_size_prev + _idx_local_idx;
    DDM_LOG_TRACE_VAR("UnorderedMapGlobIter(map,unit,lidx)", _idx);
    DDM_LOG_TRACE("UnorderedMapGlobIter(map,unit,lidx) >");
  }

  /**
   * Copy constructor.
   */
  UnorderedMapGlobIter(
    const self_t & other) = default;

  /**
   * Assignment operator.
   */
  self_t & operator=(
    const self_t & other) = default;

  /**
   * Null-pointer constructor.
   */
  UnorderedMapGlobIter(std::nullptr_t)
  : _map(nullptr),
    _idx(-1),
    _myid(DART_UNDEFINED_UNIT_ID),
    _idx_unit_id(DART_UNDEFINED_UNIT_ID),
    _idx_local_idx(-1),
    _is_nullptr(true)
  {
    DDM_LOG_TRACE("UnorderedMapGlobIter(nullptr)");
  }

  /**
   * Null-pointer assignment operator.
   */
  self_t & operator=(std::nullptr_t) noexcept
  {
    _is_nullptr = true;
    return *this;
  }

  constexpr bool operator==(std::nullptr_t) const noexcept
  {
    return _is_nullptr;
  }

  constexpr bool operator!=(std::nullptr_t) const noexcept
  {
    return !_is_nullptr;
  }

  /**
   * Random access operator.
   */
  reference operator[](index_type offset)
  {
    auto res = *this;
    res += offset;
    return *this;
  }

  /**
   * Type conversion operator to global pointer.
   *
   * \return  A global reference to the element at the iterator's position
   */
  constexpr operator pointer() const
  {
    return pointer(dart_gptr());
  }

  /**
   * Explicit conversion to \c dart_gptr_t.
   *
   * \return  A DART global pointer to the element at the iterator's
   *          position
   */
  constexpr dart_gptr_t dart_gptr() const
  {
    return _map->globmem().at(
                            _idx_unit_id,
                            _idx_local_idx)
                          .dart_gptr();
  }

  /**
   * Dereference operator.
   *
   * \return  A global reference to the element at the iterator's position.
   */
  reference operator*()
  {
    if (is_local()) {
      // To local map iterator:
      auto l_map_it = local();
      DDM_ASSERT_MSG(l_map_it != nullptr,
                      "Converting global iterator at local position to "
                      "local iterator failed");
      // To native pointer via conversion:
      return reference(static_cast<raw_pointer>(l_map_it));
    } else {
      return reference(dart_gptr());
    }
  }

  /**
   * Dereference operator.
   *
   * \return  A global reference to the element at the iterator's position.
   */
  const_reference operator*() const
  {
    if (is_local()) {
      // To local map iterator:
      auto l_map_it = local();
      DDM_ASSERT_MSG(l_map_it != nullptr,
                      "Converting global iterator at local position to "
                      "local iterator failed");
      // To native pointer via conversion:
      return reference(static_cast<raw_pointer>(l_map_it));
    } else {
      return reference(dart_gptr());
    }
  }

#if 0
  /**
   * Requires type \c pointer to provide \c operator->().
   */
  pointer operator->() const
  {
    return static_cast<pointer>(*this);
  }
#endif

  /**
   * Checks whether the element referenced by this global iterator is in
   * the calling unit's local memory.
   */
  constexpr bool is_local() const noexcept
  {
    return (_myid == _idx_unit_id);
  }

  /**
   * Conversion to local bucket iterator.
   */
  local_iterator local()
  {
    if (_myid != _idx_unit_id) {
      // Iterator position does not point to local element
      return local_iterator(nullptr);
    }
    return (_map->lbegin() + _idx_local_idx);
  }

  /**
   * Conversion to local bucket iterator.
   */
  const_local_iterator local() const
  {
    if (_myid != _idx_unit_id) {
      // Iterator position does not point to local element
      return local_iterator(nullptr);
    }
    return (_map->lbegin() + _idx_local_idx);
  }

  /**
   * Unit and local offset at the iterator's position.
   */
  inline local_index lpos() const noexcept
  {
    local_index local_pos;
    local_pos.unit  = _idx_unit_id;
    local_pos.index = _idx_local_idx;
    return local_pos;
  }

  /**
   * Map iterator to global index domain.
   */
  constexpr self_t global() const noexcept
  {
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
   * Prefix increment operator.
   */
  inline self_t & operator++()
  {
    increment(1);
    return *this;
  }

  /**
   * Prefix decrement operator.
   */
  inline self_t & operator--()
  {
    decrement(1);
    return *this;
  }

  /**
   * Postfix increment operator.
   */
  inline self_t operator++(int)
  {
    auto result = *this;
    increment(1);
    return result;
  }

  /**
   * Postfix decrement operator.
   */
  inline self_t operator--(int)
  {
    auto result = *this;
    decrement(1);
    return result;
  }

  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  constexpr bool operator==(
    const UnorderedMapGlobIter<K_, M_, H_, P_, A_> & other) const noexcept
  {
    return (this == std::addressof(other) || _idx == other._idx);
  }

  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  constexpr bool operator!=(
    const UnorderedMapGlobIter<K_, M_, H_, P_, A_> & other) const noexcept
  {
    return !(*this == other);
  }

  self_t & operator+=(index_type offset)
  {
    increment(offset);
    return *this;
  }

  self_t & operator-=(index_type offset)
  {
    decrement(offset);
    return *this;
  }

  self_t operator+(index_type offset) const
  {
    auto res = *this;
    res += offset;
    return res;
  }

  self_t operator-(index_type offset) const
  {
    auto res = *this;
    res -= offset;
    return res;
  }

  constexpr index_type operator+(
    const self_t & other) const noexcept
  {
    return _idx + other._idx;
  }

  constexpr index_type operator-(
    const self_t & other) const noexcept
  {
    return _idx - other._idx;
  }

  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  constexpr bool operator<(
    const UnorderedMapGlobIter<K_, M_, H_, P_, A_> & other) const noexcept
  {
    return (_idx < other._idx);
  }

  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  constexpr bool operator<=(
    const UnorderedMapGlobIter<K_, M_, H_, P_, A_> & other) const noexcept
  {
    return (_idx <= other._idx);
  }

  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  constexpr bool operator>(
    const UnorderedMapGlobIter<K_, M_, H_, P_, A_> & other) const noexcept
  {
    return (_idx > other._idx);
  }

  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  constexpr bool operator>=(
    const UnorderedMapGlobIter<K_, M_, H_, P_, A_> & other) const noexcept
  {
    return (_idx >= other._idx);
  }

private:
  /**
   * Advance pointer by specified position offset.
   */
  void increment(index_type offset)
  {
    DDM_LOG_TRACE("UnorderedMapGlobIter.increment()",
                   "gidx:",   _idx, "-> (",
                   "unit:",   _idx_unit_id,
                   "lidx:",   _idx_local_idx, ")",
                   "offset:", offset);
    if (offset < 0) {
      decrement(-offset);
    } else {
      // Note:
      //
      // increment(0) is not a no-op as UnorderedMapGlobIter(map, 0) should
      // reference the first existing element, not the first possible element
      // position.
      // The first existing element has gidx:0 and lidx:0 but might not be
      // located at unit 0.
      // Example:
      //
      //     unit 0    unit 1    unit 2
      //   [ (empty) | (empty) | elem_0, elem_1 ]
      //                         |
      //                         '- first element
      //
      //   --> UnorderedMapGlobIter(map, 0) -> (gidx:0, unit:2, lidx:0)
      //
      _idx           += offset;
      _idx_local_idx = _idx;
      auto & l_cumul_sizes = _map->_local_cumul_sizes;
      // Find unit at global offset:
      while (_idx >= l_cumul_sizes[_idx_unit_id] &&
             _idx_unit_id < l_cumul_sizes.size() - 1) {
        DDM_LOG_TRACE("UnorderedMapGlobIter.increment",
                       "local cumulative size of unit", _idx_unit_id, ":",
                       l_cumul_sizes[_idx_unit_id]);
        _idx_unit_id++;
      }
      if (_idx_unit_id > 0) {
        _idx_local_idx = _idx - l_cumul_sizes[_idx_unit_id-1];
      }
    }
    DDM_LOG_TRACE("UnorderedMapGlobIter.increment >", *this);
  }

  /**
   * Decrement pointer by specified position offset.
   */
  void decrement(index_type offset)
  {
    DDM_LOG_TRACE("UnorderedMapGlobIter.decrement()",
                   "gidx:",   _idx, "-> (",
                   "unit:",   _idx_unit_id,
                   "lidx:",   _idx_local_idx, ")",
                   "offset:", -offset);
    if (offset < 0) {
      increment(-offset);
    } else if (offset > 0) {
      // TODO
      _idx           -= offset;
      _idx_local_idx  = _idx;
    }
    DDM_LOG_TRACE("UnorderedMapGlobIter.decrement >", *this);
  }

private:
  /// Pointer to referenced map instance.
  map_t                  * _map           = nullptr;
  /// Current position of the iterator in global canonical index space.
  index_type               _idx           = -1;
  /// Unit id of the active unit.
  team_unit_t              _myid;
  /// Unit id at the iterator's current position.
  team_unit_t              _idx_unit_id;
  /// Logical offset in local index space at the iterator's current position.
  index_type               _idx_local_idx = -1;
  /// Whether the iterator represents a null pointer.
  bool                     _is_nullptr    = false;

}; // class UnorderedMapGlobIter

template<
  typename Key,
  typename Mapped,
  typename Hash,
  typename Pred,
  typename Alloc >
std::ostream & operator<<(
  std::ostream & os,
  const ddm::UnorderedMapGlobIter<
          Key, Mapped, Hash, Pred, Alloc> & it)
{
  std::ostringstream ss;
  ss << "ddm::UnorderedMapGlobIter<"
     << typeid(Key).name()    << ","
     << typeid(Mapped).name() << ">"
     << "("
     << "idx:"  << it._idx           << ", "
     << "unit:" << it._idx_unit_id   << ", "
     << "lidx:" << it._idx_local_idx
     << ")";
  return operator<<(os, ss.str());
}

} // namespace ddm

#endif // DDM__MAP__UNORDERED_MAP_GLOB_ITER_H__INCLUDED
