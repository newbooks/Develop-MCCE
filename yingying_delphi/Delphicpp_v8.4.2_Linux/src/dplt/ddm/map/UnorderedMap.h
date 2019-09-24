#ifndef DDM__MAP__UNORDERED_MAP_H__INCLUDED
#define DDM__MAP__UNORDERED_MAP_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/GlobRef.h"
#include "../../ddm/Team.h"
#include "../../ddm/Exception.h"
#include "../../ddm/Array.h"
#include "../../ddm/GlobDynamicMem.h"
#include "../../ddm/Allocator.h"

#include "../../ddm/atomic/GlobAtomicRef.h"

#include "../../ddm/map/UnorderedMapLocalRef.h"
#include "../../ddm/map/UnorderedMapLocalIter.h"
#include "../../ddm/map/UnorderedMapGlobIter.h"

#include <iterator>
#include <utility>
#include <limits>
#include <vector>
#include <functional>
#include <algorithm>
#include <cstddef>

namespace ddm {

template<typename Key>
class HashLocal
{
private:
  typedef ddm::default_size_t size_type;

public:
  typedef Key          argument_type;
  typedef team_unit_t result_type;

public:
  /**
   * Default constructor.
   */
  HashLocal()
  : _team(nullptr),
    _nunits(0),
    _myid(DART_UNDEFINED_UNIT_ID)
  { }

  /**
   * Constructor.
   */
  HashLocal(
    ddm::Team & team)
  : _team(&team),
    _nunits(team.size()),
    _myid(team.myid())
  { }

  result_type operator()(
    const argument_type & key) const
  {
    return _myid;
  }

private:
  ddm::Team * _team   = nullptr;
  size_type    _nunits = 0;
  team_unit_t   _myid;
}; // class HashLocal

#ifndef DOXYGEN

template<
  typename Key,
  typename Mapped,
  typename Hash    = ddm::HashLocal<Key>,
  typename Pred    = std::equal_to<Key>,
  typename Alloc   = ddm::allocator::DynamicAllocator<
                       std::pair<const Key, Mapped> > >
class UnorderedMap
{
  /**
   * The Cray compiler (as of CCE8.5.6) does not support
   * std::is_trivially_copyable.
   *
   * TODO: Remove the guard once this has been fixed by Cray.
   */
#ifndef __CRAYC
  static_assert(std::is_trivially_copyable<Key>::value,
    "Element type must be trivially copyable");
  static_assert(std::is_trivially_copyable<Mapped>::value,
    "Element type must be trivially copyable");
#endif

  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  friend class UnorderedMapLocalRef;

  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  friend class UnorderedMapGlobIter;

  template<typename K_, typename M_, typename H_, typename P_, typename A_>
  friend class UnorderedMapLocalIter;

private:
  typedef UnorderedMap<Key, Mapped, Hash, Pred, Alloc>
    self_t;

public:
  typedef Key                                                       key_type;
  typedef Mapped                                                 mapped_type;
  typedef Hash                                                        hasher;
  typedef Pred                                                     key_equal;
  typedef Alloc                                               allocator_type;

  typedef ddm::default_index_t                                   index_type;
  typedef ddm::default_index_t                              difference_type;
  typedef ddm::default_size_t                                     size_type;
  typedef std::pair<const key_type, mapped_type>                  value_type;

  typedef UnorderedMapLocalRef<Key, Mapped, Hash, Pred, Alloc>    local_type;

  typedef ddm::GlobDynamicMem<value_type, allocator_type>     glob_mem_type;

  typedef typename glob_mem_type::reference                        reference;
  typedef typename glob_mem_type::const_reference            const_reference;

  typedef typename reference::template rebind<mapped_type>::other
    mapped_type_reference;
  typedef typename const_reference::template rebind<mapped_type>::other
    const_mapped_type_reference;

  typedef typename glob_mem_type::global_iterator
    node_iterator;
  typedef typename glob_mem_type::const_global_iterator
    const_node_iterator;
  typedef typename glob_mem_type::local_iterator
    local_node_iterator;
  typedef typename glob_mem_type::const_local_iterator
    const_local_node_iterator;
  typedef typename glob_mem_type::reverse_global_iterator
    reverse_node_iterator;
  typedef typename glob_mem_type::const_reverse_global_iterator
    const_reverse_node_iterator;
  typedef typename glob_mem_type::reverse_local_iterator
    reverse_local_node_iterator;
  typedef typename glob_mem_type::const_reverse_local_iterator
    const_reverse_local_node_iterator;

  typedef typename glob_mem_type::global_iterator
    local_node_pointer;
  typedef typename glob_mem_type::const_global_iterator
    const_local_node_pointer;

  typedef UnorderedMapGlobIter<Key, Mapped, Hash, Pred, Alloc>
    iterator;
  typedef UnorderedMapGlobIter<Key, Mapped, Hash, Pred, Alloc>
    const_iterator;
  typedef typename std::reverse_iterator<iterator>
    reverse_iterator;
  typedef typename std::reverse_iterator<const_iterator>
    const_reverse_iterator;

  typedef UnorderedMapLocalIter<Key, Mapped, Hash, Pred, Alloc>
    local_pointer;
  typedef UnorderedMapLocalIter<Key, Mapped, Hash, Pred, Alloc>
    const_local_pointer;
  typedef UnorderedMapLocalIter<Key, Mapped, Hash, Pred, Alloc>
    local_iterator;
  typedef UnorderedMapLocalIter<Key, Mapped, Hash, Pred, Alloc>
    const_local_iterator;
  typedef typename std::reverse_iterator<local_iterator>
    reverse_local_iterator;
  typedef typename std::reverse_iterator<const_local_iterator>
    const_reverse_local_iterator;

  typedef typename glob_mem_type::local_reference
    local_reference;
  typedef typename glob_mem_type::const_local_reference
    const_local_reference;

  typedef ddm::Array<
            size_type, int, ddm::CSRPattern<1, ddm::ROW_MAJOR, int> >
    local_sizes_map;

private:
  /// Team containing all units interacting with the map.
  ddm::Team           * _team            = nullptr;
  /// DART id of the local unit.
  team_unit_t            _myid{DART_UNDEFINED_UNIT_ID};
  /// Global memory allocation and -access.
  glob_mem_type        * _globmem         = nullptr;
  /// Iterator to initial element in the map.
  iterator               _begin           = nullptr;
  /// Iterator past the last element in the map.
  iterator               _end             = nullptr;
  /// Number of elements in the map.
  size_type              _remote_size     = 0;
  /// Native pointer to first local element in the map.
  local_iterator         _lbegin          = nullptr;
  /// Native pointer past the last local element in the map.
  local_iterator         _lend            = nullptr;
  /// Mapping units to their number of local map elements.
  local_sizes_map        _local_sizes;
  /// Cumulative (postfix sum) local sizes of all units.
  std::vector<size_type> _local_cumul_sizes;
  /// Iterators to elements in local memory space that are marked for move
  /// to remote unit in next commit.
  std::vector<iterator>  _move_elements;
  /// Global pointer to local element in _local_sizes.
  dart_gptr_t            _local_size_gptr = DART_GPTR_NULL;
  /// Hash type for mapping of key to unit and local offset.
  hasher                 _key_hash;
  /// Predicate for key comparison.
  key_equal              _key_equal;
  /// Capacity of local buffer containing locally added node elements that
  /// have not been committed to global memory yet.
  /// Default is 4 KB.
  size_type              _local_buffer_size
                           = 4096 / sizeof(value_type);

public:
  /// Local proxy object, allows use in range-based for loops.
  local_type local;

public:
  UnorderedMap(
    size_type   nelem = 0,
    Team      & team  = ddm::Team::All())
  : _team(&team),
    _myid(team.myid()),
    _key_hash(team),
    local(this)
  {
    DDM_LOG_TRACE_VAR("UnorderedMap(nelem,team)", nelem);
    if (_team->size() > 0) {
      allocate(nelem, team);
    }
    DDM_LOG_TRACE("UnorderedMap(nelem,team) >");
  }

  UnorderedMap(
    size_type   nelem,
    size_type   nlbuf,
    Team      & team  = ddm::Team::All())
  : _team(&team),
    _myid(team.myid()),
    _key_hash(team),
    _local_buffer_size(nlbuf),
    local(this)
  {
    DDM_LOG_TRACE("UnorderedMap(nelem,nlbuf,team)",
                   "nelem:", nelem, "nlbuf:", nlbuf);
    if (_team->size() > 0) {
      allocate(nelem, team);
    }
    DDM_LOG_TRACE("UnorderedMap(nelem,nlbuf,team) >");
  }

  ~UnorderedMap()
  {
    DDM_LOG_TRACE_VAR("UnorderedMap.~UnorderedMap()", this);
    deallocate();
    DDM_LOG_TRACE_VAR("UnorderedMap.~UnorderedMap >", this);
  }

  //////////////////////////////////////////////////////////////////////////
  // Distributed container
  //////////////////////////////////////////////////////////////////////////

  inline const Team & team() const noexcept
  {
    if (_team != nullptr) {
      return *_team;
    }
    return ddm::Team::Null();
  }

  inline const glob_mem_type & globmem() const
  {
    return *_globmem;
  }

  //////////////////////////////////////////////////////////////////////////
  // Dynamic distributed memory
  //////////////////////////////////////////////////////////////////////////

  void barrier()
  {
    DDM_LOG_TRACE_VAR("UnorderedMap.barrier()", _team->dart_id());
    // Apply changes in local memory spaces to global memory space:
    if (_globmem != nullptr) {
      _globmem->commit();
    }
    // Accumulate local sizes of remote units:
    _local_sizes.barrier();
    _remote_size = 0;
    for (int u = 0; u < _team->size(); ++u) {
      size_type local_size_u;
      if (u != _myid) {
        local_size_u = _local_sizes[u];
        _remote_size += local_size_u;
      } else {
        local_size_u = _local_sizes.local[0];
      }
      _local_cumul_sizes[u] = local_size_u;
      if (u > 0) {
        _local_cumul_sizes[u] += _local_cumul_sizes[u-1];
      }
      DDM_LOG_TRACE("UnorderedMap.barrier",
                     "local size at unit", u, ":", local_size_u,
                     "cumulative size:", _local_cumul_sizes[u]);
    }
    auto new_size = size();
    DDM_LOG_TRACE("UnorderedMap.barrier", "new size:", new_size);
    DDM_ASSERT_EQ(_remote_size, new_size - _local_sizes.local[0],
                   "invalid size after global commit");
    _begin = iterator(this, 0);
    _end   = iterator(this, new_size);
    DDM_LOG_TRACE("UnorderedMap.barrier >", "passed barrier");
  }

  bool allocate(
    /// Initial global capacity of the container.
    size_type    nelem = 0,
    /// Team containing all units associated with the container.
    ddm::Team & team  = ddm::Team::All())
  {
    DDM_LOG_TRACE("UnorderedMap.allocate()");
    DDM_LOG_TRACE_VAR("UnorderedMap.allocate", nelem);
    DDM_LOG_TRACE_VAR("UnorderedMap.allocate", _local_buffer_size);
    if (_team == nullptr || *_team == ddm::Team::Null()) {
      DDM_LOG_TRACE("UnorderedMap.allocate",
                     "initializing with specified team -",
                     "team size:", team.size());
      _team = &team;
      DDM_LOG_TRACE_VAR("UnorderedMap.allocate", team.dart_id());
    } else {
      DDM_LOG_TRACE("UnorderedMap.allocate",
                     "initializing with initial team");
    }
    _local_cumul_sizes = std::vector<size_type>(_team->size(), 0);
    DDM_ASSERT_GT(_local_buffer_size, 0, "local buffer size must not be 0");
    if (nelem < _team->size() * _local_buffer_size) {
      nelem = _team->size() * _local_buffer_size;
      DDM_LOG_TRACE("UnorderedMap.allocate", "nelem increased to", nelem);
    }
    _key_hash    = hasher(*_team);
    _remote_size = 0;
    auto lcap    = ddm::math::div_ceil(nelem, _team->size());
    // Initialize members:
    _myid        = _team->myid();

    DDM_LOG_TRACE("UnorderedMap.allocate", "initialize global memory,",
                   "local capacity:", lcap);
    _globmem     = new glob_mem_type(lcap, *_team);
    DDM_LOG_TRACE("UnorderedMap.allocate", "global memory initialized");

    // Initialize local sizes with 0:
    _local_sizes.allocate(_team->size(), ddm::BLOCKED, *_team);
    _local_sizes.local[0] = 0;
    _local_size_gptr      = _local_sizes[_myid].dart_gptr();

    // Global iterators:
    _begin       = iterator(this, 0);
    _end         = _begin;
    DDM_LOG_TRACE_VAR("UnorderedMap.allocate", _begin);
    DDM_LOG_TRACE_VAR("UnorderedMap.allocate", _end);
    // Local iterators:
    _lbegin      = local_iterator(this, 0);
    _lend        = _lbegin;
    DDM_LOG_TRACE_VAR("UnorderedMap.allocate", _lbegin);
    DDM_LOG_TRACE_VAR("UnorderedMap.allocate", _lend);
    // Register deallocator of this map instance at the team
    // instance that has been used to initialized it:
    _team->register_deallocator(
             this, std::bind(&UnorderedMap::deallocate, this));
    // Assure all units are synchronized after allocation, otherwise
    // other units might start working on the map before allocation
    // completed at all units:
    if (ddm::is_initialized()) {
      DDM_LOG_TRACE("UnorderedMap.allocate",
                     "waiting for allocation of all units");
      _team->barrier();
    }
    DDM_LOG_TRACE("UnorderedMap.allocate >", "finished");
    return true;
  }

  void deallocate()
  {
    DDM_LOG_TRACE_VAR("UnorderedMap.deallocate()", this);
    // Assure all units are synchronized before deallocation, otherwise
    // other units might still be working on the map:
    if (ddm::is_initialized()) {
      barrier();
    }
    // Remove this function from team deallocator map to avoid
    // double-free:
    _team->unregister_deallocator(
      this, std::bind(&UnorderedMap::deallocate, this));
    // Deallocate map elements:
    DDM_LOG_TRACE_VAR("UnorderedMap.deallocate()", _globmem);
    if (_globmem != nullptr) {
      delete _globmem;
      _globmem = nullptr;
    }
    _local_cumul_sizes    = std::vector<size_type>(_team->size(), 0);
    _local_sizes.local[0] = 0;
    _remote_size          = 0;
    _begin                = iterator();
    _end                  = _begin;
    DDM_LOG_TRACE_VAR("UnorderedMap.deallocate >", this);
  }

  //////////////////////////////////////////////////////////////////////////
  // Global Iterators
  //////////////////////////////////////////////////////////////////////////

  inline iterator & begin() noexcept
  {
    return _begin;
  }

  inline const_iterator & begin() const noexcept
  {
    return _begin;
  }

  inline const_iterator & cbegin() const noexcept
  {
    return const_cast<const self_t *>(this)->begin();
  }

  inline iterator & end() noexcept
  {
    return _end;
  }

  inline const_iterator & end() const noexcept
  {
    return _end;
  }

  inline const_iterator & cend() const noexcept
  {
    return const_cast<const self_t *>(this)->end();
  }

  //////////////////////////////////////////////////////////////////////////
  // Local Iterators
  //////////////////////////////////////////////////////////////////////////

  inline local_iterator & lbegin() noexcept
  {
    return _lbegin;
  }

  inline const_local_iterator & lbegin() const noexcept
  {
    return _lbegin;
  }

  inline const_local_iterator & clbegin() const noexcept
  {
    return const_cast<const self_t *>(this)->lbegin();
  }

  inline local_iterator & lend() noexcept
  {
    return _lend;
  }

  inline const_local_iterator & lend() const noexcept
  {
    return _lend;
  }

  inline const_local_iterator & clend() const noexcept
  {
    return const_cast<const self_t *>(this)->lend();
  }

  //////////////////////////////////////////////////////////////////////////
  // Capacity
  //////////////////////////////////////////////////////////////////////////

  constexpr size_type max_size() const noexcept
  {
    return std::numeric_limits<key_type>::max();
  }

  inline size_type size() const noexcept
  {
    return _remote_size + _local_sizes.local[0];
  }

  inline size_type capacity() const noexcept
  {
    return _globmem->size();
  }

  inline bool empty() const noexcept
  {
    return size() == 0;
  }

  inline size_type lsize() const noexcept
  {
    return _local_sizes.local[0];
  }

  inline size_type lcapacity() const noexcept
  {
    return _globmem != nullptr
           ? _globmem->local_size()
           : 0;
  }

  //////////////////////////////////////////////////////////////////////////
  // Element Access
  //////////////////////////////////////////////////////////////////////////

  mapped_type_reference operator[](const key_type & key)
  {
    DDM_LOG_TRACE("UnorderedMap.[]()", "key:", key);
    iterator      git_value   = insert(
                                   std::make_pair(key, mapped_type()))
                                .first;
    DDM_LOG_TRACE_VAR("UnorderedMap.[]", git_value);
    dart_gptr_t   gptr_mapped = git_value.dart_gptr();
    value_type  * lptr_value  = static_cast<value_type *>(
                                  git_value.local());
    mapped_type * lptr_mapped = nullptr;

    _lptr_value_to_mapped(lptr_value, gptr_mapped, lptr_mapped);
    // Create global reference to mapped value member in element:
    mapped_type_reference mapped(gptr_mapped,
                                 lptr_mapped);
    DDM_LOG_TRACE("UnorderedMap.[] >", mapped);
    return mapped;
  }

  const_mapped_type_reference at(const key_type & key) const
  {
    DDM_LOG_TRACE("UnorderedMap.at() const", "key:", key);
    auto found = find(key);
    if (found == _end) {
      // No equivalent key in map, throw:
      DDM_THROW(
        ddm::exception::InvalidArgument,
        "No element in map for key " << key);
    }
    dart_gptr_t gptr_mapped   = iterator(this, found.pos()).dart_gptr();

    value_type  * lptr_value  = static_cast<value_type *>(
                                  found.local());
    mapped_type * lptr_mapped = nullptr;

    _lptr_value_to_mapped(lptr_value, gptr_mapped, lptr_mapped);
    // Create global reference to mapped value member in element:
    const_mapped_type_reference mapped(gptr_mapped,
                                       lptr_mapped);
    return mapped;
  }

  mapped_type_reference at(const key_type & key)
  {
    DDM_LOG_TRACE("UnorderedMap.at()", "key:", key);
    // TODO: Unoptimized, currently calls find(key) twice as operator[](key)
    //       calls insert(key).
    const_iterator git_value = find(key);
    if (git_value == _end) {
      // No equivalent key in map, throw:
      DDM_THROW(
        ddm::exception::InvalidArgument,
        "No element in map for key " << key);
    }
    auto mapped = this->operator[](key);
    DDM_LOG_TRACE("UnorderedMap.at >", mapped);
    return mapped;
  }

  //////////////////////////////////////////////////////////////////////////
  // Element Lookup
  //////////////////////////////////////////////////////////////////////////

  size_type count(const key_type & key) const
  {
    DDM_LOG_TRACE_VAR("UnorderedMap.count()", key);
    size_type nelem = 0;
    if (find(key) != _end) {
      nelem = 1;
    }
    DDM_LOG_TRACE("UnorderedMap.count >", nelem);
    return nelem;
  }

  iterator find(const key_type & key)
  {
    DDM_LOG_TRACE_VAR("UnorderedMap.find()", key);
    iterator found = std::find_if(
                       _begin, _end,
                       [&](const value_type & v) {
                         return _key_equal(v.first, key);
                       });
    DDM_LOG_TRACE("UnorderedMap.find >", found);
    return found;
  }

  const_iterator find(const key_type & key) const
  {
    DDM_LOG_TRACE_VAR("UnorderedMap.find() const", key);
    const_iterator found = std::find_if(
                             _begin, _end,
                             [&](const value_type & v) {
                               return _key_equal(v.first, key);
                             });
    DDM_LOG_TRACE("UnorderedMap.find const >", found);
    return found;
  }

  //////////////////////////////////////////////////////////////////////////
  // Modifiers
  //////////////////////////////////////////////////////////////////////////

  std::pair<iterator, bool> insert(
    /// The element to insert.
    const value_type & value)
  {
    auto && key = value.first;
    DDM_LOG_TRACE("UnorderedMap.insert()", "key:", key);

    auto result = std::make_pair(_end, false);

    DDM_ASSERT(_globmem != nullptr);
    // Look up existing element at given key:
    DDM_LOG_TRACE("UnorderedMap.insert", "element key lookup");
    const_iterator found = find(key);
    DDM_LOG_TRACE_VAR("UnorderedMap.insert", found);

    if (found != _end) {
      DDM_LOG_TRACE("UnorderedMap.insert", "key found");
      // Existing element found, no insertion:
      result.first  = iterator(this, found.pos());
      result.second = false;
    } else {
      DDM_LOG_TRACE("UnorderedMap.insert", "key not found");
      // Unit mapped to the new element's key by the hash function:
      auto unit = _key_hash(key);
      DDM_LOG_TRACE("UnorderedMap.insert", "target unit:", unit);
      // No element with specified key exists, insert new value.
      result = _insert_at(unit, value);
    }
    DDM_LOG_DEBUG("UnorderedMap.insert >",
                   (result.second ? "inserted" : "existing"), ":",
                   result.first);
    return result;
  }

  template<class InputIterator>
  void insert(
    // Iterator at first value in the range to insert.
    InputIterator first,
    // Iterator past the last value in the range to insert.
    InputIterator last)
  {
    // TODO: Calling insert() on every single element in the range could cause
    //       multiple calls of globmem.grow(_local_buffer_size).
    //       Could be optimized to allocate additional memory in a single call
    //       of globmem.grow(std::distance(first,last)).
    for (auto it = first; it != last; ++it) {
      insert(*it);
    }
  }

  iterator erase(
    const_iterator position)
  {
    DDM_THROW(
      ddm::exception::NotImplemented,
      "ddm::UnorderedMap.erase is not implemented.");
  }

  size_type erase(
    /// Key of the container element to remove.
    const key_type & key)
  {
    DDM_THROW(
      ddm::exception::NotImplemented,
      "ddm::UnorderedMap.erase is not implemented.");
  }

  iterator erase(
    /// Iterator at first element to remove.
    const_iterator first,
    /// Iterator past the last element to remove.
    const_iterator last)
  {
    DDM_THROW(
      ddm::exception::NotImplemented,
      "ddm::UnorderedMap.erase is not implemented.");
  }

  //////////////////////////////////////////////////////////////////////////
  // Bucket Interface
  //////////////////////////////////////////////////////////////////////////

  inline size_type bucket(const key_type & key) const
  {
    return _key_hash(key);
  }

  inline size_type bucket_size(size_type bucket_index) const
  {
    size_type   bsize = _local_sizes[bucket_index];
    return bsize;
  }

  //////////////////////////////////////////////////////////////////////////
  // Observers
  //////////////////////////////////////////////////////////////////////////

  inline key_equal key_eq() const
  {
    return _key_equal;
  }

  inline hasher hash_function() const
  {
    return _key_hash;
  }

private:
  /**
   * Helper to resolve address of mapped value from map entries.
   *
   * std::pair cannot be used as MPI data type directly.
   * Offset-to-member only works reliably with offsetof in the general case
   * We have to use `offsetof` as there is no instance of value_type
   * available that could be used to calculate the member offset as
   * `l_ptr_value` is possibly undefined.
   *
   * Using `std::declval()` instead (to generate a compile-time
   * pseudo-instance for member resolution) only works if Key and Mapped
   * are default-constructible.
   * 
   * Finally, the distance obtained from
   *
   *   &(lptr_value->second) - lptr_value
   *
   * had different alignment than the address obtained via offsetof in some
   * cases, depending on the combination of MPI runtime and compiler.
   * Apparently some compilers / standard libs have special treatment
   * (padding?) for members of std::pair such that
   *
   *   __builtin_offsetof(type, member)
   *
   * differs from the member-offset provided by the type system.
   * The alternative, using `offsetof` (resolves to `__builtin_offsetof`
   * automatically if needed) and manual pointer increments works, however.
   */
  void _lptr_value_to_mapped(
    // [IN]  native pointer to map entry
    value_type    * lptr_value,
    // [OUT] corresponding global pointer to mapped value
    dart_gptr_t   & gptr_mapped,
    // [OUT] corresponding native pointer to mapped value
    mapped_type * & lptr_mapped) const
  {
    // Byte offset of mapped value in element type:
    auto mapped_offs = offsetof(value_type, second);
    DDM_LOG_TRACE("UnorderedMap.lptr_value_to_mapped()",
                   "byte offset of mapped member:", mapped_offs);
    // Increment pointers to element by byte offset of mapped value member:
    if (lptr_value != nullptr) {
        if (std::is_standard_layout<value_type>::value) {
        // Convert to char pointer for byte-wise increment:
        char * b_lptr_mapped = reinterpret_cast<char *>(lptr_value);
        b_lptr_mapped       += mapped_offs;
        // Convert to mapped type pointer:
        lptr_mapped          = reinterpret_cast<mapped_type *>(b_lptr_mapped);
      } else {
        lptr_mapped = &(lptr_value->second);
      }
    }
    if (!DART_GPTR_ISNULL(gptr_mapped)) {
      DDM_ASSERT_RETURNS(
        dart_gptr_incaddr(&gptr_mapped, mapped_offs),
        DART_OK);
    }
    DDM_LOG_TRACE("UnorderedMap.lptr_value_to_mapped >",
                   "gptr to mapped:", gptr_mapped);
    DDM_LOG_TRACE("UnorderedMap.lptr_value_to_mapped >",
                   "lptr to mapped:", lptr_mapped);
  }

  /**
   * Insert value at specified unit.
   */
  std::pair<iterator, bool> _insert_at(
    team_unit_t        unit,
    /// The element to insert.
    const value_type & value)
  {
    DDM_LOG_TRACE("UnorderedMap._insert_at()",
                   "unit:",   unit,
                   "key:",    value.first);
    auto result = std::make_pair(_end, false);
    // Increase local size first to reserve storage for the new element.
    // Use atomic increment to prevent hazard when other units perform
    // remote insertion at the local unit:
    size_type old_local_size   = GlobRef<Atomic<size_type>>(
                                    _local_size_gptr
                                 ).fetch_add(1);
    size_type new_local_size   = old_local_size + 1;
    size_type local_capacity   = _globmem->local_size();
    _local_cumul_sizes[unit]  += 1;
    DDM_LOG_TRACE_VAR("UnorderedMap._insert_at", local_capacity);
    DDM_LOG_TRACE_VAR("UnorderedMap._insert_at", _local_buffer_size);
    DDM_LOG_TRACE_VAR("UnorderedMap._insert_at", old_local_size);
    DDM_LOG_TRACE_VAR("UnorderedMap._insert_at", new_local_size);
    DDM_LOG_TRACE_VAR("UnorderedMap._insert_at", _local_cumul_sizes[_myid]);
    DDM_ASSERT_GT(new_local_size, 0, "new local size is 0");
    // Pointer to new element:
    value_type * lptr_insert = nullptr;
    // Acquire target pointer of new element:
    if (new_local_size > local_capacity) {
      DDM_LOG_TRACE("UnorderedMap._insert_at",
                     "globmem.grow(", _local_buffer_size, ")");
      lptr_insert = static_cast<value_type *>(
                      _globmem->grow(_local_buffer_size));
    } else {
      lptr_insert = static_cast<value_type *>(
                      _globmem->lbegin() + old_local_size);
    }
    // Assign new value to insert position.
    DDM_LOG_TRACE("UnorderedMap._insert_at", "value target address:",
                   lptr_insert);
    DDM_ASSERT(lptr_insert != nullptr);
    // Using placement new to avoid assignment/copy as value_type is
    // const:
    new (lptr_insert) value_type(value);
    // Convert local iterator to global iterator:
    DDM_LOG_TRACE("UnorderedMap._insert_at", "converting to global iterator",
                   "unit:", unit, "lidx:", old_local_size);
    result.first  = iterator(this, unit, old_local_size);
    result.second = true;

    if (unit != _myid) {
      DDM_LOG_TRACE("UnorderedMap.insert", "remote insertion");
      // Mark inserted element for move to remote unit in next commit:
      _move_elements.push_back(result.first);
    }

    // Update iterators as global memory space has been changed for the
    // active unit:
    auto new_size = size();
    DDM_LOG_TRACE("UnorderedMap._insert_at", "new size:", new_size);
    DDM_LOG_TRACE("UnorderedMap._insert_at", "updating _begin");
    _begin        = iterator(this, 0);
    DDM_LOG_TRACE("UnorderedMap._insert_at", "updating _end");
    _end          = iterator(this, new_size);
    DDM_LOG_TRACE_VAR("UnorderedMap._insert_at", _begin);
    DDM_LOG_TRACE_VAR("UnorderedMap._insert_at", _end);
    DDM_LOG_DEBUG("UnorderedMap._insert_at >",
                   (result.second ? "inserted" : "existing"), ":",
                   result.first);
    return result;
  }

}; // class UnorderedMap

#endif // ifndef DOXYGEN

} // namespace ddm

#endif // DDM__MAP__UNORDERED_MAP_H__INCLUDED
