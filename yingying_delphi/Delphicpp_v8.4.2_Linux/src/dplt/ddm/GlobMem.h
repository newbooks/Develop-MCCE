#ifndef DDM__GLOBMEM_H_
#define DDM__GLOBMEM_H_

#include "dart-impl/dart.h"

#include "../ddm/Types.h"
#include "../ddm/GlobPtr.h"
#include "../ddm/Allocator.h"
#include "../ddm/Team.h"
#include "../ddm/Onesided.h"

#include "../ddm/internal/Logging.h"

namespace ddm {

/**
 * \defgroup  DDMGlobalMemoryConcept  Global Memory Concept
 * Concept of distributed global memory space shared by units in a specified
 * team.
 *
 * \ingroup DDMConcept
 * \{
 * \par Description
 *
 * An abstraction of global memory that provides sequential iteration and
 * random access to local and global elements to units in a specified team.
 * The C++ STL does not specify a counterpart of this concept as it only
 * considers local memory that is implicitly described by the random access
 * pointer interface.
 *
 * The model of global memory represents a single, virtual global address
 * space partitioned into the local memory spaces of its associated units.
 * The global memory concept depends on the allocator concept that specifies
 * allocation of physical memory.
 *
 * Local pointers are usually, but not necessarily represented as raw native
 * pointers as returned by \c malloc.
 *
 * \see DDMAllocatorConcept
 *
 * \par Types
 *
 * Type Name            | Description                                            |
 * -------------------- | ------------------------------------------------------ |
 * \c GlobalRAI         | Random access iterator on global address space         |
 * \c LocalRAI          | Random access iterator on a single local address space |
 *
 *
 * \par Methods
 *
 * Return Type          | Method             | Parameters                         | Description                                                                                                |
 * -------------------- | ------------------ | ---------------------------------- | ---------------------------------------------------------------------------------------------------------- |
 * <tt>GlobalRAI</tt>   | <tt>begin</tt>     | &nbsp;                             | Global pointer to the initial address of the global memory space                                           |
 * <tt>GlobalRAI</tt>   | <tt>end</tt>       | &nbsp;                             | Global pointer past the final element in the global memory space                                           |
 * <tt>LocalRAI</tt>    | <tt>lbegin</tt>    | &nbsp;                             | Local pointer to the initial address in the local segment of the global memory space                       |
 * <tt>LocalRAI</tt>    | <tt>lbegin</tt>    | <tt>unit u</tt>                    | Local pointer to the initial address in the local segment at unit \c u of the global memory space          |
 * <tt>LocalRAI</tt>    | <tt>lend</tt>      | &nbsp;                             | Local pointer past the final element in the local segment of the global memory space                       |
 * <tt>LocalRAI</tt>    | <tt>lend</tt>      | <tt>unit u</tt>                    | Local pointer past the final element in the local segment at unit \c u of the global memory space          |
 * <tt>GlobalRAI</tt>   | <tt>at</tt>        | <tt>index gidx</tt>                | Global pointer to the element at canonical global offset \c gidx in the global memory space                |
 * <tt>void</tt>        | <tt>put_value</tt> | <tt>value & v_in, index gidx</tt>  | Stores value specified in parameter \c v_in to address in global memory at canonical global offset \c gidx |
 * <tt>void</tt>        | <tt>get_value</tt> | <tt>value * v_out, index gidx</tt> | Loads value from address in global memory at canonical global offset \c gidx into local address \c v_out   |
 * <tt>void</tt>        | <tt>barrier</tt>   | &nbsp;                             | Blocking synchronization of all units associated with the global memory instance                           |
 *
 * \}
 */


/**
 * Global memory with address space of static size.
 *
 * \concept{DDMGlobalMemoryConcept}
 */
template<
  /// Type of elements maintained in the global memory space
  typename ElementType,
  /// Type implementing the DDM allocator concept used to allocate and
  /// deallocate physical memory
  class    AllocatorType =
             ddm::allocator::CollectiveAllocator<ElementType> >
class GlobMem
{
private:
  typedef GlobMem<ElementType, AllocatorType>
    self_t;

public:
  typedef AllocatorType                                    allocator_type;
  typedef ElementType                                          value_type;
  typedef typename allocator_type::size_type                    size_type;
  typedef typename allocator_type::difference_type        difference_type;
  typedef typename allocator_type::difference_type             index_type;
  typedef typename allocator_type::pointer                        pointer;
  typedef typename allocator_type::void_pointer              void_pointer;
  typedef typename allocator_type::const_pointer            const_pointer;
  typedef typename allocator_type::const_void_pointer  const_void_pointer;
  typedef          value_type *                             local_pointer;
  typedef const    value_type *                       local_const_pointer;

public:
  /**
   * Constructor, collectively allocates the given number of elements in
   * local memory of every unit in a team.
   *
   * \note Must not lead to implicit barrier:
   *       Synchronization depends on underlying allocator.
   *       For example, \c ddm::LocalAllocator is used in \c ddm::Shared
   *       and only called at owner unit.
   */
  explicit GlobMem(
    /// Number of local elements to allocate in global memory space
    size_type   n_local_elem,
    /// Team containing all units operating on the global memory region
    Team      & team = ddm::Team::All())
  : _allocator(team),
    _team(team),
    _nlelem(n_local_elem),
    _nunits(team.size())
  {
    DDM_LOG_TRACE("GlobMem(nlocal,team)",
                   "number of local values:", _nlelem,
                   "team size:",              team.size());
    _begptr = _allocator.allocate(_nlelem);
    DDM_ASSERT_MSG(!DART_GPTR_ISNULL(_begptr), "allocation failed");

    // Use id's of team all
    update_lbegin();
    update_lend();
    DDM_LOG_TRACE("GlobMem(nlocal,team) >");
  }

  /**
   * Constructor, collectively allocates the given number of elements in
   * local memory of every unit in a team.
   *
   * \note Must not lead to implicit barrier:
   *       Synchronization depends on underlying allocator.
   *       For example, \c ddm::LocalAllocator is used in \c ddm::Shared
   *       and only called at owner unit.
   */
  explicit GlobMem(
    /// Local elements to allocate in global memory space
    std::initializer_list<value_type>   local_elements,
    /// Team containing all units operating on the global memory region
    Team                              & team = ddm::Team::All())
  : _allocator(team),
    _team(team),
    _nlelem(local_elements.size()),
    _nunits(team.size())
  {
    DDM_LOG_DEBUG("GlobMem(lvals,team)",
                   "number of local values:", _nlelem,
                   "team size:",              team.size());
    _begptr = _allocator.allocate(_nlelem);
    DDM_ASSERT_MSG(!DART_GPTR_ISNULL(_begptr), "allocation failed");

    // Use id's of team all
    update_lbegin();
    update_lend();
    DDM_ASSERT_EQ(std::distance(_lbegin, _lend), local_elements.size(),
                   "Capacity of local memory range differs from number "
                   "of specified local elements");

    // Initialize allocated local elements with specified values:
    auto copy_end = std::copy(local_elements.begin(),
                              local_elements.end(),
                              _lbegin);
    DDM_ASSERT_EQ(_lend, copy_end,
                   "Initialization of specified local values failed");

    if (_nunits > 1) {
      // Wait for initialization of local values at all units.
      // Barrier synchronization is okay here as multiple units are
      // involved in initialization of values in global memory:
      //
      // TODO: Should depend on allocator trait
      //         ddm::allocator_traits<Alloc>::is_collective()
      DDM_LOG_DEBUG("GlobMem(lvals,team)", "barrier");
      team.barrier();
    }

    DDM_LOG_DEBUG("GlobMem(lvals,team) >",
                   "_lbegin:", _lbegin, "_lend:", _lend);
  }

  /**
   * Destructor, collectively frees underlying global memory.
   */
  ~GlobMem()
  {
    DDM_LOG_TRACE_VAR("GlobMem.~GlobMem()", _begptr);
    _allocator.deallocate(_begptr);
    DDM_LOG_TRACE("GlobMem.~GlobMem >");
  }

  /**
   * Copy constructor.
   */
  GlobMem(const self_t & other)
    = default;

  /**
   * Assignment operator.
   */
  self_t & operator=(const self_t & rhs)
    = default;

  /**
   * Equality comparison operator.
   */
  inline bool operator==(const self_t & rhs) const
  {
    return (_begptr == rhs._begptr &&
            _team   == rhs._team   &&
            _nunits == rhs._nunits &&
            _nlelem == rhs._nlelem &&
            _lbegin == rhs._lbegin &&
            _lend   == rhs._lend);
  }

  /**
   * Inequality comparison operator.
   */
  constexpr bool operator!=(const self_t & rhs) const
  {
    return !(*this == rhs);
  }

  /**
   * Global pointer of the initial address of the global memory.
   */
  inline const GlobPtr<value_type> begin() const
  {
    return GlobPtr<value_type>(_begptr);
  }

  /**
   * Global pointer of the initial address of the global memory.
   */
  inline GlobPtr<value_type> begin()
  {
    return GlobPtr<value_type>(_begptr);
  }

  /**
   * Native pointer of the initial address of the local memory of
   * the unit that initialized this GlobMem instance.
   */
  constexpr local_const_pointer lbegin() const
  {
    return _lbegin;
  }

  /**
   * Native pointer of the initial address of the local memory of
   * the unit that initialized this GlobMem instance.
   */
  inline local_pointer lbegin()
  {
    return _lbegin;
  }


  /**
   * Native pointer of the initial address of the local memory of
   * the unit that initialized this GlobMem instance.
   */
  constexpr local_const_pointer lend() const
  {
    return _lend;
  }

  /**
   * Native pointer of the initial address of the local memory of
   * the unit that initialized this GlobMem instance.
   */
  inline local_pointer lend()
  {
    return _lend;
  }

  /**
  * Get the size of the GlobMem
  */
  inline size_type size()
  {
      return  _nlelem;
  }


  /**
   * Write value to global memory at given offset.
   *
   * \see  ddm::put_value
   */
  template<typename ValueType = value_type>
  inline void put_value(
    const ValueType & newval,
    index_type        global_index)
  {
    DDM_LOG_TRACE("GlobMem.put_value(newval, gidx = %d)", global_index);
    ddm::put_value(newval, GlobPtr<ValueType>(_begptr) + global_index);
  }

  /**
   * Read value from global memory at given offset.
   *
   * \see  ddm::get_value
   */
  template<typename ValueType = value_type>
  inline void get_value(
    ValueType  * ptr,
    index_type   global_index) const
  {
    DDM_LOG_TRACE("GlobMem.get_value(newval, gidx = %d)", global_index);
    ddm::get_value(ptr, GlobPtr<ValueType>(_begptr) + global_index);
  }

  /**
   * Synchronize all units associated with this global memory instance.
   */
  void barrier() const
  {
    _team.barrier();
  }

  /**
   * Complete all outstanding asynchronous operations on the referenced
   * global memory on all units.
   */
  void flush()
  {
    dart_flush(_begptr);
  }

  /**
   * Complete all outstanding asynchronous operations on the referenced
   * global memory on all units.
   */
  void flush_all()
  {
    dart_flush_all(_begptr);
  }

  void flush_local()
  {
    dart_flush_local(_begptr);
  }

  void flush_local_all()
  {
    dart_flush_local_all(_begptr);
  }

  /**
   * Resolve the global pointer from an element position in a unit's
   * local memory.
   */
  template<typename IndexType>
  ddm::GlobPtr<value_type> at(
    /// The unit id
    team_unit_t unit,
    /// The unit's local address offset
    IndexType   local_index) const
  {
    DDM_LOG_DEBUG("GlobMem.at(unit,l_idx)", unit, local_index);
    if (_nunits == 0 || DART_GPTR_ISNULL(_begptr)) {
      DDM_LOG_DEBUG("GlobMem.at(unit,l_idx) >",
                     "global memory not allocated");
      return ddm::GlobPtr<value_type>(nullptr);
    }
    // Initialize with global pointer to start address:
    dart_gptr_t gptr = _begptr;
    // Resolve global unit id
    DDM_LOG_TRACE_VAR("GlobMem.at (=g_begptr)", gptr);
    DDM_LOG_TRACE_VAR("GlobMem.at", gptr.unitid);
    team_unit_t lunit{gptr.unitid};
    DDM_LOG_TRACE_VAR("GlobMem.at", lunit);
    lunit = (lunit + unit) % _nunits;
    DDM_LOG_TRACE_VAR("GlobMem.at", lunit);
    // Apply global unit to global pointer:
    dart_gptr_setunit(&gptr, lunit);
    // Apply local offset to global pointer:
    ddm::GlobPtr<value_type> res_gptr(gptr);
    res_gptr += local_index;
    DDM_LOG_DEBUG("GlobMem.at (+g_unit) >", res_gptr);
    return res_gptr;
  }

private:

  /**
   * Native pointer of the initial address of the local memory of
   * a unit.
   * \param team_unit_id id of unit in \c ddm::Team::All()
   */
  void update_lbegin()
  {
    void *addr;
    dart_gptr_t gptr = _begptr;
    DDM_LOG_TRACE_VAR("GlobMem.update_lbegin",
                       GlobPtr<value_type>((dart_gptr_t)gptr));
    DDM_ASSERT_RETURNS(
      dart_gptr_setunit(&gptr, _team.myid()),
      DART_OK);
    DDM_ASSERT_RETURNS(
      dart_gptr_getaddr(gptr, &addr),
      DART_OK);
    DDM_LOG_TRACE_VAR("GlobMem.update_lbegin >", addr);
    _lbegin = static_cast<local_pointer>(addr);
  }

  /**
   * Native pointer of the final address of the local memory of
   * a unit.
   */
  void update_lend()
  {
    void *addr;
    dart_gptr_t gptr = _begptr;
    DDM_ASSERT_RETURNS(
      dart_gptr_setunit(&gptr, _team.myid()),
      DART_OK);
    DDM_ASSERT_RETURNS(
      dart_gptr_incaddr(&gptr, _nlelem * sizeof(value_type)),
      DART_OK);
    DDM_ASSERT_RETURNS(
      dart_gptr_getaddr(gptr, &addr),
      DART_OK);
    _lend = static_cast<local_pointer>(addr);
  }

private:
  allocator_type          _allocator;
  dart_gptr_t             _begptr     = DART_GPTR_NULL;
  ddm::Team&             _team;
  size_type               _nlelem     = 0;
  size_type               _nunits     = 0;
  local_pointer           _lbegin     = nullptr;
  local_pointer           _lend       = nullptr;
};

template<typename T>
GlobPtr<T> memalloc(size_t nelem)
{
  dart_gptr_t gptr;
  dart_storage_t ds = dart_storage<T>(nelem);
  if (dart_memalloc(ds.nelem, ds.dtype, &gptr) != DART_OK) {
    return GlobPtr<T>(nullptr);
  }
  return GlobPtr<T>(gptr);
}

template<typename T>
void memfree(GlobPtr<T> ptr)
{
  dart_memfree(ptr.dart_gptr());
}

} // namespace ddm

#endif // DDM__GLOBMEM_H_
