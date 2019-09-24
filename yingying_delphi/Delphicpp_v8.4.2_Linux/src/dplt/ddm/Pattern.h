#ifndef DDM__PATTERN_H_
#define DDM__PATTERN_H_

#ifdef DOXYGEN

namespace ddm {

/**
 * \defgroup  DDMPatternConcept  Pattern Concept
 * Concept of a distribution pattern of n-dimensional containers to units in
 * a team.
 *
 * \ingroup DDMConcept
 * \{
 * \par Description
 *
 * A pattern realizes a projection of a global index range to a
 * local view:
 *
 * Team Spec                    | Container                     |
 * ---------------------------- | ----------------------------- |
 * <tt>[ unit 0 : unit 1 ]</tt> | <tt>[ 0  1  2  3  4  5 ]</tt> |
 * <tt>[ unit 1 : unit 0 ]</tt> | <tt>[ 6  7  8  9 10 11 ]</tt> |
 *
 * This pattern would assign local indices to teams like this:
 *
 * Team            | Local indices                 |
 * --------------- | ----------------------------- |
 * <tt>unit 0</tt> | <tt>[ 0  1  2  9 10 11 ]</tt> |
 * <tt>unit 1</tt> | <tt>[ 3  4  5  6  7  8 ]</tt> |
 *
 * \par Methods
 *
 * Return Type              | Method                     | Parameters                      | Description                                                                                                |
 * ------------------------ | -------------------------- | ------------------------------- | ---------------------------------------------------------------------------------------------------------- |
 * <tt>index</tt>           | <tt>local_at</tt>          | <tt>index[d] lp</tt>            | Linear local offset of the local point <i>lp</i> in local memory.                                          |
 * <tt>index</tt>           | <tt>global_at</tt>         | <tt>index[d] gp</tt>            | Global offset of the global point <i>gp</i> in the pattern's iteration order.                                 |
 * <tt>team_unit_t</tt>     | <tt>unit_at</tt>           | <tt>index[d] gp</tt>            | Unit id mapped to the element at global point <i>p</i>                                                     |
 * <b>global to local</b>   | &nbsp;                     | &nbsp;                          | &nbsp;                                                                                                     |
 * <tt>{unit,index}</tt>    | <tt>local</tt>             | <tt>index gi</tt>               | Unit and linear local offset at the global index <i>gi</i>                                                 |
 * <tt>{unit,index[d]}</tt> | <tt>local</tt>             | <tt>index[d] gp</tt>            | Unit and local coordinates at the global point <i>gp</i>                                                   |
 * <tt>{unit,index}</tt>    | <tt>local_index</tt>       | <tt>index[d] gp</tt>            | Unit and local linear offset at the global point <i>gp</i>                                                 |
 * <tt>point[d]</tt>        | <tt>local_coords</tt>      | <tt>index[d] gp</tt>            | Local coordinates at the global point <i>gp</i>                                                            |
 * <b>local to global</b>   | &nbsp;                     | &nbsp;                          | &nbsp;                                                                                                     |
 * <tt>index[d]</tt>        | <tt>global</tt>            | <tt>unit u, index[d] lp</tt>    | Local coordinates <i>lp</i> of unit <i>u</i> to global coordinates                                         |
 * <tt>index</tt>           | <tt>global_index</tt>      | <tt>unit u, index[d] lp</tt>    | Local coordinates <i>lp</i> of unit <i>u</i> to global index                                               |
 * <tt>index[d]</tt>        | <tt>global</tt>            | <tt>index[d] lp</tt>            | Local coordinates <i>lp</i> of active unit to global coordinates                                           |
 * <tt>index</tt>           | <tt>global</tt>            | <tt>unit u, index li</tt>       | Local offset <i>li</i> of unit <i>u</i> to global index                                                    |
 * <tt>index</tt>           | <tt>global</tt>            | <tt>index li</tt>               | Local offset <i>li</i> of active unit to global index                                                      |
 * <b>blocks</b>            | &nbsp;                     | &nbsp;                          | &nbsp;                                                                                                     |
 * <tt>size[d]</tt>         | <tt>blockspec</tt>         | &nbsp;                          | Number of blocks in all dimensions.                                                                        |
 * <tt>index</tt>           | <tt>block_at</tt>          | <tt>index[d] gp</tt>            | Global index of block at global coordinates <i>gp</i>                                                      |
 * <tt>viewspec</tt>        | <tt>block</tt>             | <tt>index gbi</tt>              | Offset and extent in global cartesian space of block at global block index <i>gbi</i>                      |
 * <tt>viewspec</tt>        | <tt>local_block</tt>       | <tt>index lbi</tt>              | Offset and extent in global cartesian space of block at local block index <i>lbi</i>                       |
 * <tt>viewspec</tt>        | <tt>local_block_local</tt> | <tt>index lbi</tt>              | Offset and extent in local cartesian space of block at local block index <i>lbi</i>                        |
 * <b>locality test</b>     | &nbsp;                     | &nbsp;                          | &nbsp;                                                                                                     |
 * <tt>bool</tt>            | <tt>is_local</tt>          | <tt>index gi, unit u</tt>       | Whether the global index <i>gi</i> is mapped to unit <i>u</i>                                              |
 * <tt>bool</tt>            | <tt>is_local</tt>          | <tt>dim d, index o, unit u</tt> | (proposed) Whether any element in dimension <i>d</i> at global offset <i>o</i> is local to unit <i>u</i>   |
 * <b>size</b>              | &nbsp;                     | &nbsp;                          | &nbsp;                                                                                                     |
 * <tt>size</tt>            | <tt>capacity</tt>          | &nbsp;                          | Maximum number of elements in the pattern in total                                                         |
 * <tt>size</tt>            | <tt>local_capacity</tt>    | &nbsp;                          | Maximum number of elements assigned to a single unit                                                       |
 * <tt>size</tt>            | <tt>size</tt>              | &nbsp;                          | Number of elements indexed in the pattern                                                                  |
 * <tt>size</tt>            | <tt>local_size</tt>        | &nbsp;                          | Number of elements local to the calling unit                                                               |
 * <tt>size[d]</tt>         | <tt>extents</tt>           | &nbsp;                          | Number of elements in the pattern in all dimension                                                         |
 * <tt>size</tt>            | <tt>extent</tt>            | <tt>dim d</tt>                  | Number of elements in the pattern in dimension <i>d</i>                                                    |
 * <tt>size[d]</tt>         | <tt>local_extents</tt>     | <tt>unit u</tt>                 | Number of elements local to the given unit, by dimension                                                   |
 * <tt>size</tt>            | <tt>local_extent</tt>      | <tt>dim d</tt>                  | Number of elements local to the calling unit in dimension <i>d</i>                                         |
 *
 * \}
 *
 *
 * Defines how a list of global indices is mapped to single units
 * within a Team.
 *
 * \tparam  NumDimensions  The number of dimensions of the pattern
 * \tparam  Arrangement    The memory order of the pattern (ROW_MAJOR
 *                         or COL_MAJOR), defaults to ROW_MAJOR.
 *                         Memory order defines how elements in the
 *                         pattern will be iterated predominantly
 *                         \see MemArrange
 */
template<
  dim_t      NumDimensions,
  MemArrange Arrangement   = ROW_MAJOR,
  typename   IndexType     = ddm::default_index_t
>
class PatternConcept
{
  typedef typename std::make_unsigned<IndexType>::type      SizeType;
  typedef ViewSpec<NumDimensions, IndexType>            ViewSpecType;

public:

  static constexpr char const * PatternName = "TheConcretePatternTypeName";

public:

  typedef IndexType        index_type;
  typedef SizeType          size_type;
  typedef ViewSpecType  viewspec_type;

  typedef struct {
    team_unit_t unit;
    index_type  index;
  } local_index_type;

  typedef struct {
    team_unit_t unit;
    std::array<index_type, NumDimensions> coords;
  } local_coords_type;

public:

  /**
   * Copy constructor using non-const lvalue reference parameter.
   *
   * Introduced so variadic constructor is not a better match for
   * copy-construction.
   */
  PatternConcept(self_t & other);

  /**
   * Equality comparison operator.
   */
  bool operator==(
    /// Pattern instance to compare for equality
    const self_t & other) const;

  /**
   * Inquality comparison operator.
   */
  bool operator!=(
    /// Pattern instance to compare for inequality
    const self_t & other) const;

  /**
   * Assignment operator.
   */
  PatternConcept & operator=(const PatternConcept & other);

  /**
   * Resolves the global index of the first local element in the pattern.
   *
   * \see DDMPatternConcept
   */
  index_type lbegin() const;

  /**
   * Resolves the global index past the last local element in the pattern.
   *
   * \see DDMPatternConcept
   */
  index_type lend() const;

  ////////////////////////////////////////////////////////////////////////////
  /// unit_at
  ////////////////////////////////////////////////////////////////////////////

  /**
   * Convert given point in pattern to its assigned unit id.
   *
   * \see DDMPatternConcept
   */
  team_unit_t unit_at(
    /// Absolute coordinates of the point
    const std::array<index_type, NumDimensions> & coords,
    /// View specification (offsets) to apply on \c coords
    const ViewSpec_t & viewspec) const;

  /**
   * Convert given coordinate in pattern to its assigned unit id.
   *
   * \see DDMPatternConcept
   */
  team_unit_t unit_at(
    const std::array<index_type, NumDimensions> & coords) const;

  /**
   * Convert given global linear index to its assigned unit id.
   *
   * \see DDMPatternConcept
   */
  team_unit_t unit_at(
    /// Global linear element offset
    index_type global_pos,
    /// View to apply global position
    const ViewSpec_t & viewspec) const;

  /**
   * Convert given global linear index to its assigned unit id.
   *
   * \see  blocksize()
   * \see  blockspec()
   * \see  blocksizespec()
   *
   * \see DDMPatternConcept
   */
  team_unit_t unit_at(
    /// Global linear element offset
    index_type global_pos) const;

  ////////////////////////////////////////////////////////////////////////////
  /// extent
  ////////////////////////////////////////////////////////////////////////////

  /**
   * The number of elements in this pattern in the given dimension.
   *
   * \see  blocksize()
   * \see  local_size()
   * \see  local_extent()
   *
   * \see  DDMPatternConcept
   */
  index_type extent(dim_t dim) const;

  /**
   * The actual number of elements in this pattern that are local to the
   * calling unit in the given dimension.
   *
   * \see  local_extents()
   * \see  blocksize()
   * \see  local_size()
   * \see  extent()
   *
   * \see  DDMPatternConcept
   */
  index_type local_extent(dim_t dim) const;

  /**
   * The actual number of elements in this pattern that are local to the
   * active unit, by dimension.
   *
   * \see  local_extent()
   * \see  blocksize()
   * \see  local_size()
   * \see  extent()
   *
   * \see  DDMPatternConcept
   */
  std::array<size_type, NumDimensions> local_extents() const;

  /**
   * The actual number of elements in this pattern that are local to the
   * given unit, by dimension.
   *
   * \see  local_extent()
   * \see  blocksize()
   * \see  local_size()
   * \see  extent()
   *
   * \see  DDMPatternConcept
   */
  std::array<size_type, NumDimensions> local_extents(
    team_unit_t unit) const;

  ////////////////////////////////////////////////////////////////////////////
  /// local
  ////////////////////////////////////////////////////////////////////////////

  /**
   * Convert given local coordinates and viewspec to linear local offset
   * (index).
   *
   * \see DDMPatternConcept
   */
  index_type local_at(
    /// Point in local memory
    const std::array<index_type, NumDimensions> & local_coords,
    /// View specification (offsets) to apply on \c coords
    const ViewSpec_t & viewspec) const;

  /**
   * Convert given local coordinates to linear local offset (index).
   *
   * \see DDMPatternConcept
   */
  index_type local_at(
    /// Point in local memory
    const std::array<index_type, NumDimensions> & local_coords) const;

  /**
   * Converts global coordinates to their associated unit and its respective
   * local coordinates.
   *
   * \see  DDMPatternConcept
   */
  local_coords_t local(
    const std::array<index_type, NumDimensions> & global_coords) const;

  /**
   * Converts global index to its associated unit and respective local index.
   *
   * \see  DDMPatternConcept
   */
  local_index_t local(
    index_type g_index) const;

  /**
   * Converts global coordinates to their associated unit's respective
   * local coordinates.
   *
   * \see  DDMPatternConcept
   */
  std::array<index_type, NumDimensions> local_coords(
    const std::array<index_type, NumDimensions> & global_coords) const;

  /**
   * Resolves the unit and the local index from global coordinates.
   *
   * \see  DDMPatternConcept
   */
  local_index_t local_index(
    const std::array<index_type, NumDimensions> & global_coords) const;

  ////////////////////////////////////////////////////////////////////////////
  /// global
  ////////////////////////////////////////////////////////////////////////////

  /**
   * Converts local coordinates of a given unit to global coordinates.
   *
   * \see  DDMPatternConcept
   */
  std::array<index_type, NumDimensions> global(
    team_unit_t unit,
    const std::array<index_type, NumDimensions> & local_coords) const;

  /**
   * Converts local coordinates of active unit to global coordinates.
   *
   * \see  DDMPatternConcept
   */
  std::array<index_type, NumDimensions> global(
    const std::array<index_type, NumDimensions> & local_coords) const;

  /**
   * Resolve an element's linear global index from the calling unit's local
   * index of that element.
   *
   * \see  at  Inverse of global()
   *
   * \see  DDMPatternConcept
   */
  index_type global(
    index_type local_index) const;

  /**
   * Resolve an element's linear global index from a given unit's local
   * coordinates of that element.
   *
   * \see  at
   *
   * \see  DDMPatternConcept
   */
  index_type global_index(
    team_unit_t unit,
    const std::array<index_type, NumDimensions> & local_coords) const;

  /**
   * Global coordinates and viewspec to global position in the pattern's
   * iteration order.
   *
   * \see  at
   * \see  local_at
   *
   * \see  DDMPatternConcept
   */
  index_type global_at(
    const std::array<index_type, NumDimensions> & view_coords,
    const ViewSpec_t                           & viewspec) const;

  /**
   * Global coordinates to global position in the pattern's iteration order.
   *
   * NOTE:
   * Expects extent[d] to be a multiple of blocksize[d] * nunits[d]
   * to ensure the balanced property.
   *
   * \see  at
   * \see  local_at
   *
   * \see  DDMPatternConcept
   */
  index_type global_at(
    const std::array<index_type, NumDimensions> & global_coords) const;

  ////////////////////////////////////////////////////////////////////////////
  /// at
  ////////////////////////////////////////////////////////////////////////////

  /**
   * Global coordinates to local index.
   *
   * Convert given global coordinates in pattern to their respective
   * linear local index.
   *
   * \see  DDMPatternConcept
   */
  index_type at(
    const std::array<index_type, NumDimensions> & global_coords) const;

  /**
   * Global coordinates and viewspec to local index.
   *
   * Convert given global coordinate in pattern to its linear local index.
   *
   * \see  DDMPatternConcept
   */
  index_type at(
    const std::array<index_type, NumDimensions> & global_coords,
    const ViewSpec_t & viewspec) const;

  /**
   * Global coordinates to local index.
   *
   * Convert given coordinate in pattern to its linear local index.
   *
   * \see  DDMPatternConcept
   */
  template<typename ... Values>
  index_type at(index_type value, Values ... values) const;

  ////////////////////////////////////////////////////////////////////////////
  /// is_local
  ////////////////////////////////////////////////////////////////////////////

  /**
   * Whether the given global index is local to the specified unit.
   *
   * \see  DDMPatternConcept
   */
  bool is_local(
    index_type index,
    team_unit_t unit) const;

  /**
   * Whether the given global index is local to the unit that created
   * this pattern instance.
   *
   * \see  DDMPatternConcept
   */
  bool is_local(
    index_type index) const;

  ////////////////////////////////////////////////////////////////////////////
  /// block
  ////////////////////////////////////////////////////////////////////////////

  /**
   * Cartesian arrangement of pattern blocks.
   */
  const BlockSpec_t & blockspec() const;

  /**
   * Index of block at given global coordinates.
   *
   * \see  DDMPatternConcept
   */
  index_type block_at(
    /// Global coordinates of element
    const std::array<index_type, NumDimensions> & g_coords) const;

  /**
   * View spec (offset and extents) of block at global linear block index in
   * cartesian element space.
   */
  ViewSpec_t block(
    index_type global_block_index) const;

  /**
   * View spec (offset and extents) of block at local linear block index in
   * global cartesian element space.
   */
  ViewSpec_t local_block(
    index_type local_block_index) const;

  /**
   * View spec (offset and extents) of block at local linear block index in
   * local cartesian element space.
   */
  ViewSpec_t local_block_local(
    index_type local_block_index) const;

  /**
   * Maximum number of elements in a single block in the given dimension.
   *
   * \return  The blocksize in the given dimension
   *
   * \see     DDMPatternConcept
   */
  size_type blocksize(
    /// The dimension in the pattern
    dim_t dimension) const;

  /**
   * Maximum number of elements in a single block in all dimensions.
   *
   * \return  The maximum number of elements in a single block assigned to
   *          a unit.
   *
   * \see     DDMPatternConcept
   */
  size_type max_blocksize() const;

  /**
   * Maximum number of elements assigned to a single unit in total,
   * equivalent to the local capacity of every unit in this pattern.
   *
   * \see  DDMPatternConcept
   */
  inline size_type local_capacity(
    team_unit_t unit = UNDEFINED_TEAM_UNIT_ID) const;

  /**
   * The actual number of elements in this pattern that are local to the
   * calling unit in total.
   *
   * \see  blocksize()
   * \see  local_extent()
   * \see  local_capacity()
   *
   * \see  DDMPatternConcept
   */
  inline size_type local_size(
    team_unit_t unit = UNDEFINED_TEAM_UNIT_ID) const;

  /**
   * The maximum number of elements arranged in this pattern.
   *
   * \see  DDMPatternConcept
   */
  inline index_type capacity() const;

  /**
   * The number of elements arranged in this pattern.
   *
   * \see  DDMPatternConcept
   */
  inline index_type size() const;

  /**
   * The Team containing the units to which this pattern's elements are
   * mapped.
   */
  inline ddm::Team & team() const;

  /**
   * Distribution specification of this pattern.
   */
  const DistributionSpec_t & distspec() const;

  /**
   * Size specification of the index space mapped by this pattern.
   *
   * \see DDMPatternConcept
   */
  SizeSpec_t sizespec() const;

  /**
   * Size specification of the index space mapped by this pattern.
   *
   * \see DDMPatternConcept
   */
  const std::array<size_type, NumDimensions> & extents() const;

  /**
   * Cartesian index space representing the underlying memory model of the
   * pattern.
   *
   * \see DDMPatternConcept
   */
  const MemoryLayout_t & memory_layout() const;

  /**
   * Cartesian index space representing the underlying local memory model
   * of this pattern for the calling unit.
   * Not part of DDM Pattern concept.
   */
  const LocalMemoryLayout_t & local_memory_layout() const;

  /**
   * Cartesian arrangement of the Team containing the units to which this
   * pattern's elements are mapped.
   *
   * \see DDMPatternConcept
   */
  const TeamSpec_t & teamspec() const;

  /**
   * Convert given global linear offset (index) to global cartesian
   * coordinates.
   *
   * \see DDMPatternConcept
   */
  std::array<index_type, NumDimensions> coords(
    index_type index) const;

  /**
   * Memory order followed by the pattern.
   */
  constexpr static MemArrange memory_order();

  /**
   * Number of dimensions of the cartesian space partitioned by the pattern.
   */
  constexpr static dim_t ndim();

  /**
   * Number of elements missing in the overflow block of given dimension
   * compared to the regular blocksize (\see blocksize(d)), with
   * 0 <= \c underfilled_blocksize(d) < blocksize(d).
   */
  size_type underfilled_blocksize(
    dim_t dimension) const;

};

} // namespace ddm

#endif // DOXYGEN

// Static regular pattern types:
#include "../ddm/pattern/BlockPattern.h"
#include "../ddm/pattern/TilePattern.h"
#include "../ddm/pattern/ShiftTilePattern.h"
#include "../ddm/pattern/SeqTilePattern.h"

// Static irregular pattern types:
#include "../ddm/pattern/CSRPattern.h"
#include "../ddm/pattern/LoadBalancePattern.h"

#include "../ddm/Types.h"
#include "../ddm/Distribution.h"

namespace ddm {

/**
 * Template alias for ddm::Pattern with the same default template
 * arguments
 *
 * \see BlockPattern
 */
template<
  dim_t      NumDimensions,
  MemArrange Arrangement   = ROW_MAJOR,
  typename   IndexType     = ddm::default_index_t
>
using Pattern = ddm::BlockPattern<NumDimensions, Arrangement, IndexType>;

} // namespace ddm

#include "../ddm/pattern/PatternIterator.h"
#include "../ddm/pattern/PatternProperties.h"
#include "../ddm/pattern/MakePattern.h"

#endif // DDM__PATTERN_H_
