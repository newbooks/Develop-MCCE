#ifndef DDM__MATRIX_H_INCLUDED
#define DDM__MATRIX_H_INCLUDED

#include "dart-impl/dart.h"

#include "../ddm/Team.h"
#include "../ddm/Pattern.h"
#include "../ddm/GlobRef.h"
#include "../ddm/GlobMem.h"
#include "../ddm/Allocator.h"
#include "../ddm/HView.h"

#include "../ddm/iterator/GlobIter.h"

#include "../ddm/matrix/MatrixRefView.h"
#include "../ddm/matrix/MatrixRef.h"
#include "../ddm/matrix/LocalMatrixRef.h"

#include <type_traits>


/**
 * \defgroup  DDMMatrixConcept  Matrix Concept
 * Concept for a distributed n-dimensional matrix.
 *
 * Extends concepts \c DDMContainerConcept and \c DDMArrayConcept.
 *
 * \see DDMContainerConcept
 * \see DDMArrayConcept
 * \see DDMViewConcept
 *
 * \ingroup DDMContainerConcept
 * \{
 * \par Description
 *
 * The Matrix concept extends the n-dimensional Array by
 * operations that are prevalent in linear algebra, such
 * as projection to lower dimensions (e.g. rows and columns)
 * or slices.
 *
 * \par Types
 *
 * As defined in \c DDMContainerConcept.
 *
 * Type name                       | Description
 * ------------------------------- | -------------------------------------------------------------------------------------------------------------------
 * <tt>value_type</tt>             | Type of the container elements.
 * <tt>difference_type</tt>        | Integer type denoting a distance in cartesian index space.
 * <tt>index_type</tt>             | Integer type denoting an offset/coordinate in cartesian index space.
 * <tt>size_type</tt>              | Integer type denoting an extent in cartesian index space.
 * <tt>iterator</tt>               | Iterator on container elements in global index space.
 * <tt>const_iterator</tt>         | Iterator on const container elements in global index space.
 * <tt>reverse_iterator</tt>       | Reverse iterator on container elements in global index space.
 * <tt>const_reverse_iterator</tt> | Reverse iterator on const container elements in global index space.
 * <tt>reference</tt>              | Reference on container elements in global index space.
 * <tt>const_reference</tt>        | Reference on const container elements in global index space.
 * <tt>local_pointer</tt>          | Native pointer on local container elements.
 * <tt>const_local_pointer</tt>    | Native pointer on const local container elements.
 * <tt>view_type</tt>              | View specifier on container elements, model of \c DDMViewConcept.
 * <tt>local_type</tt>             | Reference to local element range, allows range-based iteration.
 * <tt>pattern_type</tt>           | Concrete model of the Pattern concept that specifies the container's data distribution and cartesian access pattern.
 *
 *
 * \par Methods
 *
 * Return Type              | Method                | Parameters                      | Description
 * ------------------------ | --------------------- | ------------------------------- | -------------------------------------------------------------------------------------------------------------
 * <tt>view_type<d></tt>    | <tt>block</tt>        | <tt>index_type    bi</tt>       | Matrix proxy object representing a view specifier on the matrix block at canonical block index <tt>bi</tt>.
 * <tt>view_type<d></tt>    | <tt>block</tt>        | <tt>index_type[d] bp</tt>       | Matrix proxy object representing a view specifier on the matrix block at block coordinate <tt>bc</tt>.
 *
 * As defined in \c DDMContainerConcept:
 *
 * Return Type              | Method                | Parameters                                            | Description
 * ------------------------ | --------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------
 * <tt>local_type</tt>      | <tt>local</tt>        | &nbsp;                                                | Container proxy object representing a view specifier on the container's local elements.
 * <tt>pattern_type</tt>    | <tt>pattern</tt>      | &nbsp;                                                | Object implementing the Pattern concept specifying the container's data distribution and iteration pattern.
 * <tt>iterator</tt>        | <tt>begin</tt>        | &nbsp;                                                | Iterator referencing the first container element.
 * <tt>iterator</tt>        | <tt>end</tt>          | &nbsp;                                                | Iterator referencing the element past the last container element.
 * <tt>Element *</tt>       | <tt>lbegin</tt>       | &nbsp;                                                | Native pointer referencing the first local container element, same as <tt>local().begin()</tt>.
 * <tt>Element *</tt>       | <tt>lend</tt>         | &nbsp;                                                | Native pointer referencing the element past the last local container element, same as <tt>local().end()</tt>.
 * <tt>size_type</tt>       | <tt>size</tt>         | &nbsp;                                                | Number of elements in the container.
 * <tt>size_type</tt>       | <tt>local_size</tt>   | &nbsp;                                                | Number of local elements in the container, same as <tt>local().size()</tt>.
 * <tt>bool</tt>            | <tt>is_local</tt>     | <tt>index_type gi</tt>                                | Whether the element at the given linear offset in global index space <tt>gi</tt> is local.
 * <tt>bool</tt>            | <tt>allocate</tt>     | <tt>size_type n, DistributionSpec<DD> ds, Team t</tt> | Allocation of <tt>n</tt> container elements distributed in Team <tt>t</tt> as specified by distribution spec <tt>ds</tt>
 * <tt>void</tt>            | <tt>deallocate</tt>   | &nbsp;                                                | Deallocation of the container and its elements.
 *
 * \}
 */

namespace ddm {

/// Forward-declaration
template <
  typename T,
  dim_t NumDimensions,
  typename IndexT,
  class PatternT >
class Matrix;
/// Forward-declaration
template <
  typename T,
  dim_t NumDimensions,
  dim_t CUR,
  class PatternT >
class MatrixRef;
/// Forward-declaration
template <
  typename T,
  dim_t NumDimensions,
  dim_t CUR,
  class PatternT >
class LocalMatrixRef;

/**
 * An n-dimensional array supporting subranges and sub-dimensional
 * projection.
 *
 * \concept{DDMMatrixConcept}
 *
 * \todo
 * Projection order matrix.sub().local() is not fully implemented yet.
 * Currently only matrix.local().sub() is supported.
 *
 * \note
 * Roughly follows the design presented in
 *   "The C++ Programming Language" (Bjarne Stroustrup)
 *   Chapter 29: A Matrix Design
 *
 */
template<
  typename ElementT,
  dim_t    NumDimensions,
  typename IndexT         = ddm::default_index_t,
  class    PatternT       = TilePattern<NumDimensions, ROW_MAJOR, IndexT> >
class Matrix
{
  /**
   * The Cray compiler (as of CCE8.5.6) does not support
   * std::is_trivially_copyable.
   *
   * TODO: Remove the guard once this has been fixed by Cray.
   */
#ifndef __CRAYC
  static_assert(std::is_trivially_copyable<ElementT>::value,
    "Element type must be trivially copyable");
#endif
  static_assert(std::is_same<IndexT, typename PatternT::index_type>::value,
    "Index type IndexT must be the same for Matrix and specified pattern");

private:
  typedef Matrix<ElementT, NumDimensions, IndexT, PatternT>
    self_t;
  typedef typename std::make_unsigned<IndexT>::type
    SizeType;
  typedef MatrixRefView<ElementT, NumDimensions, PatternT>
    MatrixRefView_t;
  typedef LocalMatrixRef<ElementT, NumDimensions, NumDimensions, PatternT>
    LocalRef_t;
  typedef LocalMatrixRef<
            const ElementT, NumDimensions, NumDimensions, PatternT>
    LocalRef_const_t;
  typedef PatternT
    Pattern_t;
  typedef GlobMem<ElementT, ddm::allocator::CollectiveAllocator<ElementT>>
    GlobMem_t;
  typedef DistributionSpec<NumDimensions>
    DistributionSpec_t;
  typedef SizeSpec<NumDimensions, typename PatternT::size_type>
    SizeSpec_t;
  typedef TeamSpec<NumDimensions, typename PatternT::index_type>
    TeamSpec_t;
  typedef std::array<typename PatternT::size_type, NumDimensions>
    Extents_t;
  typedef std::array<typename PatternT::index_type, NumDimensions>
    Offsets_t;

public:
  template<
    typename T_,
    dim_t NumDimensions1,
    dim_t NumDimensions2,
    class PatternT_ >
  friend class MatrixRef;
  template<
    typename T_,
    dim_t NumDimensions1,
    dim_t NumDimensions2,
    class PatternT_ >
  friend class LocalMatrixRef;

public:
  typedef ElementT                                              value_type;
  typedef typename PatternT::size_type                           size_type;
  typedef typename PatternT::index_type                    difference_type;
  typedef typename PatternT::index_type                         index_type;

  typedef GlobIter<      value_type, Pattern_t>                   iterator;
  typedef GlobIter<const value_type, Pattern_t>             const_iterator;

  typedef std::reverse_iterator<iterator>                 reverse_iterator;
  typedef std::reverse_iterator<const_iterator>     const_reverse_iterator;

  typedef GlobRef<      value_type>                              reference;
  typedef GlobRef<const value_type>                        const_reference;

  typedef GlobIter<      value_type, Pattern_t>                    pointer;
  typedef GlobIter<const value_type, Pattern_t>              const_pointer;

  typedef       ElementT *                                   local_pointer;
  typedef const ElementT *                             const_local_pointer;

// Public types as required by ddm container concept
public:
  /// Type specifying the view on local matrix elements.
  typedef LocalMatrixRef<
            ElementT, NumDimensions, NumDimensions, PatternT>
    local_type;

  /// Type specifying the view on const local matrix elements.
  typedef LocalMatrixRef<
            const ElementT, NumDimensions, NumDimensions, PatternT>
    const_local_type;

  /// The type of the pattern specifying linear iteration order and how
  /// elements are distribute to units.
  typedef PatternT                                            pattern_type;

  /// Type of views on matrix elements such as sub-matrices, row- and
  /// column vectors.
  template <dim_t NumViewDim>
    using view_type =
          MatrixRef<ElementT, NumDimensions, NumViewDim, PatternT>;

  /// Type of views on matrix elements such as sub-matrices, row- and
  /// column vectors.
  template <dim_t NumViewDim>
    using const_view_type =
          MatrixRef<const ElementT, NumDimensions, NumViewDim, PatternT>;

public:
  /// Local view proxy object.
  local_type local;

public:

  typedef std::integral_constant<dim_t, NumDimensions>                rank;

  static constexpr dim_t ndim() {
    return NumDimensions;
  }

public:
  /**
   * Default constructor, for delayed allocation.
   *
   * Sets the associated team to DART_TEAM_NULL for global matrix instances
   * that are declared before \c ddm::Init().
   */
  Matrix(
    Team & team = ddm::Team::Null());

  /**
   * Constructor, creates a new instance of Matrix.
   */
  Matrix(
    const SizeSpec_t         & ss,
    const DistributionSpec_t & ds  = DistributionSpec_t(),
    Team                     & t   = ddm::Team::All(),
    const TeamSpec_t         & ts  = TeamSpec_t());

  /**
   * Constructor, creates a new instance of Matrix from a pattern instance.
   */
  Matrix(
    const PatternT & pat);

  /**
   * Constructor, creates a new instance of Matrix.
   */
  Matrix(
    /// Number of elements
    size_t nelem,
    /// Team containing all units operating on the Matrix instance
    Team & t = ddm::Team::All())
  : Matrix(PatternT(nelem, t))
  { }

  /**
   * Constructor, creates a new instance of Matrix
   * of given extents.
   */
  template<typename... Args>
  Matrix(SizeType arg, Args... args)
  : Matrix(PatternT(arg, args... ))
  { }

  /**
   * Copy-constructor, deleted.
   */
  Matrix(const self_t &) = delete;

  /**
   * Destructor, frees underlying memory.
   */
  ~Matrix();

  /**
   * View at block at given global block coordinates.
   */
  view_type<NumDimensions> block(
    const std::array<index_type, NumDimensions> & block_gcoords);

  /**
   * View at block at given global block offset.
   */
  view_type<NumDimensions> block(
    index_type block_gindex);

#if 0
  /**
   * View at block at given global block offset with halo region.
   */
  halo_view_type<NumDimensions> block(
    index_type                            block_gindex,
    const ddm::HaloSpec<NumDimensions> & halospec);
#endif

  /**
   * Explicit allocation of matrix elements, used for delayed allocation
   * of default-constructed Matrix instance.
   *
   * \see  DDMContainerConcept
   */
  bool allocate(
    const SizeSpec_t         & sizespec,
    const DistributionSpec_t & distribution,
    const TeamSpec_t         & teamspec,
    ddm::Team               & team = ddm::Team::All()
  );

  /**
   * Allocation and distribution of matrix elements as specified by a given
   * Pattern instance.
   */
  bool allocate(
    const PatternT & pattern
  );

  /**
   * Explicit deallocation of matrix elements, called implicitly in
   * destructor and team deallocation.
   *
   * \see  DDMContainerConcept
   */
  void deallocate();

  Team                      & team();

  constexpr size_type         size()                const noexcept;
  constexpr size_type         local_size()          const noexcept;
  constexpr size_type         local_capacity()      const noexcept;
  constexpr size_type         extent(dim_t dim)     const noexcept;
  constexpr Extents_t         extents()             const noexcept;
  constexpr index_type        offset(dim_t dim)     const noexcept;
  constexpr Offsets_t         offsets()             const noexcept;
  constexpr bool              empty()               const noexcept;

  /**
   * Synchronize units associated with the matrix.
   *
   * \see  DDMContainerConcept
   */
  void                       barrier()              const;

  /**
   * The pattern used to distribute matrix elements to units in its
   * associated team.
   *
   * \see  DDMContainerConcept
   */
  constexpr const Pattern_t & pattern()             const;

  /**
   * Iterator referencing first matrix element in global index space.
   *
   * \see  DDMContainerConcept
   */
                  iterator    begin()        noexcept;

  /**
   * Iterator referencing first matrix element in global index space.
   *
   * \see  DDMContainerConcept
   */
  constexpr const_iterator    begin()  const noexcept;

  /**
   * Iterator referencing past the last matrix element in global index
   * space.
   *
   * \see  DDMContainerConcept
   */
                  iterator    end()          noexcept;

  /**
   * Iterator referencing past the last matrix element in global index
   * space.
   *
   * \see  DDMContainerConcept
   */
  constexpr const_iterator    end()    const noexcept;

  /**
   * Pointer to first element in local range.
   *
   * \see  DDMContainerConcept
   */
                  ElementT *  lbegin()       noexcept;

  /**
   * Pointer to first element in local range.
   *
   * \see  DDMContainerConcept
   */
  constexpr const ElementT *  lbegin() const noexcept;

  /**
   * Pointer past final element in local range.
   *
   * \see  DDMContainerConcept
   */
                  ElementT *  lend()         noexcept;

  /**
   * Pointer past final element in local range.
   *
   * \see  DDMContainerConcept
   */
  constexpr const ElementT *  lend()   const noexcept;

  /**
   * Subscript operator, returns a submatrix reference at given offset
   * in global element range.
   */
  constexpr const_view_type<NumDimensions-1> operator[](
    size_type n       ///< Offset in highest matrix dimension.
  ) const;

  /**
   * Subscript operator, returns a submatrix reference at given offset
   * in global element range.
   */
  view_type<NumDimensions-1> operator[](
    size_type n       ///< Offset in highest matrix dimension.
  );

  template<dim_t SubDimension>
  view_type<NumDimensions> sub(
    size_type n,      ///< Offset of the sub-range.
    size_type range   ///< Width of the sub-range.
  );

  /**
   * Projection to given offset in a sub-dimension.
   *
   * \see  row
   * \see  col
   */
  template<dim_t SubDimension>
  view_type<NumDimensions-1> sub(
    size_type n       ///< Offset in selected dimension.
  );

  /**
   * Local proxy object representing a view consisting of matrix elements
   * that are located in the active unit's local memory.
   */
  local_type sub_local() noexcept;

  /**
   * Projection to given offset in first sub-dimension (column), same as
   * \c sub<0>(n).
   *
   * \returns  A \c MatrixRef object representing the nth column
   *
   * \see  sub
   * \see  row
   */
  view_type<NumDimensions-1> col(
    size_type n       ///< Column offset.
  );

  /**
   * Projection to given offset in second sub-dimension (rows), same as
   * \c sub<1>(n).
   *
   * \returns  A \c MatrixRef object representing the nth row
   *
   * \see  sub
   * \see  col
   */
  view_type<NumDimensions-1> row(
    size_type n       ///< Row offset.
  );

  /**
   * Create a view representing the matrix slice within a column
   * range.
   * Same as \c sub<1>(offset, extent).
   *
   * \returns  A matrix view
   *
   * \see  sub
   */
  view_type<NumDimensions> cols(
    size_type offset, ///< Offset of first column in range.
    size_type range   ///< Number of columns in the range.
  );

  /**
   * Create a view representing the matrix slice within a row
   * range.
   * Same as \c sub<0>(offset, extent).
   *
   * \returns  A matrix view
   *
   * \see  sub
   */
  view_type<NumDimensions> rows(
    size_type n,      ///< Offset of first row in range.
    size_type range   ///< Number of rows in the range.
  );

  /**
   * Fortran-style subscript operator.
   * As an example, the operation \c matrix(i,j) is equivalent to
   * \c matrix[i][j].
   *
   * \returns  A global reference to the element at the given global
   *           coordinates.
   */
  template<typename ... Args>
  reference at(
    Args... args      ///< Global coordinates
  );

  /**
   * Fortran-style subscript operator, alias for \c at().
   * As an example, the operation \c matrix(i,j) is equivalent to
   * \c matrix[i][j].
   *
   * \returns  A global reference to the element at the given global
   *           coordinates.
   * \see  at
   */
  template<typename... Args>
  reference operator()(
    Args... args      ///< Global coordinates
  );

  /**
   * Whether the element at a global, canonical offset in the matrix
   * is local to the active unit.
   *
   * \see  DDMContainerConcept
   */
  constexpr bool is_local(
    size_type g_pos   ///< Canonical offset in global index space.
  ) const;

  /**
   * Whether the element at a global, canonical offset in a specific
   * dimension of the matrix is local to the active unit.
   */
  template<dim_t Dimension>
  constexpr bool is_local(
    size_type g_pos   ///< Linear offset in the selected dimension.
  ) const;

  template <int level>
  ddm::HView<self_t, level> hview();

  /**
   * Conversion operator to type \c MatrixRef.
   */
  operator
    MatrixRef<ElementT, NumDimensions, NumDimensions, PatternT> ();
    
private:
  /// Team containing all units that collectively instantiated the
  /// Matrix instance
  ddm::Team                 * _team = nullptr;
  /// Capacity (total number of elements) of the matrix
  size_type                    _size;
  /// Number of local elements in the array
  size_type                    _lsize;
  /// Number allocated local elements in the array
  size_type                    _lcapacity;
  /// Global pointer to initial element in the array
  pointer                      _begin;
  /// The matrix elements' distribution pattern
  Pattern_t                    _pattern;
  /// Global memory allocation and -access
  GlobMem_t                  * _glob_mem;
  /// Native pointer to first local element in the array
  ElementT                   * _lbegin;
  /// Native pointer past last local element in the array
  ElementT                   * _lend;
  /// Proxy instance for applying a view, e.g. in subscript operator
  view_type<NumDimensions>     _ref;
  
 public:

    bool operator< (const Matrix& MatrixObj) const
    {
        if(MatrixObj._team < this->_team)
            return true;
    };
};

/**
 * Template alias for ddm::Matrix with the same default template
 * arguments
 *
 * \see Matrix
 */
template <
  typename T,
  dim_t    NumDimensions,
  typename IndexT   = ddm::default_index_t,
  class    PatternT = Pattern<NumDimensions, ROW_MAJOR, IndexT> >
using NArray = ddm::Matrix<T, NumDimensions, IndexT, PatternT>;

}  // namespace ddm

#include "../ddm/matrix/internal/Matrix-inl.h"

#endif  // DDM__MATRIX_H_INCLUDED
