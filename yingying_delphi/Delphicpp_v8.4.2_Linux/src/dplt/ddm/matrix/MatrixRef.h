#ifndef DDM__MATRIX__MATRIX_REF_H__
#define DDM__MATRIX__MATRIX_REF_H__

#include "../../ddm/Matrix.h"

#include "../dart-impl/dart.h"

#include "../../ddm/Team.h"
#include "../../ddm/Pattern.h"
#include "../../ddm/GlobRef.h"
#include "../../ddm/HView.h"

#include "../../ddm/view/Sub.h"

#include "../../ddm/iterator/GlobIter.h"
#include "../../ddm/iterator/GlobViewIter.h"


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
  dim_t NumViewDim,
  class PatternT >
class MatrixRef;
/// Forward-declaration
template <
  typename T,
  dim_t NumDimensions,
  dim_t NumViewDim,
  class PatternT >
class LocalMatrixRef;

/**
 * A view on a referenced \c Matrix object, such as a dimensional
 * projection returned by \c Matrix::sub.
 *
 * TODO:
 * Projection order matrix.sub().local() is not fully implemented yet.
 * Currently only matrix.local().sub() is supported.
 *
 * \see DDMMatrixConcept
 *
 * \ingroup Matrix
 */
template <
  typename ElementT,
  dim_t NumDimensions,
  dim_t NumViewDim = NumDimensions,
  class PatternT =
    TilePattern<NumDimensions, ROW_MAJOR, ddm::default_index_t> >
class MatrixRef
{
 private:
  typedef MatrixRef<ElementT, NumDimensions, NumViewDim, PatternT>
    self_t;
  typedef PatternT
    Pattern_t;
  typedef typename PatternT::index_type
    Index_t;
  typedef typename PatternT::size_type
    Size_t;
  typedef MatrixRefView<ElementT, NumDimensions, PatternT>
    MatrixRefView_t;
  typedef CartesianIndexSpace<
            NumViewDim,
            PatternT::memory_order(),
            typename PatternT::index_type >
    IndexSpace_t;
  typedef std::array<typename PatternT::size_type, NumDimensions>
    Extents_t;
  typedef std::array<typename PatternT::index_type, NumDimensions>
    Offsets_t;

 public:
  typedef PatternT                                        pattern_type;
  typedef typename PatternT::index_type                     index_type;
  typedef ElementT                                          value_type;
  typedef IndexSpace_t                                index_space_type;

  typedef typename PatternT::size_type                       size_type;
  typedef typename PatternT::index_type                difference_type;

  typedef GlobViewIter<      value_type, PatternT>            iterator;
  typedef GlobViewIter<const value_type, PatternT>      const_iterator;

  typedef std::reverse_iterator<iterator>             reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  typedef GlobRef<      value_type>                          reference;
  typedef GlobRef<const value_type>                    const_reference;

  typedef GlobViewIter<      value_type, PatternT>             pointer;
  typedef GlobViewIter<const value_type, PatternT>       const_pointer;

  typedef LocalMatrixRef<
            ElementT, NumDimensions, NumDimensions, PatternT>
                                                            local_type;
  typedef LocalMatrixRef<
            const ElementT, NumDimensions, NumDimensions, PatternT>
                                                      const_local_type;

  typedef LocalMatrixRef<
            ElementT, NumDimensions, NumDimensions, PatternT>
                                                  local_reference_type;
  typedef LocalMatrixRef<
            const ElementT, NumDimensions, NumDimensions, PatternT>
                                            const_local_reference_type;

 public:
  template<
    typename T_,
    dim_t NumDimensions_,
    typename IndexT_,
    class PatternT_ >
  friend class Matrix;
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

  inline operator
    MatrixRef<ElementT, NumDimensions, NumViewDim-1, PatternT> && ();

public:
  typedef std::integral_constant<dim_t, NumViewDim>
    rank;

  static constexpr dim_t ndim() {
    return NumViewDim;
  }

public:

  MatrixRef<ElementT, NumDimensions, NumViewDim, PatternT>()
  {
    DDM_LOG_TRACE_VAR("MatrixRef<T,D,C>()", NumDimensions);
    DDM_LOG_TRACE_VAR("MatrixRef<T,D,C>()", NumViewDim);
  }

  template <class T_>
  MatrixRef<ElementT, NumDimensions, NumViewDim, PatternT>(
    const MatrixRef<T_, NumDimensions, NumViewDim+1, PatternT> & prev,
    index_type coord);

  constexpr const Team      & team()                const noexcept;

  constexpr size_type         size()                const noexcept;
  constexpr size_type         local_size()          const noexcept;
  constexpr size_type         local_capacity()      const noexcept;
  constexpr size_type         extent(dim_t dim)     const noexcept;
  constexpr Extents_t         extents()             const noexcept;
  constexpr index_type        offset(dim_t dim)     const noexcept;
  constexpr Offsets_t         offsets()             const noexcept;
  constexpr bool              empty()               const noexcept;

  inline    void              barrier()             const;

  /**
   * The pattern used to distribute matrix elements to units in its
   * associated team.
   */
  constexpr const PatternT  & pattern()             const noexcept;

  constexpr const_pointer     data()                const noexcept;
                  pointer     data()                      noexcept;

            iterator          begin()                     noexcept;
  constexpr const_iterator    begin()               const noexcept;

            iterator          end()                       noexcept;
  constexpr const_iterator    end()                 const noexcept;

  /// View representing elements in the active unit's local memory.
  inline          local_type   sub_local()                noexcept;

  /// Pointer to first element in local range.
  inline          ElementT   * lbegin()                   noexcept;
  /// Pointer to first element in local range.
  constexpr const ElementT   * lbegin()             const noexcept;
  /// Pointer past final element in local range.
  inline          ElementT   * lend()                     noexcept;
  /// Pointer past final element in local range.
  constexpr const ElementT   * lend()               const noexcept;

  /**
   * Subscript operator, returns a submatrix reference at given offset
   * in global element range.
   */
  MatrixRef<ElementT, NumDimensions, NumViewDim-1, PatternT>
    operator[](index_type n);

  /**
   * Subscript operator, returns a submatrix reference at given offset
   * in global element range.
   */
  constexpr MatrixRef<const ElementT, NumDimensions, NumViewDim-1, PatternT>
    operator[](index_type n) const;

  template<dim_t NumSubDimensions>
  MatrixRef<ElementT, NumDimensions, NumDimensions-1, PatternT>
  sub(size_type n);

  MatrixRef<ElementT, NumDimensions, NumDimensions-1, PatternT>
  col(size_type n);

  MatrixRef<ElementT, NumDimensions, NumDimensions-1, PatternT>
  row(size_type n);

  template<dim_t SubDimension>
  MatrixRef<ElementT, NumDimensions, NumDimensions, PatternT>
  sub(
    size_type n,
    size_type range);

  /**
   * Create a view representing the matrix slice within a row
   * range.
   * Same as \c sub<0>(offset, extent).
   *
   * \returns  A matrix view
   *
   * \see  sub
   */
  MatrixRef<ElementT, NumDimensions, NumDimensions, PatternT>
  rows(
    /// Offset of first row in range
    size_type n,
    /// Number of rows in the range
    size_type range);

  /**
   * Create a view representing the matrix slice within a column
   * range.
   * Same as \c sub<1>(offset, extent).
   *
   * \returns  A matrix view
   *
   * \see  sub
   */
  MatrixRef<ElementT, NumDimensions, NumDimensions, PatternT>
  cols(
    /// Offset of first column in range
    size_type offset,
    /// Number of columns in the range
    size_type range);

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
    /// Global coordinates
    Args... args);

  /**
   * Fortran-style subscript operator.
   * As an example, the operation \c matrix(i,j) is equivalent to
   * \c matrix[i][j].
   *
   * \returns  A global reference to the element at the given global
   *           coordinates.
   */
  reference at(
    /// Global coordinates
    const ::std::array<index_type, NumDimensions> & coords);

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
    /// Global coordinates
    Args... args);

  constexpr bool is_local(index_type n) const;

  template<dim_t Dimension>
  constexpr bool is_local(index_type n) const;

  constexpr const ViewSpec<NumDimensions, index_type> & viewspec() const {
    return _refview._viewspec;
  }

  template <int level>
  ddm::HView<Matrix<ElementT, NumDimensions, Index_t, PatternT>, level>
  inline hview();

 private:
  MatrixRefView<ElementT, NumDimensions, PatternT> _refview;
};

/**
 * A view on a referenced \c Matrix object, such as a dimensional
 * projection returned by \c Matrix::sub.
 * Partial Specialization for value deferencing.
 *
 * \ingroup Matrix
 */
template <
  typename ElementT,
  dim_t NumDimensions,
  class PatternT >
class MatrixRef< ElementT, NumDimensions, 0, PatternT >
{
 private:
   typedef MatrixRef<ElementT, NumDimensions, 0, PatternT> self_t;

 public:
  template<
    typename T_,
    dim_t NumDimensions_,
    typename IndexT_,
    class PatternT_ >
  friend class Matrix;
  template<
    typename T_,
    dim_t NumDimensions1,
    dim_t NumDimensions2,
    class PatternT_ >
  friend class MatrixRef;

 public:
  typedef PatternT                       pattern_type;
  typedef typename PatternT::index_type  index_type;
  typedef ElementT                       value_type;

  /**
   * Default constructor.
   */
  MatrixRef<ElementT, NumDimensions, 0, PatternT>()
  {
    DDM_LOG_TRACE_VAR("MatrixRef<T,D,0>()", NumDimensions);
  }

  /**
   * Copy constructor.
   */
  MatrixRef<ElementT, NumDimensions, 0, PatternT>(
    const self_t & other)
  : _refview(other._refview) {
    DDM_LOG_TRACE_VAR("MatrixRef<T,D,0>(other)", NumDimensions);
  }

  MatrixRef<ElementT, NumDimensions, 0, PatternT>(
    const MatrixRef<ElementT, NumDimensions, 1, PatternT> & previous,
    index_type coord);

  /**
   * TODO[TF] The following two functions don't seem to be implemented.
   */
  constexpr GlobRef<const ElementT> local_at(
    team_unit_t unit,
    index_type   elem) const;

  inline GlobRef<ElementT> local_at(
    team_unit_t unit,
    index_type   elem);

  constexpr bool is_local() const;

  constexpr ViewSpec<NumDimensions, index_type> &
  viewspec() const {
    return _refview._viewspec;
  }

  inline ViewSpec<NumDimensions, index_type> &
  viewspec() {
    return _refview._viewspec;
  }

  constexpr operator ElementT() const;
  constexpr operator GlobPtr<ElementT, PatternT>() const;

  /**
   * Assignment operator.
   */
  inline ElementT operator= (const ElementT & value);
  inline ElementT operator+=(const ElementT & value);
  inline ElementT operator+ (const ElementT & value);
  inline ElementT operator-=(const ElementT & value);
  inline ElementT operator- (const ElementT & value);
  inline ElementT operator*=(const ElementT & value);
  inline ElementT operator* (const ElementT & value);
  inline ElementT operator/=(const ElementT & value);
  inline ElementT operator/ (const ElementT & value);

 private:
  MatrixRefView<ElementT, NumDimensions, PatternT> _refview;
};

} // namespace ddm

#include "../../ddm/matrix/internal/MatrixRef-inl.h"

#endif // DDM__MATRIX__MATRIX_REF_H__
