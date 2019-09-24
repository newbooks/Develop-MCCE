#ifndef DDM__MATRIX__MATRIX_REF_VIEW_H_INCLUDED
#define DDM__MATRIX__MATRIX_REF_VIEW_H_INCLUDED

#include "../dart-impl/dart.h"

#include "../../ddm/Team.h"
#include "../../ddm/Pattern.h"
#include "../../ddm/GlobRef.h"
#include "../../ddm/HView.h"

#include "../../ddm/iterator/GlobIter.h"


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
 * Stores information needed by subscripting and subdim selection.
 * A new \c MatrixRefView instance is created once for every dimension in
 * multi-subscripting.
 *
 * \ingroup DDMMatrixConcept
 */
template <
  typename T,
  dim_t NumDimensions,
  class PatternT =
    TilePattern<NumDimensions, ROW_MAJOR, ddm::default_index_t> >
class MatrixRefView
{
 public:
  typedef typename PatternT::index_type             index_type;

 private:
  /// The view's next unspecified dimension, initialized with 0.
  dim_t                                             _dim       = 0;
  /// The matrix referenced by the view.
  Matrix<T, NumDimensions, index_type, PatternT>  * _mat;
  /// Coordinates of a single referenced element if view references fully
  /// specified coordinates.
  ::std::array<index_type, NumDimensions>           _coord     = {{  }};
  /// View offset and extents in global index range.
  ViewSpec<NumDimensions, index_type>               _viewspec;
  /// View offset and extents in local index range.
  ViewSpec<NumDimensions, index_type>               _l_viewspec;

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
    class PatternT_ >
  friend class MatrixRefView;
  template<
    typename T_,
    dim_t NumDimensions1,
    dim_t NumDimensions2,
    class PatternT_ >
  friend class LocalMatrixRef;
  template<
    typename T_,
    dim_t NumDimensions1,
    typename IndexT_,
    class PatternT_ >
  friend class Matrix;

  MatrixRefView<T, NumDimensions, PatternT>();

  template <class T_>
  MatrixRefView<T, NumDimensions, PatternT>(
    const MatrixRefView<T_, NumDimensions, PatternT> & other);

  template <class T_>
  MatrixRefView<T, NumDimensions, PatternT>(
    Matrix<T_, NumDimensions, index_type, PatternT> * matrix);

  GlobRef<T>       global_reference();
  GlobRef<const T> global_reference() const;

  GlobRef<T>       global_reference(
    const ::std::array<typename PatternT::index_type, NumDimensions> & coords
  );
  GlobRef<const T> global_reference(
    const ::std::array<typename PatternT::index_type, NumDimensions> & coords
  ) const;
};

} // namespace ddm

#include "../../ddm/matrix/internal/MatrixRefView-inl.h"

#endif  // DDM__MATRIX__MATRIX_REF_VIEW_H_INCLUDED
