#ifndef DDM__VIEW_H__INCLUDED
#define DDM__VIEW_H__INCLUDED

/**
 * \defgroup  DDMViewConcept  Multidimensional View Concept
 *
 * \ingroup DDMNDimConcepts
 * \{
 * \par Description
 *
 * Definitions for multidimensional view expressions.
 * A view expression consumes a view object (its origin) and returns a view
 * object that applies the expression's modification on the consumed origin.
 *
 * The result of a view expression satisfies the multidimensional Range
 * concept.
 *
 * \see DDMDimensionalConcept
 * \see DDMRangeConcept
 * \see DDMIteratorConcept
 * \see \c ddm::view_traits
 *
 * \par Terminology
 *
 * A \b View is a mapping between two index sets, from a \b Domain space to
 * an \b Image space in the view's codomain.
 * Views can be chained such that the image obtained from the application of
 * a view expression can again act as the domain of other views.
 * In effect, a view expression can be understood as a composite function on
 * an index set.
 * The \b Origin of a view is the first domain in the view chain that has not
 * been created by a view expression; simply put, a view origin is usually a
 * container.
 *
 * \par Expressions
 *
 * View Specifier            | Synopsis
 * ------------------------- | --------------------------------------------------
 * <tt>ddm::sub</tt>        | Subrange of domain in a specified dimension
 * <tt>ddm::intersect</tt>  | View from intersection of two domains
 * <tt>ddm::difference</tt> | View from difference of two domains
 * <tt>ddm::combine</tt>    | Composite view of two possibply unconnected domains
 * <tt>ddm::local</tt>      | Local subspace of domain
 * <tt>ddm::global</tt>     | Maps subspace to elements in global domain
 * <tt>ddm::apply</tt>      | Obtain image of domain view (inverse of \c domain)
 * <tt>ddm::domain</tt>     | Obtain domain of view image (inverse of \c apply)
 * <tt>ddm::origin</tt>     | Obtain the view origin (root domain)
 * <tt>ddm::blocks</tt>     | Decompose domain into blocks
 * <tt>ddm::block</tt>      | Subspace of decomposed domain in a specific block
 * <tt>ddm::index</tt>      | Returns a view's index set
 *
 * \par Examples
 *
 * \code
 * auto matrix_rect = ddm::sub<0>(10,20,
 *                    ddm::sub<1>(30,40,
 *                    matrix));
 *
 * auto matrix_rect_size = ddm::size(matrix_rect);
 * // -> 10x10 = 100
 *
 * auto matrix_rect_begin_gidx = ddm::index(ddm::begin(matrix_rect));
 * auto matrix_rect_end_gidx   = ddm::index(ddm::end(matrix_rect));
 *
 * for (auto elem : matrix_rect) {
 *   // ...
 * }
 *
 * \endcode
 *
 * \}
 */

#include "../ddm/view/Apply.h"
#include "../ddm/view/Block.h"
#include "../ddm/view/Chunked.h"
#include "../ddm/view/Domain.h"
#include "../ddm/view/Origin.h"
#include "../ddm/view/Global.h"
#include "../ddm/view/Local.h"
#include "../ddm/view/Remote.h"
#include "../ddm/view/Apply.h"
#include "../ddm/view/Sub.h"

#include "../ddm/view/SetUnion.h"
#include "../ddm/view/SetIntersect.h"

#include "../ddm/view/IndexSet.h"
#include "../ddm/view/IndexRange.h"
#include "../ddm/view/CartView.h"

#include "../ddm/view/MultiView.h"
#include "../ddm/view/StridedView.h"

#include "../ddm/view/ViewTraits.h"

#include "../ddm/view/ViewMod.h"
#include "../ddm/view/ViewBlocksMod.h"

#include "../ddm/Range.h"


namespace ddm {

template <
  class Iterator,
  class Sentinel >
class IteratorViewOrigin;

// Currently only supporting
//  - global iterators
//  - 1-dimensional IteratorViewOrigin types.

template <
  class Iterator,
  class Sentinel >
struct view_traits<IteratorViewOrigin<Iterator, Sentinel> > {
  typedef IteratorViewOrigin<Iterator, Sentinel>               domain_type;
  typedef IteratorViewOrigin<Iterator, Sentinel>               origin_type;
  typedef IteratorViewOrigin<Iterator, Sentinel>                image_type;

  // Uses container::local_type directly, e.g. ddm::LocalArrayRef:
//typedef typename ddm::view_traits<domain_type>::local_type   local_type;
  // Uses ViewLocalMod wrapper on domain, e.g. ViewLocalMod<ddm::Array>:
  typedef ViewLocalMod<domain_type>                             local_type;
  typedef ViewGlobalMod<domain_type>                           global_type;

  typedef typename Iterator::index_type                         index_type;
  typedef ddm::IndexSetIdentity< 
            IteratorViewOrigin<Iterator, Sentinel> >        index_set_type;

  typedef std::integral_constant<bool, false>                is_projection;
  typedef std::integral_constant<bool, true>                 is_view;
  typedef std::integral_constant<bool, true>                 is_origin;
  typedef std::integral_constant<bool, false>                is_local;
};

template <
  class Iterator,
  class Sentinel
//class DomainType = ddm::IteratorViewOriginOrigin<1>
>
class IteratorViewOrigin
: public ddm::IteratorRange<Iterator, Sentinel>
{
public:
  typedef typename Iterator::index_type                         index_type;
private:
  typedef IteratorViewOrigin<Iterator, Sentinel>  self_t;
  typedef IteratorRange<Iterator, Sentinel>       base_t;
public:
  typedef self_t                                               domain_type;
  typedef self_t                                               origin_type;
  typedef self_t                                                image_type;
  // Alternative: IteratorLocalView<self_t, Iterator, Sentinel>
  typedef typename view_traits<domain_type>::local_type         local_type;
  typedef typename view_traits<domain_type>::global_type       global_type;

  typedef typename Iterator::pattern_type                     pattern_type;
public:
  constexpr IteratorViewOrigin(Iterator begin, Iterator end)
  : base_t(std::move(begin), std::move(end)) {
  }

  constexpr const pattern_type & pattern() const {
    return this->begin().pattern();
  }

  constexpr local_type local() const {
    // for local_type: IteratorLocalView<self_t, Iterator, Sentinel>
    // return local_type(this->begin(), this->end());
    return local_type(*this);
  }

};


template <class Iterator, class Sentinel>
constexpr ddm::IteratorViewOrigin<Iterator, Sentinel>
make_view(Iterator begin, Sentinel end) {
  return ddm::IteratorViewOrigin<Iterator, Sentinel>(
           std::move(begin),
           std::move(end));
}

} // namespace ddm


#endif // DDM__VIEW_H__INCLUDED
