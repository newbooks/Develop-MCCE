#ifndef DDM__VIEW__VIEW_MOD_H__INCLUDED
#define DDM__VIEW__VIEW_MOD_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"
#include "../../ddm/Iterator.h"

#include "../../ddm/util/UniversalMember.h"

#include "../../ddm/view/IndexSet.h"
#include "../../ddm/view/ViewTraits.h"

#include "../../ddm/view/Local.h"
#include "../../ddm/view/Global.h"
#include "../../ddm/view/Origin.h"
#include "../../ddm/view/Domain.h"
#include "../../ddm/view/Apply.h"


namespace ddm {

/**
 * \defgroup  DDMViewExpressionConcept  Multidimensional View Expressions
 *
 * \ingroup DDMViewConcept
 * \{
 * \par Description
 *
 * Implementing view modifier chain as combination of command pattern
 * and chain of responsibility pattern.
 * For now, only compile-time projections/slices are supported such as:
 *
 * \code
 *   sub<0>(10,20).sub<1>(30,40)
 * \endcode
 *
 * but not run-time projections/slices like:
 *
 * \code
 *   sub(0, { 10,20 }).sub(1, { 30,40 })
 * \endcode
 *
 * \par Implementation Notes
 *
 * A view composition is a chained application of view modifier types
 * that depend on the type of their predecessor in the chain.
 *
 * Example:
 *
 * \code
 *  sub<0>(2).sub<1>(3,4)
 *  :         :
 *  |         |
 *  |         '--> ViewSubMod<0, ViewSubMod<-1, ViewOrigin> >
 *  |                            '------------.-----------'
 *  |                                         '--> parent
 *  '--> ViewSubMod<-1, ViewOrigin >
 *                      '----.---'
 *                           '--> parent
 * \endcode
 *
 * Consequently, specific ViewMod types are defined for every modifier
 * category.
 *
 * \}
 *
 *
 * \note
 *
 * As an alternative, all view modifications could be stored in command
 * objects of a single ViewMod type. Expressions then could not be
 * evalated at compile-time, however.
 *
 * However, View modifier types should subclass a common ViewMod base
 * class - or vice versa, following the policy pattern with the
 * operation specified as policy:
 *
 * \code
 *   template <dim_t DimDiff, class DomainType>
 *   class ViewMod : DomainType
 *   {
 *      // ...
 *   }
 * \endcode
 *
 * or:
 *
 * \code
 *   template <dim_t DimDiff, class ViewModOperation>
 *   class ViewMod : ViewModOperation
 *   {
 *      // ...
 *   }
 *
 *   class ViewModSubOperation;
 *   // defines
 *   // - sub<N>(...)
 *   // - view_mod_op() { return sub<N>(...); }
 *
 *   ViewMod<0, ViewModSubOperation> view_sub(initializer_list);
 *   // - either calls view_mod_op(initializer_list) in constructor
 *   // - or provides method sub<N>(...) directly
 * \endcode
 *
 * \todo Eventually, these probably are not public definitions.
 *       Move to namespace internal.
 *       Define dereference operator*() for view types, delegating to
 *       domain::operator* recursively.
 */


#ifndef DOXYGEN

// ------------------------------------------------------------------------
// Forward-declarations
// ------------------------------------------------------------------------

template <
  dim_t NDim = 1>
class ViewOrigin;

template <
  class ViewModType,
  class DomainType >
class ViewModBase;

template <
  class DomainType = ViewOrigin<1> >
class ViewLocalMod;

template <
  class DomainType = ViewOrigin<1>,
  dim_t SubDim     = 0 >
class ViewSubMod;

template <
  class DomainType = ViewOrigin<1> >
class ViewGlobalMod;


// --------------------------------------------------------------------
// ViewOrigin
// --------------------------------------------------------------------

/**
 * Monotype for the logical symbol that represents a view origin.
 */
template <dim_t NDim>
class ViewOrigin
{
  typedef ViewOrigin self_t;

 public:
  typedef ddm::default_index_t                                 index_type;
  typedef self_t                                               domain_type;
  typedef IndexSetIdentity<self_t>                          index_set_type;

 public:
  typedef std::integral_constant<bool, false>   is_local;
  typedef std::integral_constant<dim_t, NDim>   rank;

 private:
  std::array<index_type, NDim>                  _extents    = { };
  index_set_type                                _index_set;
 public:
  constexpr ViewOrigin()               = delete;
  constexpr ViewOrigin(self_t &&)      = default;
  constexpr ViewOrigin(const self_t &) = default;
  ~ViewOrigin()                        = default;
  self_t & operator=(self_t &&)        = default;
  self_t & operator=(const self_t &)   = default;

  constexpr explicit ViewOrigin(
      std::initializer_list<index_type> extents)
  : _extents(extents)
  , _index_set(*this)
  { }

  constexpr const domain_type & domain() const {
    return *this;
  }

  inline domain_type & domain() {
    return *this;
  }

  constexpr const index_set_type & index_set() const {
    return _index_set;
  }

  constexpr bool operator==(const self_t & rhs) const {
    return (this == &rhs);
  }
  
  constexpr bool operator!=(const self_t & rhs) const {
    return !(*this == rhs);
  }

  template <dim_t ExtentDim = 0>
  constexpr index_type extent() const {
    return _extents[ExtentDim];
  }
  
  constexpr index_type size() const {
    return _size<0>();
  }

 private:
  template <dim_t SizeDim = 0>
  constexpr index_type _size() const {
    return extent<SizeDim>() +
             (SizeDim < NDim
               ? _size<SizeDim + 1>()
               : 0);
  }
};

template <dim_t NDim>
struct view_traits<ViewOrigin<NDim>> {
  typedef ViewOrigin<NDim>                                       origin_type;
  typedef ViewOrigin<NDim>                                       domain_type;
  typedef ViewOrigin<NDim>                                        image_type;
  typedef typename ViewOrigin<NDim>::index_type                   index_type;
  typedef typename ViewOrigin<NDim>::index_set_type           index_set_type;

  typedef std::integral_constant<bool, false>                  is_projection;
  typedef std::integral_constant<bool, true>                   is_view;
  typedef std::integral_constant<bool, true>                   is_origin;
  typedef std::integral_constant<bool, false>                  is_local;

  typedef std::integral_constant<dim_t, NDim>                  rank;
};


// ------------------------------------------------------------------------
// ViewModBase
// ------------------------------------------------------------------------

template <
  class ViewModType,
  class DomainType >
class ViewModBase {
  typedef ViewModBase<ViewModType, DomainType> self_t;
//typedef typename std::remove_reference<DomainType>::type domain_value_type;
 public:
  typedef DomainType                                             domain_type;
  typedef typename view_traits<domain_type>::origin_type         origin_type;
  typedef typename view_traits<domain_type>::index_type           index_type;
  typedef typename origin_type::value_type                        value_type;

  typedef std::integral_constant<dim_t, domain_type::rank::value>       rank;

 protected:
  ddm::UniversalMember<domain_type> _domain;

  ViewModType & derived() {
    return static_cast<ViewModType &>(*this);
  }
  const ViewModType & derived() const {
    return static_cast<const ViewModType &>(*this);
  }

  /**
   * Constructor, creates a view on a given domain.
   */
  constexpr explicit ViewModBase(domain_type && domain)
  : _domain(std::forward<domain_type>(domain))
  { }

  /**
   * Constructor, creates a view on a given domain.
   */
  constexpr explicit ViewModBase(const domain_type & domain)
  : _domain(domain)
  { }

  constexpr ViewModBase()               = delete;
  constexpr ViewModBase(const self_t &) = delete;
  self_t & operator=(const self_t &)    = delete;

 public:
  constexpr ViewModBase(self_t &&)      = default;
  self_t & operator=(self_t &&)         = default;

  constexpr const domain_type & domain() const {
    return _domain;
  }

  constexpr bool operator==(const ViewModType & rhs) const {
    return &derived() == &rhs;
  }
  
  constexpr bool operator!=(const ViewModType & rhs) const {
    return !(derived() == rhs);
  }

  constexpr bool is_local() const {
    return view_traits<ViewModType>::is_local::value;
  }

  constexpr index_type size() const {
    return ddm::index(derived()).size();
  }
};


// ------------------------------------------------------------------------
// ViewLocalMod
// ------------------------------------------------------------------------

template <
  class DomainType >
struct view_traits<ViewLocalMod<DomainType> > {
  typedef DomainType                                           domain_type;
  typedef typename view_traits<domain_type>::origin_type       origin_type;
  typedef typename view_traits<domain_type>::pattern_type     pattern_type;
  typedef typename domain_type::local_type                      image_type;
  typedef ViewLocalMod<DomainType>                              local_type;
  typedef domain_type                                          global_type;

  typedef typename DomainType::index_type                       index_type;
  typedef ddm::IndexSetLocal< ViewLocalMod<DomainType> >   index_set_type;

  typedef std::integral_constant<bool, false>                is_projection;
  typedef std::integral_constant<bool, true>                 is_view;
  typedef std::integral_constant<bool, false>                is_origin;
  typedef std::integral_constant<bool, true>                 is_local;

  typedef std::integral_constant<dim_t, DomainType::rank::value> rank;
};

template <
  class DomainType >
class ViewLocalMod
: public ViewModBase< ViewLocalMod<DomainType>, DomainType > {
 public:
  typedef DomainType                                           domain_type;
  typedef typename view_traits<DomainType>::origin_type        origin_type;
  typedef typename domain_type::local_type                      image_type;
  typedef typename DomainType::index_type                       index_type;
 private:
  typedef ViewLocalMod<DomainType>                                  self_t;
  typedef ViewModBase< ViewLocalMod<DomainType>, DomainType >       base_t;
 public:
  typedef ddm::IndexSetLocal< ViewLocalMod<DomainType> >   index_set_type;
  typedef self_t                                                local_type;
  typedef typename domain_type::global_type                    global_type;

  typedef std::integral_constant<bool, true>                      is_local;

  typedef decltype(ddm::begin(ddm::local(
                std::declval<
                  typename std::add_lvalue_reference<origin_type>::type >()
              )))
    iterator;

 private:
  index_set_type  _index_set;
 public:
  constexpr ViewLocalMod()               = delete;
  constexpr ViewLocalMod(self_t &&)      = default;
  constexpr ViewLocalMod(const self_t &) = delete;
  ~ViewLocalMod()                        = default;
  self_t & operator=(self_t &&)          = default;
  self_t & operator=(const self_t &)     = delete;

  /**
   * Constructor, creates a view on a given domain.
   */
  constexpr explicit ViewLocalMod(
    domain_type && domain)
  : base_t(std::forward<domain_type>(domain))
  , _index_set(*this)
  { }

  /**
   * Constructor, creates a view on a given domain.
   */
  constexpr explicit ViewLocalMod(
    const DomainType & domain)
  : base_t(domain)
  , _index_set(*this)
  { }

  constexpr bool operator==(const self_t & rhs) const {
    return (this == &rhs ||
             ( base_t::operator==(rhs) &&
               _index_set == rhs._index_set ) );
  }

  constexpr bool operator!=(const self_t & rhs) const {
    return not (*this == rhs);
  }

  constexpr iterator begin() const {
    return ddm::begin(
             ddm::local(
               ddm::origin(
                 *this
               )
             )
           )
#if 1
         + _index_set.pre()[
             _index_set.first()
           ];
#else
         + ddm::index(ddm::local(ddm::domain(*this))).pre()[
             *ddm::begin(ddm::index(ddm::local(ddm::domain(*this))))
           ];
#endif
  }

  constexpr iterator end() const {
    return ddm::begin(
             ddm::local(
               ddm::origin(
                 *this
               )
             )
           )
#if 1
         + _index_set.pre()[
             _index_set.last()
           ] + 1;
#else
         + ddm::index(ddm::local(ddm::domain(*this))).pre()[
             *(ddm::begin(ddm::index(ddm::local(ddm::domain(*this))))
                 + ddm::index(ddm::local(ddm::domain(*this))).size()
                 - 1 )
           ] + 1;
#endif
  }

  constexpr auto operator[](int offset) const
  -> decltype(*(ddm::begin(
                ddm::local(ddm::origin(
                  std::declval<
                    typename std::add_lvalue_reference<domain_type>::type >()
                ))))) {
    return *(this->begin() + offset);
  }

  constexpr const local_type & local() const {
    return *this;
  }

  inline local_type & local() {
    return *this;
  }

  constexpr const global_type & global() const {
    return ddm::global(ddm::domain(*this));
  }

  inline global_type & global() {
    return ddm::global(ddm::domain(*this));
  }

  constexpr const index_set_type & index_set() const {
    return _index_set;
  }
};


// ------------------------------------------------------------------------
// ViewSubMod
// ------------------------------------------------------------------------

template <
  class DomainType,
  dim_t SubDim >
struct view_traits<ViewSubMod<DomainType, SubDim> > {
  typedef DomainType                                           domain_type;
  typedef typename view_traits<domain_type>::origin_type       origin_type;
  typedef typename view_traits<domain_type>::pattern_type     pattern_type;
  typedef ViewSubMod<DomainType, SubDim>                        image_type;
  typedef ViewSubMod<DomainType, SubDim>                        local_type;
  typedef ViewSubMod<DomainType, SubDim>                       global_type;

  typedef typename DomainType::index_type                       index_type;
  typedef ddm::IndexSetSub<ViewSubMod<DomainType, SubDim>> index_set_type;

  typedef std::integral_constant<bool, false>                is_projection;
  typedef std::integral_constant<bool, true>                 is_view;
  typedef std::integral_constant<bool, false>                is_origin;
  typedef std::integral_constant<bool,
    view_traits<domain_type>::is_local::value >              is_local;

  typedef std::integral_constant<dim_t, DomainType::rank::value> rank;
};


template <
  class DomainType,
  dim_t SubDim >
class ViewSubMod
: public ViewModBase<
           ViewSubMod<DomainType, SubDim>,
           DomainType >
{
 public:
  typedef DomainType                                             domain_type;
  typedef typename view_traits<DomainType>::index_type            index_type;
 private:
  typedef ViewSubMod<DomainType, SubDim>                              self_t;
  typedef ViewModBase< ViewSubMod<DomainType, SubDim>, DomainType >   base_t;
 public:
  typedef ddm::IndexSetSub< ViewSubMod<DomainType, SubDim> > index_set_type;
  typedef ViewLocalMod<self_t>                                    local_type;
  typedef self_t                                                 global_type;

  typedef std::integral_constant<bool, false>                       is_local;

  typedef decltype(ddm::begin(
                     std::declval<
                       typename std::add_lvalue_reference<domain_type>::type
                     >() ))
    iterator;

 private:
  index_type     _begin_idx;
  index_type     _end_idx;
  index_set_type _index_set;

 public:
  constexpr ViewSubMod()               = delete;
  constexpr ViewSubMod(self_t &&)      = default;
  constexpr ViewSubMod(const self_t &) = delete;
  ~ViewSubMod()                        = default;
  self_t & operator=(self_t &&)        = default;
  self_t & operator=(const self_t &)   = delete;

  constexpr ViewSubMod(
    domain_type && domain,
    index_type     begin,
    index_type     end)
  : base_t(std::forward<domain_type>(domain))
  , _begin_idx(begin)
  , _end_idx(end)
  , _index_set(*this, begin, end)
  { }

  constexpr ViewSubMod(
    domain_type  & domain,
    index_type     begin,
    index_type     end)
  : base_t(domain)
  , _begin_idx(begin)
  , _end_idx(end)
  , _index_set(*this, begin, end)
  { }

  constexpr iterator begin() const {
    return ddm::begin(ddm::domain(*this)) +
             *ddm::begin(ddm::index(*this));
  }

  constexpr iterator end() const {
    return ddm::begin(ddm::domain(*this)) +
             *ddm::end(ddm::index(*this));
  }

  constexpr auto operator[](int offset) const
  -> decltype(*(ddm::begin(
                  std::declval<
                    typename std::add_lvalue_reference<domain_type>::type
                  >() ))) {
    return *(this->begin() + offset);
  }

  constexpr const index_set_type & index_set() const {
    return _index_set;
  }

  constexpr local_type local() const {
    return local_type(*this);
  }
};


// ------------------------------------------------------------------------
// ViewGlobalMod
// ------------------------------------------------------------------------

template <
  class DomainType >
struct view_traits<ViewGlobalMod<DomainType> > {
  typedef DomainType                                           domain_type;
  typedef typename view_traits<domain_type>::origin_type       origin_type;
  typedef typename view_traits<domain_type>::pattern_type     pattern_type;
  typedef typename domain_type::global_type                     image_type;
  typedef typename domain_type::local_type                      local_type;
  typedef ViewGlobalMod<DomainType>                            global_type;

  typedef typename DomainType::index_type                       index_type;
  typedef ddm::IndexSetLocal< ViewLocalMod<DomainType> >   index_set_type;

  typedef std::integral_constant<bool, false>                is_projection;
  typedef std::integral_constant<bool, true>                 is_view;
  typedef std::integral_constant<bool, false>                is_origin;
  typedef std::integral_constant<bool, false>                is_local;
};

template <
  class DomainType >
class ViewGlobalMod
: public ViewModBase< ViewGlobalMod<DomainType>, DomainType >
{
 public:
  typedef DomainType                                           domain_type;
  typedef typename view_traits<DomainType>::origin_type        origin_type;
  typedef typename domain_type::global_type                     image_type;
  typedef typename DomainType::index_type                       index_type;
 private:
  typedef ViewGlobalMod<DomainType>                                 self_t;
  typedef ViewModBase< ViewLocalMod<DomainType>, DomainType >       base_t;
 public:
  typedef ddm::IndexSetGlobal< ViewGlobalMod<DomainType> > index_set_type;
  typedef self_t                                               global_type;
  typedef typename domain_type::local_type                      local_type;

  typedef std::integral_constant<bool, false>                     is_local;

 private:
  index_set_type  _index_set;
 public:
  constexpr ViewGlobalMod()               = delete;
  constexpr ViewGlobalMod(self_t &&)      = default;
  constexpr ViewGlobalMod(const self_t &) = default;
  ~ViewGlobalMod()                        = default;
  self_t & operator=(self_t &&)           = default;
  self_t & operator=(const self_t &)      = default;

  /**
   * Constructor, creates a view on a given domain.
   */
  constexpr explicit ViewGlobalMod(
    const domain_type & domain)
  : base_t(domain)
  , _index_set(*this)
  { }

  /**
   * Constructor, creates a view on a given domain.
   */
  constexpr explicit ViewGlobalMod(
    domain_type && domain)
  : base_t(std::forward<domain_type>(domain))
  , _index_set(*this)
  { }

  constexpr auto begin() const
  -> decltype(ddm::begin(ddm::global(ddm::domain(*this)))) {
    return ddm::begin(
             ddm::global(
               ddm::domain(
                 *this)));
  }

  constexpr auto end() const
  -> decltype(ddm::end(ddm::global(ddm::domain(*this)))) {
    return ddm::begin(
             ddm::global(
               ddm::domain(
                 *this)))
           + *ddm::end(ddm::index(ddm::domain(*this)));
  }

  constexpr auto operator[](int offset) const
  -> decltype(*(ddm::begin(
                 ddm::global(ddm::domain(*this))))) {
    return *(this->begin() + offset);
  }

  constexpr const local_type & local() const {
    // if any parent domain is local, it will return *this
    // and in effect eliminate ddm::global( ... ddm::local( ... ))
    return ddm::local(ddm::domain(*this));
  }

  inline local_type & local() {
    // if any parent domain is local, it will return *this
    // and in effect eliminate ddm::global( ... ddm::local( ... ))
    return ddm::local(ddm::domain(*this));
  }

  constexpr const global_type & global() const {
    return *this;
  }

  inline global_type & global() {
    return *this;
  }

  constexpr const index_set_type & index_set() const {
    return _index_set;
  }
};

#endif // DOXYGEN

} // namespace ddm

#endif // DDM__VIEW__VIEW_MOD_H__INCLUDED
