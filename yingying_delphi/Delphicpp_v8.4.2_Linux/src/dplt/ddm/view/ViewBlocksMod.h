#ifndef DDM__VIEW__VIEW_BLOCKS_MOD_H__INCLUDED
#define DDM__VIEW__VIEW_BLOCKS_MOD_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Range.h"
#include "../../ddm/Iterator.h"

#include "../../ddm/view/IndexSet.h"
#include "../../ddm/view/ViewTraits.h"
#include "../../ddm/view/ViewMod.h"

#include "../../ddm/view/Local.h"
#include "../../ddm/view/Global.h"
#include "../../ddm/view/Origin.h"
#include "../../ddm/view/Domain.h"
#include "../../ddm/view/Apply.h"
#include "../../ddm/view/Chunked.h"
#include "../../ddm/view/Sub.h"


namespace ddm {

#ifndef DOXYGEN

// ------------------------------------------------------------------------
// Forward-declarations
// ------------------------------------------------------------------------
//
template <
  class DomainType = ViewOrigin<1> >
class ViewBlocksMod;

// ------------------------------------------------------------------------
// ViewBlockMod
// ------------------------------------------------------------------------

template <
  class DomainType >
struct view_traits<ViewBlockMod<DomainType> > {
  typedef DomainType                                           domain_type;
  typedef typename view_traits<domain_type>::origin_type       origin_type;
  typedef typename view_traits<domain_type>::pattern_type     pattern_type;
  typedef ViewBlockMod<DomainType>                              image_type;
  typedef ViewBlockMod<DomainType>                              local_type;
  typedef ViewBlockMod<DomainType>                             global_type;

  typedef typename DomainType::index_type                       index_type;
  typedef ddm::IndexSetSub<ViewBlockMod<DomainType>>       index_set_type;

  typedef std::integral_constant<bool, false>                is_projection;
  typedef std::integral_constant<bool, true>                 is_view;
  typedef std::integral_constant<bool, false>                is_origin;
  typedef std::integral_constant<bool,
    view_traits<domain_type>::is_local::value >              is_local;

  typedef std::integral_constant<dim_t, 1>                            rank;
};


template <
  class DomainType >
class ViewBlockMod
// Actually just an adapter for block_idx -> sub(begin_idx, end_idx),
// should sublass
//
//   public ViewSubMod<DomainType, 0>
//
: public ViewModBase <
           ViewBlockMod<DomainType>,
           DomainType >
{
 public:
  typedef DomainType                                           domain_type;
  typedef typename view_traits<DomainType>::index_type          index_type;
 private:
  typedef ViewBlockMod<DomainType>                                  self_t;
  typedef ViewModBase< ViewBlockMod<DomainType>, DomainType >
                                                                    base_t;
 public:
//typedef ddm::IndexSetSub< ViewBlockMod<DomainType> >     index_set_type;
  typedef ddm::IndexSetSub< DomainType >                   index_set_type;
  typedef ViewLocalMod<self_t>                                  local_type;
  typedef self_t                                               global_type;

  typedef std::integral_constant<bool, false>                     is_local;

  typedef decltype(ddm::begin(
                std::declval<
                  typename std::add_lvalue_reference<domain_type>::type
                >() ))
    iterator;

 private:
  index_set_type _index_set;

 public:
  constexpr ViewBlockMod()               = delete;
  constexpr ViewBlockMod(self_t &&)      = default;
  constexpr ViewBlockMod(const self_t &) = delete;
  ~ViewBlockMod()                        = default;
  self_t & operator=(self_t &&)          = default;
  self_t & operator=(const self_t &)     = delete;

  constexpr ViewBlockMod(
    const domain_type & domain,
    index_type    block_idx)
  : base_t(domain)
  , _index_set(domain,
               block_first_gidx(domain, block_idx),
               block_final_gidx(domain, block_idx))
  { }

  /**
   * Constructor, creates a view on a block in the specified domain.
   */
  constexpr ViewBlockMod(
    domain_type     && domain,
    index_type         block_idx)
  : base_t(std::forward<domain_type>(domain))
  , _index_set(domain,
               block_first_gidx(domain, block_idx),
               block_final_gidx(domain, block_idx))
  { }

  constexpr auto begin() const
  -> decltype(ddm::begin(
                std::declval<
                  typename std::add_lvalue_reference<domain_type>::type
                >() )) {
    return this->domain().begin() +
              _index_set[0];
  }

  constexpr auto end() const
  -> decltype(ddm::begin(
                std::declval<
                  typename std::add_lvalue_reference<domain_type>::type
                >() )) {
    return this->domain().begin() +
               _index_set.last() + 1;
  }

  constexpr auto operator[](int offset) const
  -> decltype(*(ddm::begin(
                  std::declval<
                    typename std::add_lvalue_reference<domain_type>::type
                  >() ))) {
    return begin()[offset];
  }

  constexpr const index_set_type & index_set() const {
    return _index_set;
  }

  constexpr local_type local() const {
    return local_type(*this);
  }

 private:
  /// Block index of first element in view
  ///
  constexpr index_type block_first_gidx(
      const DomainType & vdomain,
      index_type         block_idx) const {
    //
    // TODO: If domain is local, use pattern().local_block(block_idx)
    //
    return std::max(
             ( // block viewspec (extents, offsets)
               ddm::index(vdomain)
                 .pattern().block(block_idx).offsets()[0]
             ),
             ddm::index(vdomain).first()
           )
           - ddm::index(vdomain).first();
  }

  /// Index past block index of last element in view:
  ///
  constexpr index_type block_final_gidx(
      const DomainType & vdomain,
      index_type         block_idx) const {
    //
    // TODO: If domain is local, use pattern().local_block(block_idx)
    //
    return std::min<index_type>(
             ddm::index(vdomain).last() + 1,
             ( // block viewspec (extents, offsets)
               ddm::index(vdomain)
                 .pattern().block(block_idx).offsets()[0]
             + ddm::index(vdomain)
                 .pattern().block(block_idx).extents()[0]
             )
           )
           - ddm::index(vdomain).first();
  }
};

// ------------------------------------------------------------------------
// ViewBlocksMod
// ------------------------------------------------------------------------

template <class ViewType>
constexpr ViewBlocksMod<ViewType>
blocks(const ViewType & domain) {
  return ViewBlocksMod<ViewType>(domain);
}

template <
  class    ViewType,
  typename ViewValueT =
             typename std::remove_const<
               typename std::remove_reference<ViewType>::type
             >::type
>
constexpr ViewBlocksMod<ViewValueT>
blocks(ViewType && domain) {
  return ViewBlocksMod<ViewValueT>(std::forward<ViewType>(domain));
}

template <
  class DomainType >
struct view_traits<ViewBlocksMod<DomainType> > {
  typedef DomainType                                           domain_type;
  typedef typename view_traits<domain_type>::origin_type       origin_type;
  typedef typename view_traits<domain_type>::pattern_type     pattern_type;
  typedef ViewBlocksMod<DomainType>                             image_type;
  typedef typename domain_type::local_type                      local_type;
  typedef ViewBlocksMod<DomainType>                            global_type;

  typedef typename DomainType::index_type                       index_type;
  typedef ddm::IndexSetBlocks<ViewBlocksMod<DomainType>>   index_set_type;

  typedef std::integral_constant<bool, false>                is_projection;
  typedef std::integral_constant<bool, true>                 is_view;
  typedef std::integral_constant<bool, false>                is_origin;
  typedef std::integral_constant<bool, false>                is_local;
};

template <
  class DomainType >
class ViewBlocksMod
: public ViewModBase< ViewBlocksMod<DomainType>, DomainType > {
 public:
  typedef DomainType                                           domain_type;
  typedef typename view_traits<DomainType>::origin_type        origin_type;
  typedef typename view_traits<DomainType>::index_type          index_type;
 private:
  typedef ViewBlocksMod<DomainType>                                 self_t;
  typedef ViewModBase<ViewBlocksMod<DomainType>, DomainType>        base_t;
  typedef ViewBlockMod<DomainType>                              block_type;
 public:
  typedef ddm::IndexSetBlocks<ViewBlocksMod<DomainType>>   index_set_type;
  typedef self_t                                               global_type;
  typedef typename domain_type::local_type                      local_type;

  typedef std::integral_constant<bool, false>                     is_local;

 private:
  index_set_type  _index_set;

 public:
  template <class ViewBlocksModType>
  class block_iterator
  : public internal::IndexIteratorBase<
             block_iterator<ViewBlocksModType>,
             block_type,
             index_type,
             std::nullptr_t,
             block_type > {
   private:
    typedef internal::IndexIteratorBase<
              block_iterator<ViewBlocksModType>,
              block_type,     // value type
              index_type,     // difference type
              std::nullptr_t, // pointer type
              block_type >    // reference type
      iterator_base_t;
   private:
//  ViewBlocksModType & _blocks_view;
    const self_t & _blocks_view;
   public:
    constexpr block_iterator()                         = delete;
    constexpr block_iterator(block_iterator &&)        = default;
    constexpr block_iterator(const block_iterator &)   = delete;
    ~block_iterator()                                  = default;
    block_iterator & operator=(block_iterator &&)      = default;
    block_iterator & operator=(const block_iterator &) = default;

    constexpr block_iterator(
      const block_iterator    & other,
      index_type                position)
    : iterator_base_t(position)
    , _blocks_view(other._blocks_view)
    { }

    constexpr block_iterator(
      const ViewBlocksModType & blocks_view,
      index_type          position)
    : iterator_base_t(position)
    , _blocks_view(blocks_view)
    { }

    constexpr block_type dereference(index_type idx) const {
      // Dereferencing block iterator returns block at block index
      // with iterator position.
      // Note that block index is relative to the domain and is
      // translated to global block index in IndexSetBlocks.
      return ViewBlockMod<DomainType>(
               ddm::domain(_blocks_view), idx);
    }
  };

 public:
  typedef block_iterator<self_t>                  iterator;
  typedef block_iterator<const self_t>      const_iterator;

 public:
  constexpr ViewBlocksMod()               = delete;
  constexpr ViewBlocksMod(const self_t &) = default;
  constexpr ViewBlocksMod(self_t &&)      = default;
  ~ViewBlocksMod()                        = default;
  self_t & operator=(self_t &&)           = default;
  self_t & operator=(const self_t &)      = default;

  /**
   * Constructor, creates a view on a given domain.
   */
  constexpr explicit ViewBlocksMod(
    const domain_type & domain)
  : base_t(domain)
  , _index_set(*this)
  { }

  /**
   * Constructor, creates a view on a given domain.
   */
  constexpr explicit ViewBlocksMod(
    domain_type && domain)
  : base_t(std::forward<domain_type>(domain))
  , _index_set(*this)
  { }

  constexpr iterator begin() const {
    return const_iterator(*this, _index_set.first());
  }
  inline iterator begin() {
    return iterator(*this, _index_set.first());
  }

  constexpr iterator end() const {
    return const_iterator(*this, _index_set.last() + 1);
  }
  inline iterator end() {
    return iterator(*this, _index_set.last() + 1);
  }

  constexpr block_type operator[](int offset) const {
    return *iterator(*this, _index_set[offset]);
  }
  inline block_type operator[](int offset) {
    return *iterator(*this, _index_set[offset]);
  }

  constexpr auto local() const
//constexpr const local_type local() const {
    -> decltype(ddm::local(
                  std::declval<
                    typename std::add_lvalue_reference<domain_type>::type
                  >() )) {
    return ddm::local(this->domain());
  }

  inline auto local()
//inline local_type local() {
    -> decltype(ddm::local(
                  std::declval<
                    typename std::add_lvalue_reference<domain_type>::type
                  >() )) {
    return ddm::local(this->domain());
  }

  constexpr const global_type global() const {
    return ddm::global(ddm::domain(*this));
  }

  inline global_type global() {
    return ddm::global(ddm::domain(*this));
  }

  constexpr const index_set_type & index_set() const {
    return _index_set;
  }
};

#endif // DOXYGEN

} // namespace ddm

#endif // DDM__VIEW__VIEW_BLOCKS_MOD_H__INCLUDED
