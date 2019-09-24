#ifndef DDM__RANGES_H__INCLUDED
#define DDM__RANGES_H__INCLUDED

/**
 * \defgroup  DDMRangeConcept  Multidimensional Range Concept
 *
 * \ingroup DDMNDimConcepts
 * \{
 * \par Description
 *
 * Definitions for multidimensional range expressions.
 *
 * \see DDMDimensionalConcept
 * \see DDMViewConcept
 * \see DDMIteratorConcept
 * \see \c ddm::view_traits
 *
 * Variables used in the following:
 *
 * - \c r instance of a range model type
 * - \c o index type, representing element offsets in the range and their
 *        distance
 * - \c i iterator referencing elements in the range
 *
 * \par Types
 *
 * \par Expressions
 *
 * Expression               | Returns | Effect | Precondition | Postcondition
 * ------------------------ | ------- | ------ | ------------ | -------------
 * <tt>*ddm::begin(r)</tt> |         |        |              | 
 * <tt>r[o]</tt>            |         |        |              | 
 *
 * \par Functions
 *
 * - \c ddm::begin
 * - \c ddm::end
 * - \c ddm::distance
 * - \c ddm::size
 *
 * \par Metafunctions
 *
 * - \c ddm::is_range<X>
 *
 * \}
 */


#include "../ddm/Types.h"

#include <type_traits>


namespace ddm {


#ifndef DOXYGEN


// Related: boost::range
//
// https://github.com/boostorg/range/tree/develop/include/boost/range
//

template <typename ViewT>
struct view_traits;

// Forward-declaration
template <typename ViewType>
class IndexSetIdentity;

// Forward-declaration
template <typename Iterator, typename Sentinel = Iterator>
class IteratorRange;

// Forward-declaration
template <typename Iterator, typename Sentinel>
class IteratorRange<Iterator *, Sentinel *>;

#endif



/**
 * \concept{DDMRangeConcept}
 */
template <typename RangeType>
constexpr auto begin(const RangeType & range)
  -> decltype(range.begin()) {
  return range.begin();
}

/**
 * \concept{DDMRangeConcept}
 */
template <class RangeType>
constexpr auto end(const RangeType & range)
  -> decltype(range.end()) {
  return range.end();
}

/**
 * \concept{DDMRangeConcept}
 */
template <class RangeType>
constexpr auto
size(const RangeType & r)
  -> decltype(r.size()) {
  return r.size();
}


namespace detail {

template<typename T>
struct _is_range_type
{
private:
  typedef char yes;
  typedef long no;

  typedef typename std::remove_reference<
            typename std::remove_const<T>::type
          >::type
    ValueT;

#ifdef __TODO__
private:
  // Test if ddm::begin(x) is valid expression:
  template <typename C> static yes has_ddm_begin(
                                     decltype(
                                       ddm::begin(
                                         std::move(std::declval<T>())
                                       )
                                     ) * );
  template <typename C> static no  has_ddm_begin(...);    

  // Test if ddm::end(x) is valid expression:
  template <typename C> static yes has_ddm_end(
                                     decltype(
                                       ddm::end(
                                         std::move(std::declval<T>())
                                       )
                                     ) * );
  template <typename C> static no  has_ddm_end(...);    

public:
  enum { value = (
          sizeof(has_ddm_begin(static_cast<ValueT*>(nullptr))) == sizeof(yes)
       && sizeof(has_ddm_end(static_cast<ValueT*>(nullptr)))   == sizeof(yes)
         ) };

  //template<typename C, typename begin_decl =
  //                                   decltype(
  //                                     ddm::begin(
  //                                       std::move(std::declval<T>())
  //                                     )) >
  //static yes has_ddm_begin(C *);
#endif
  // Test if x.begin() is valid expression and type x::iterator is
  // defined:
  template<typename C, typename C::iterator (C::*)() const = &C::begin >
  static yes has_begin(C *);
  template<typename C, typename C::iterator (C::*)() = &C::begin >
  static yes has_begin(C *);
  static no  has_begin(...);

  // Test if x.end() is valid expression and type x::iterator is
  // defined:
  template<typename C, typename C::iterator (C::*)() const = &C::end >
  static yes has_end(C *);
  template<typename C, typename C::iterator (C::*)() = &C::end >
  static yes has_end(C *);
  static no  has_end(...);

public:
  enum { value = (
              sizeof(has_begin(static_cast<ValueT*>(nullptr))) == sizeof(yes)
           && sizeof(has_end(static_cast<ValueT*>(nullptr)))   == sizeof(yes)
         ) };
};

} // namespace detail

/**
 * Type trait for testing if `ddm::begin<T>` and `ddm::end<T>`
 * are defined.
 *
 * In the current implementation, range types must specify the
 * return type of `ddm::begin<T>` and `ddm::end<T>` as type
 * definition `iterator`.
 *
 * This requirement will become obsolete in the near future.
 *
 *
 * Example:
 *
 * \code
 *   bool g_array_is_range = ddm::is_range<
 *                                    ddm::Array<int>
 *                                 >::value;
 *   // -> true
 *
 *   bool l_array_is_range = ddm::is_range<
 *                                    typename ddm::Array<int>::local_type
 *                                 >::value;
 *   // -> true
 *
 *   struct inf_range { 
 *     typedef int           * iterator;
 *     typedef std::nullptr_t  sentinel;
 *
 *     iterator begin() { ... }
 *     sentinel end()   { ... }
 *   };
 *
 *   bool inf_range_is_range = ddm::is_range<inf_range>::value;
 *   // -> false
 *   //    because of missing definition
 *   //      iterator ddm::end<inf_range> -> iterator
 *
 *   Currently requires specialization as workaround:
 *   template <>
 *   struct is_range<inf_range> : std::integral_value<bool, true> { };
 * \endcode
 */
template <class RangeType>
struct is_range : ddm::detail::_is_range_type<RangeType> { };

template <
  typename RangeType,
  typename Iterator,
  typename Sentinel = Iterator >
class RangeBase {
public:
  typedef Iterator                iterator;
  typedef Sentinel                sentinel;
  typedef ddm::default_index_t   index_type;
protected:
  RangeType & derived() {
    return static_cast<RangeType &>(*this);
  }
  const RangeType & derived() const {
    return static_cast<const RangeType &>(*this);
  }
};


/**
 * Specialization of \c ddm::view_traits for IteratorRange.
 */
template <
  typename IteratorT,
  typename SentinelT >
struct view_traits<ddm::IteratorRange<IteratorT, SentinelT> > {
private:
  typedef IteratorRange<IteratorT, SentinelT> RangeT;
public:
  typedef RangeT                                               origin_type;
  typedef RangeT                                               domain_type;
  typedef RangeT                                                image_type;
  typedef RangeT                                               global_type;
  typedef typename RangeT::local_type                           local_type;
  typedef typename RangeT::index_type                           index_type;
  typedef typename RangeT::index_set_type                   index_set_type;

  /// Whether the view type is a projection (has less dimensions than the
  /// view's domain type).
  typedef std::integral_constant<bool, false>                is_projection;
  typedef std::integral_constant<bool, false>                is_view;
  /// Whether the view is the origin domain.
  typedef std::integral_constant<bool, true>                 is_origin;
  /// Whether the view / container type is a local view.
  /// \note A container type is local if it is identical to its
  ///       \c local_type
  typedef std::integral_constant<bool, std::is_same<
                                 RangeT,
                                 typename RangeT::local_type
                                >::value >                   is_local;
};

/**
 * Adapter template for range concept, wraps `begin` and `end` iterators
 * in range type.
 */
template <
  typename Iterator,
  typename Sentinel >
class IteratorRange
: public RangeBase< IteratorRange<Iterator, Sentinel>,
                    Iterator,
                    Sentinel >
{
  typedef IteratorRange<Iterator, Sentinel> self_t;

  Iterator & _begin;
  Sentinel & _end;

public:
  typedef Iterator                                              iterator;
  typedef Sentinel                                              sentinel;
  typedef ddm::default_index_t                               index_type;
  typedef typename iterator::pattern_type                   pattern_type;
  typedef ddm::IndexSetIdentity<self_t>                  index_set_type;
  typedef typename iterator::value_type                       value_type;

  typedef typename
            std::conditional<
              std::is_pointer<iterator>::value,
              iterator,
              typename iterator::local_type
            >::type
    local_iterator;

  typedef typename
            std::conditional<
              std::is_pointer<sentinel>::value,
              iterator,
              typename sentinel::local_type
            >::type
    local_sentinel;
            
  typedef IteratorRange<local_iterator, local_sentinel>       local_type;


public:
  template <class Container>
  constexpr explicit IteratorRange(Container && c)
  : _begin(c.begin())
  , _end(c.end())
  { }

  constexpr IteratorRange(iterator & begin, sentinel & end)
  : _begin(begin)
  , _end(end)
  { }

  constexpr iterator begin() const { return _begin; }
  constexpr iterator end()   const { return _end;   }

  constexpr const local_type local() const {
    return local_type(
     //      ddm::local(_begin),
     //      ddm::local(_end)
             _begin.local(),
             _end.local()
           );
  }

  constexpr const pattern_type & pattern() const {
    return _begin.pattern();
  }

  constexpr index_set_type index_set() const {
    return index_set_type(*this);
  }
};


/**
 * Adapter template for range concept, wraps `begin` and `end` iterators
 * in range type.
 */
template <
  typename LocalIterator,
  typename LocalSentinel >
class IteratorRange<LocalIterator *, LocalSentinel *>
: public RangeBase<
           IteratorRange<LocalIterator *, LocalSentinel *>,
           LocalIterator *,
           LocalSentinel * >
{
  typedef IteratorRange<LocalIterator, LocalSentinel> self_t;

  LocalIterator * _begin;
  LocalSentinel * _end;

public:
  typedef LocalIterator *                                       iterator;
  typedef LocalSentinel *                                       sentinel;
  typedef ddm::default_index_t                               index_type;
  typedef ddm::IndexSetIdentity<self_t>                  index_set_type;
  typedef LocalIterator                                       value_type;

  typedef iterator local_iterator;
  typedef sentinel local_sentinel;
            
  typedef IteratorRange<local_iterator, local_sentinel>       local_type;

public:
  template <class Container>
  constexpr explicit IteratorRange(Container && c)
  : _begin(c.begin())
  , _end(c.end())
  { }

  constexpr IteratorRange(iterator & begin, sentinel & end)
  : _begin(begin)
  , _end(end)
  { }

  constexpr iterator begin() const { return _begin; }
  constexpr iterator end()   const { return _end;   }

  constexpr const local_type & local() const {
    return *this;
  }

  constexpr index_set_type index_set() const {
    return index_set_type(*this);
  }
};

/**
 * Adapter utility function.
 * Wraps `begin` and `end` iterators in range type.
 */
template <class Iterator, class Sentinel>
constexpr ddm::IteratorRange<const Iterator, const Sentinel>
make_range(
  const Iterator & begin,
  const Sentinel & end) {
  return ddm::IteratorRange<const Iterator, const Sentinel>(
           begin,
           end);
}

template <class Iterator, class Sentinel>
constexpr ddm::IteratorRange<Iterator *, Sentinel *>
make_range(
  Iterator * begin,
  Sentinel * end) {
  return ddm::IteratorRange<Iterator *, Sentinel *>(
           begin,
           end);
}

template <class Iterator, class Sentinel>
ddm::IteratorRange<Iterator, Sentinel>
make_range(
  Iterator & begin,
  Sentinel & end) {
  return ddm::IteratorRange<Iterator, Sentinel>(
           begin,
           end);
}

} // namespace ddm

#endif // DDM__RANGES_H__INCLUDED
