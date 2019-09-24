#ifndef DDM__VIEW__CART_VIEW_H__INCLUDED
#define DDM__VIEW__CART_VIEW_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Dimensional.h"
#include "../../ddm/Cartesian.h"

#include <iterator>


namespace ddm {

/**
 * Base class for a cartesian view, i.e. an n-dimensional view with 
 * cartesian coordinates.
 */
template<
  typename   Iter,
  dim_t      NumDimensions,
  MemArrange Arrangement = ROW_MAJOR,
  typename SizeType      = ddm::default_size_t >
class CartViewBase { 
public: 
  typedef typename std::iterator_traits<Iter>::value_type value_type;
  typedef typename std::iterator_traits<Iter>::reference  reference;
  
private:
  CartesianIndexSpace<NumDimensions, Arrangement, SizeType> m_cart;
  Iter                                            m_begin;
  
public:
  // construct from iterator
  template<typename ... Args>
  CartViewBase(Iter it, Args... args) : 
    m_cart { (SizeType)args ... },
    m_begin { it } {  
  }

  // construct from container
  template<typename Container, typename... Args>
  CartViewBase(Container & container, Args... args) : 
    m_cart { args ... },
    m_begin { container.begin() } {  
  }

  constexpr SizeType rank() const {
    return m_cart.rank();
  }

  constexpr SizeType size() const {
    return m_cart.size();
  }

  constexpr SizeType extent(dim_t dim) const {
    return m_cart.extent(dim);
  }

  template<typename ... Args>
  reference at(Args ... args) const {
    Iter it = m_begin;
    std::advance(it, m_cart.at(args...));
    return *it;
  }

  // x(), y(), z() accessors 
  // enabled only for the appropriate sizes
  template<dim_t U=NumDimensions>
  constexpr typename std::enable_if<(U>0),SizeType>::type
  x(SizeType offs) const {
    return m_cart.x(offs);
  }
  
  template<dim_t U=NumDimensions>
  constexpr typename std::enable_if<(U>1),SizeType>::type
  y(SizeType offs) const {
    return m_cart.y(offs);
  }
  
  template<dim_t U=NumDimensions>
  constexpr typename std::enable_if<(U>2),SizeType>::type
  z(SizeType offs) const {
    return m_cart.z(offs);
  }

};

/**
 * Cartesian view class.
 */
template<
  typename   Iter,
  dim_t      NumDimensions,
  MemArrange Arrangement = ROW_MAJOR,
  typename SizeType      = ddm::default_size_t >
struct CartView 
: public CartViewBase<Iter, NumDimensions, Arrangement, SizeType> {
public:
  template<typename... Args>
  CartView(
    Iter it,
    Args... args) 
  : CartViewBase<Iter, NumDimensions, Arrangement, SizeType>(
      it,
      args...) { }

  template<typename Container, typename... Args>
  CartView(
    Container & cont,
    Args... args)
  : CartViewBase<Iter, NumDimensions, Arrangement, SizeType>(
      cont,
      args...) { }
};

} // namespace ddm

#endif // DDM__VIEW__CART_VIEW_H__INCLUDED
