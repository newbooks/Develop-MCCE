#ifndef DDM__VIEW__UNION_H__INCLUDED
#define DDM__VIEW__UNION_H__INCLUDED

#include "../../ddm/Types.h"
#include "../../ddm/Dimensional.h"
#include "../../ddm/Cartesian.h"

#include <vector>


namespace ddm {


template <class ComponentViewType>
class CompositeView
{

public:
  CompositeView(std::initializer_list<ComponentViewType> views)
  : _views(views)
  { }

  CompositeView(const std::vector<ComponentViewType> & views)
  : _views(views)
  { }

private:
  std::vector<ComponentViewType> _views;
};


template <class ComponentViewType>
constexpr CompositeView<ComponentViewType>
set_union(
  const std::vector<ComponentViewType> & views) {
  return CompositeView<ComponentViewType>(views);
}

template <class ComponentViewType>
constexpr CompositeView<ComponentViewType>
set_union(
  std::initializer_list<ComponentViewType> views) {
  return CompositeView<ComponentViewType>(views);
}


} // namespace ddm

#endif // DDM__VIEW__UNION_H__INCLUDED
