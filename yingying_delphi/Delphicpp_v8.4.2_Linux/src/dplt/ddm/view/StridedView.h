#ifndef DDM__VIEW__STRIDED_VIEW_H__INCLUDED
#define DDM__VIEW__STRIDED_VIEW_H__INCLUDED

#include "../../ddm/Types.h"

#include "../../ddm/view/SetUnion.h"
#include "../../ddm/view/MultiView.h"

#include <vector>


namespace ddm {

template <dim_t NDim>
class StridedView;

template <>
class StridedView<0>;


template <dim_t NDim>
class StridedView
: public ddm::CompositeView<
           ddm::MultiView<NDim-1>
         >
{

};

} // namespace ddm

#endif // DDM__VIEW__STRIDED_VIEW_H__INCLUDED
