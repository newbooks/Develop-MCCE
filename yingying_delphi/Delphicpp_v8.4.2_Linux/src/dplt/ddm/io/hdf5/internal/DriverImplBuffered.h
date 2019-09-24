#ifndef DDM__IO__HDF5__INTERNAL_IMPL_BUFFERED_H__
#define DDM__IO__HDF5__INTERNAL_IMPL_BUFFERED_H__

#include <hdf5.h>
#include <hdf5_hl.h>

namespace ddm {
namespace io {
namespace hdf5 {

template <class Container_t>
void StoreHDF::_write_dataset_impl_buffered(Container_t& container,
                                            const hid_t& h5dset,
                                            const hid_t& internal_type) {
  // TODO
}

#endif  // DDM__IO__HDF5__INTERNAL_IMPL_BUFFERED_H__

}  // namespace hdf5
}  // namespace io
}  // namespace ddm