#ifndef DDM__IO__HDF5__INTERNAL__INPUT_STREAM_INL_H__INCLUDED
#define DDM__IO__HDF5__INTERNAL__INPUT_STREAM_INL_H__INCLUDED

#include "../../../../ddm/io/hdf5/InputStream.h"
#include "../../../../ddm/io/hdf5/StorageDriver.h"

namespace ddm {
namespace io {
namespace hdf5 {

template <typename Container_t>
inline InputStream& operator>>(InputStream& is, Container_t& container) {
  if (is._launch_policy == ddm::launch::async) {
    is._load_object_impl_async(container);
  } else {
    is._load_object_impl(container);
  }
  return is;
}

}  // namespace hdf5
}  // namespace io
}  // namespace ddm

#endif  // DDM__IO__HDF5__INTERNAL__INPUT_STREAM_INL_H__INCLUDED
