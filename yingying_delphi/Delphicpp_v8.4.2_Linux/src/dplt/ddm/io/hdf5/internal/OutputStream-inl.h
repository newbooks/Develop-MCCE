#ifndef DDM__IO__INTERNAL__HDF5__OUTPUT_STREAM_INL_H__INCLUDED
#define DDM__IO__INTERNAL__HDF5__OUTPUT_STREAM_INL_H__INCLUDED

#include "../../../../ddm/io/hdf5/OutputStream.h"
#include "../../../../ddm/io/hdf5/StorageDriver.h"

namespace ddm {
namespace io {
namespace hdf5 {

template <typename Container_t>
inline OutputStream& operator<<(OutputStream& os, Container_t& container) {
  if (os._launch_policy == ddm::launch::async) {
    os._store_object_impl_async(container);
  } else {
    os._store_object_impl(container);
  }
  // Append future data in this stream
  os._foptions.overwrite_file = false;
  return os;
}

}  // namespace hdf5
}  // namespace io
}  // namespace ddm

#endif  // DDM__IO__INTERNAL__HDF5__OUTPUT_STREAM_INL_H__INCLUDED
