#ifndef DDM__IO__HDF5__IOSTREAM_h
#define DDM__IO__HDF5__IOSTREAM_h

#include "../../../ddm/LaunchPolicy.h"
#include "../../../ddm/io/IOStream.h"

namespace ddm {
namespace io {
namespace hdf5 {

// No subclassing necessary
using DeviceMode = ddm::io::IOSBaseMode;
using StreamMode = ddm::io::IOStreamMode<DeviceMode>;

}  // namespace hdf5
}  // namespace io
}  // namespace ddm

#include "../../../ddm/io/hdf5/IOManip.h"
#include "../../../ddm/io/hdf5/InputStream.h"
#include "../../../ddm/io/hdf5/OutputStream.h"

#endif  // DDM__IO__HDF5__IOSTREAM_h
