#ifndef DDM__IO__HDF5__INTERNAL__STORAGEDRIVER_H__INCLUDED
#define DDM__IO__HDF5__INTERNAL__STORAGEDRIVER_H__INCLUDED

namespace ddm {
namespace io {
namespace hdf5 {

/// Datatype is not set yet
struct UNSPECIFIED_TYPE;

/** Pseudo type traits to map the native c datatype to an hdf5 type.
 * This has to be implemented using a function, as the H5T_NATIVE_*
 * is a macro that expands to a non constant function
 */
template <typename T>
inline hid_t get_h5_datatype() {
  DDM_THROW(ddm::exception::InvalidArgument, "Datatype not supported");
  // To avoid compiler warning:
  return -1;
}

template <>
inline hid_t get_h5_datatype<int>() {
  return H5T_NATIVE_INT;
}

template <>
inline hid_t get_h5_datatype<long>() {
  return H5T_NATIVE_LONG;
}

template <>
inline hid_t get_h5_datatype<float>() {
  return H5T_NATIVE_FLOAT;
}

template <>
inline hid_t get_h5_datatype<double>() {
  return H5T_NATIVE_DOUBLE;
}
}
}
}
#endif  // DDM__IO__HDF5__INTERNAL__STORAGEDRIVER_H__INCLUDED
