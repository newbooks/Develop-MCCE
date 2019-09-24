#ifndef DDM__IO__HDF5__OUTPUT_STREAM_H__
#define DDM__IO__HDF5__OUTPUT_STREAM_H__

#ifdef DDM_ENABLE_HDF5

#include <string>
#include <future>

#include "../../../ddm/Matrix.h"
#include "../../../ddm/Array.h"

#include "../../../ddm/LaunchPolicy.h"

#include <chrono>
#include <thread>

namespace ddm {
namespace io {
namespace hdf5 {

/**
 * DDM stream API to store a ddm container or view
 * in an HDF5 file using parallel IO.
 *
 * All operations are collective.
 */

class OutputStream : public ::ddm::io::IOSBase<StreamMode> {

  typedef OutputStream self_t;
  typedef ddm::io::IOSBase<StreamMode> base_t;
  typedef StreamMode mode_t;

 private:
  std::string _filename;
  std::string _dataset;
  type_converter _converter;
  hdf5_options _foptions = hdf5_options();
  bool _use_cust_conv = false;
  ddm::launch _launch_policy;

  std::vector<std::shared_future<void> > _async_ops;

 public:
  /**
   * Creates an HDF5 output stream using a launch policy
   *
   * Support of \cddm::launch::async is still highly experimental and requires
   * thread support in MPI. If multi-threaded access is not supported,
   * blocking I/O is used as fallback. To wait for outstanding IO operations use
   * \cflush(). Until the stream is not flushed, no write accesses to the
   * container, as well as no barriers are allowed.
   * Otherwise the behavior is undefined.
   */
  OutputStream(
      ///
      ddm::launch lpolicy, std::string filename,
      /// device opening flags: \see ddm::io::IOSBaseMode
      mode_t open_mode = DeviceMode::no_flags)
      : _filename(filename), _dataset("data"), _launch_policy(lpolicy) {
    if ((open_mode & DeviceMode::app)) {
      _foptions.overwrite_file = false;
    }
    if (_launch_policy == ddm::launch::async && !ddm::is_multithreaded()) {
      _launch_policy = ddm::launch::sync;
      DDM_LOG_WARN(
          "Requested ASIO but DART does not support "
          "multi-threaded access. Blocking IO is used"
          "as fallback");
    }
  }

  /**
   * Creates an HDF5 output stream using blocking IO.
   *
   * The stream takes an arbitrary number of modifiers and objects,
   * where the objects are stored in the order of passing it to the stream.
   *
   * The interface follows roughly the STL stream concept.
   *
   * Example:
   * \code
   *  ddm::Array<double> array_a(1000);
   *  ddm::Array<double> array_b(500);
   *
   *  OutputStream os("file.hdf5");
   *  os << dataset("dataset") << array_a
   *     << dataset("seconddata") << array_b;
   * \endcode
   */
  OutputStream(std::string filename, mode_t open_mode = DeviceMode::no_flags)
      : OutputStream(ddm::launch::sync, filename, open_mode) {}

  ~OutputStream() {
    if (!_async_ops.empty()) {
      _async_ops.back().wait();
    }
  }

  /**
   * Synchronizes with the data sink.
   * If \cddm::launch::async is used, waits until all data is written
   */
  OutputStream flush() {
    DDM_LOG_DEBUG("flush output stream", _async_ops.size());
    if (!_async_ops.empty()) {
      _async_ops.back().wait();
    }
    DDM_LOG_DEBUG("output stream flushed");
    return *this;
  }

  // IO Manipulators

  /// set name of dataset
  friend OutputStream& operator<<(OutputStream& os, const dataset tbl) {
    os._dataset = tbl._dataset;
    return os;
  }

  /// set metadata key at which the pattern will be stored
  friend OutputStream& operator<<(OutputStream& os, const setpattern_key pk) {
    os._foptions.pattern_metadata_key = pk._key;
    return os;
  }

  /// set whether pattern layout should be stored as metadata
  friend OutputStream& operator<<(OutputStream& os, const store_pattern sp) {
    os._foptions.store_pattern = sp._store;
    return os;
  }

  /// modify an already existing dataset instead of overwriting it
  friend OutputStream& operator<<(OutputStream& os, const modify_dataset md) {
    os._foptions.modify_dataset = md._modify;
    return os;
  }

  /// custom type converter function to convert native type to HDF5 type
  friend OutputStream& operator<<(OutputStream& os, const type_converter conv) {
    os._converter = conv;
    os._use_cust_conv = true;
    return os;
  }

  /// kicker which stores an container using the specified stream properties.
  template <typename Container_t>
  friend OutputStream& operator<<(OutputStream& os, Container_t& container);

 private:
  template <typename Container_t>
  void _store_object_impl(Container_t& container) {
    if (_use_cust_conv) {
      StoreHDF::write(container, _filename, _dataset, _foptions, _converter);
    } else {
      StoreHDF::write(container, _filename, _dataset, _foptions);
    }
  }

  template <typename Container_t>
  void _store_object_impl_async(Container_t& container) {
    auto pos = _async_ops.size();

    // copy state of stream
    auto s_filename = _filename;
    auto s_dataset = _dataset;
    auto s_foptions = _foptions;
    auto s_use_cust_conv = _use_cust_conv;
    type_converter_fun_type s_converter = _converter;

    // pass pos by value as it might be out of scope when function is called
    std::shared_future<void> fut = std::async(
        std::launch::async, [&, pos, s_filename, s_dataset, s_foptions,
                             s_converter, s_use_cust_conv]() {
          if (pos != 0) {
            // wait for previous tasks
            auto last_task = _async_ops[pos - 1];
            DDM_LOG_DEBUG("waiting for future", pos);
            last_task.wait();
          }
          DDM_LOG_DEBUG("execute async io task");

          if (s_use_cust_conv) {
            StoreHDF::write(container, s_filename, s_dataset, s_foptions,
                            s_converter);
          } else {
            StoreHDF::write(container, s_filename, s_dataset, s_foptions);
          }
          DDM_LOG_DEBUG("execute async io task done");
        });
    _async_ops.push_back(fut);
  }
};

}  // namespace hdf5
}  // namespace io
}  // namespace ddm

#include "../../../ddm/io/hdf5/internal/OutputStream-inl.h"

#endif  // DDM_ENABLE_HDF5

#endif  // DDM__IO__HDF5__OUTPUT_STREAM_H__
