#ifndef DDM__UTIL__TIMESTAMP_H_
#define DDM__UTIL__TIMESTAMP_H_

#include "../../ddm/internal/Config.h"
#include <limits.h>

// OS X
#if defined(DDM__PLATFORM__OSX)
#define DDM__UTIL__TIMER_OSX
#endif
// HPUX / Sun
#if defined(DDM__PLATFORM__UX)
#define DDM__UTIL__TIMER_UX
#endif
// POSIX
#if defined(DDM__PLATFORM__POSIX)
#define DDM__UTIL__TIMER_POSIX
#endif
// Linux
#if defined(__linux__)
#define DDM__UTIL__TIMER_LINUX
#endif
// FreeBSD
#if defined(__FreeBSD__)
#define DDM__UTIL__TIMER_FREEBSD
#endif

#if defined(DDM_ENABLE_PAPI)
#define DDM__UTIL__TIMER_PAPI
#endif

namespace ddm {
namespace util {

class Timestamp {
 public:
  typedef unsigned long long counter_t;

 public:
  virtual ~Timestamp() { };
  virtual const counter_t & Value() const = 0;

  static double FrequencyScaling();
  static double FrequencyPrescale();
  static const char * VariantName();
  inline static counter_t TimestampInfinity() {
    return LLONG_MAX;
  }
  inline static counter_t TimestampNegInfinity() {
    return 0;
  }
};

} // namespace util
} // namespace ddm

#endif // DDM__UTIL__TIMESTAMP_H_
