#ifndef DDM__UTIL__INTERNAL__TIMESTAMP_CLOCK_POSIX_H_
#define DDM__UTIL__INTERNAL__TIMESTAMP_CLOCK_POSIX_H_

#include "../../../ddm/util/Timer.h"
#include "../../../ddm/util/Timestamp.h"

#include <climits>

#include <unistd.h>   // POSIX flags
#include <time.h>     // clock_gettime(), time()
#include <sys/time.h> // gethrtime(), gettimeofday() fallbacks

#include <utility>
#include <iostream>
#include <stdexcept>

// OSX
#if defined(__MACH__) && defined(__APPLE__)
typedef clock_t clockid_t;
#endif



namespace ddm {
namespace util {
namespace internal {

/**
 * Timestamp counter (TSC) for POSIX platforms.
 */
class TimestampClockPosix : public Timestamp
{
 private:
  Timestamp::counter_t value;

 public:
  /**
   * Clock modes, ordered by preference.
   */
  typedef enum {
    UNDEFINED_CLOCK = 0, ///< Undefined (unresolved) mode
    GENERIC_CLOCK,       ///< Unknown (unresolvable) mode

    UX,                  ///< HP/UX mode (gethrtime)
    MACH,                ///< MACH/OS X mode

    THREAD_CPUTIME,      ///< POSIX clock mode CLOCK_THREAD_CPUTIME_ID
    PROC_CPUTIME,        ///< POSIX clock mode CLOCK_PROCESS_CPUTIME_ID
    HIGHRES,             ///< POSIX clock mode CLOCK_HIGHRES
    MONOTONIC_PRECISE,   ///< POSIX clock mode CLOCK_MONOTONIC_PRECISE
    MONOTONIC_RAW,       ///< POSIX clock mode CLOCK_MONOTONIC_RAW
    MONOTONIC,           ///< POSIX clock mode CLOCK_MONOTONIC
    REALTIME,            ///< POSIX clock mode CLOCK_REALTIME

    NUM_CLOCK_MODES      ///< Number of clock modes in total
  } ClockMode;

 private:
  typedef std::pair<ClockMode, clockid_t> ClockDef;

 private:
  /**
   * List of clock modes as readable strings.
   */
  static const char * clockModeNames[];
  /**
   * Mapping of clock modes to clockid types, only contains mappings for
   * clock types available on the platform.
   */
  static TimestampClockPosix::ClockDef availableModes[];
  /**
   * Currently selected clock mode.
   */
  static ClockMode clockMode;
  /**
   * Currently selected clockid_t corresponding to clockMode.
   */
  static clockid_t clockId;

 public:
  static Timestamp::counter_t frequencyScaling;

 public:
  static void Calibrate(unsigned int mode = 0);

  inline TimestampClockPosix(
    const Timestamp::counter_t & counterValue)
  : value(counterValue)
  { }

  TimestampClockPosix();

  inline TimestampClockPosix(
    const TimestampClockPosix & other)
  : value(other.value)
  { }

  inline TimestampClockPosix & operator=(
    const TimestampClockPosix & rhs) {
    if (this != &rhs) {
      value = rhs.value;
    }
    return *this;
  }

  inline const Timestamp::counter_t & Value() const {
    return value;
  }

  inline static double FrequencyScaling() {
    return 1.0f;
  }

  inline static double FrequencyPrescale() {
    return 1.0f;
  }

  inline static const char * TimerName() {
    return clockModeNames[clockMode];
  }

};

} // namespace internal
} // namespace util
} // namespace ddm

#endif // DDM__UTIL__INTERNAL__TIMESTAMP_CLOCK_POSIX_H_
