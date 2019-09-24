#include "../ddm/util/Timer.h"

#if defined(DDM__UTIL__TIMER_POSIX) || \
    defined(DDM__UTIL__TIMER_OSX) || \
    defined(DDM__UTIL__TIMER_UX)

#include "../ddm/util/internal/TimestampClockPosix.h"
#include "../ddm/internal/Logging.h"

#include <stdexcept>

// POSIX standard header
#include <unistd.h>   // POSIX flags
#include <time.h>     // clock_gettime(), time(), needs librt.a
#include <sys/time.h> // gethrtime(), gettimeofday() fallbacks

#ifndef CLOCK_MONOTONIC_RAW
// Fallback
#  define CLOCK_MONOTONIC_RAW 4
#endif

#if defined(__MACH__) && defined(__APPLE__)
// OS X / Mach
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

namespace ddm {
namespace util {
namespace internal {

TimestampClockPosix::TimestampClockPosix()
{
#if defined(DDM__UTIL__TIMER_UX)
// HP-UX, Solaris
  value = static_cast<Timestamp::counter_t>(gethrtime());

#elif defined(DDM__UTIL__TIMER_OSX)
// OS X
  static double timeConvert = 0.0;
  if (timeConvert == 0.0)
  {
    mach_timebase_info_data_t timeBase;
    (void)mach_timebase_info(&timeBase);
    timeConvert = (double)timeBase.numer /
                  (double)timeBase.denom;
  }
  value = static_cast<Timestamp::counter_t>(
            static_cast<double>(mach_absolute_time() * timeConvert));

#elif defined(DDM__UTIL__TIMER_POSIX)
// POSIX
#if defined(_POSIX_TIMERS) && (_POSIX_TIMERS > 0)
// POSIX clock_gettime
  struct timespec ts;
  if (clockId != (clockid_t)-1 && clock_gettime(clockId, &ts) != -1) {
    value = static_cast<Timestamp::counter_t>(
              static_cast<double>(ts.tv_sec * 1000000) +
              static_cast<double>(ts.tv_nsec / 1000));
    return;
  }

#endif
// No support for POSIX clock_gettime
  // Imprecise fallback to wall-clock time
  struct timeval tm;
  gettimeofday(&tm, NULL);
  value = static_cast<Timestamp::counter_t>(
            static_cast<double>(tm.tv_sec * 1000000) +
            static_cast<double>(tm.tv_usec));
#else
// No POSIX-compliant time mechanism found.
  throw std::runtime_error("Could not resolve timer");
#endif
}

const char *
TimestampClockPosix::clockModeNames[] =
{
  "Default",
  "GENERIC",
  "UX",
  "MACH",
  "THREAD_CPUTIME",
  "PROC_CPUTIME",
  "HIGHRES",
  "MONOTONIC_PRECISE",
  "MONOTONIC_RAW",
  "MONOTONIC",
  "REALTIME"
};

TimestampClockPosix::ClockDef
TimestampClockPosix::availableModes[] =
{
  TimestampClockPosix::ClockDef(
      UNDEFINED_CLOCK,   static_cast<clockid_t>(-1)),
#if defined(DDM__UTIL__TIMER_UX)
  TimestampClockPosix::ClockDef(
      UX,                static_cast<clockid_t>(-1)),
#endif
#if defined(DDM__UTIL__TIMER_OSX)
  TimestampClockPosix::ClockDef(
      MACH,              static_cast<clockid_t>(-1)),
#endif
#if defined(_CLOCK_HIGHRES)
  TimestampClockPosix::ClockDef(
      HIGHRES,           CLOCK_HIGHRES),
#endif
#if defined(_POSIX_MONOTONIC_CLOCK)
# if defined(DDM__UTIL__TIMER_FREEBSD)
  TimestampClockPosix::ClockDef(
      MONOTONIC_PRECISE, CLOCK_MONOTONIC_PRECISE),
  ///< FreeBSD specific, correcponds to MONOTONIC on Linux
# endif // FreeBSD
# if defined(DDM__UTIL__TIMER_LINUX)
  TimestampClockPosix::ClockDef(
      MONOTONIC_RAW,     CLOCK_MONOTONIC_RAW),
  ///< Linux specific, like MONOTONIC but without NTP adjustment
# endif // Linux
  TimestampClockPosix::ClockDef(
      MONOTONIC,         CLOCK_MONOTONIC),
#endif // _POSIX_MONOTONIC_CLOCK
#if defined(_POSIX_THREAD_CPUTIME)
  TimestampClockPosix::ClockDef(
      THREAD_CPUTIME,    CLOCK_THREAD_CPUTIME_ID),
#endif
#if defined(_POSIX_CPUTIME)
  TimestampClockPosix::ClockDef(
      PROC_CPUTIME,      CLOCK_PROCESS_CPUTIME_ID),
#endif
  TimestampClockPosix::ClockDef(
      REALTIME,          CLOCK_REALTIME),
  TimestampClockPosix::ClockDef(
      GENERIC_CLOCK,     static_cast<clockid_t>(-1))
};

void TimestampClockPosix::Calibrate(
  unsigned int mode)
{
  const unsigned int lastAvMode = static_cast<unsigned int>(
    TimestampClockPosix::GENERIC_CLOCK);
  // Default to second index in available modes, which
  // is the first and preferred clock type
  unsigned int selectedModeIndex  = 1;
  // Iterate over all available clock types:
  DDM_LOG_DEBUG("TimestampClockPosix::Calibrate(mode)", mode);
  DDM_LOG_TRACE("TimestampClockPosix::Calibrate", "Available modes:");
  for (unsigned int avModeIdx = 1;
      TimestampClockPosix::availableModes[avModeIdx].first != lastAvMode;
       ++avModeIdx) {
    unsigned int modeNum =
      TimestampClockPosix::availableModes[avModeIdx].first;
    if (modeNum == mode) {
      // Selected mode id is contained in available modes
      selectedModeIndex = avModeIdx;
    }
    DDM_LOG_TRACE("TimestampClockPosix::Calibrate",
                   "mode:",    TimestampClockPosix::clockModeNames[modeNum],
                   "mode id:", modeNum);
  }

  clockMode = TimestampClockPosix::availableModes[selectedModeIndex].first;
  clockId   = TimestampClockPosix::availableModes[selectedModeIndex].second;
  // Print mode that finally has been activated:
#ifdef DDM_ENABLE_LOGGING
  unsigned int modeNum = static_cast<unsigned int>(clockMode);
  DDM_LOG_DEBUG("TimestampClockPosix::Calibrate",
                 "active mode:", TimestampClockPosix::clockModeNames[modeNum],
                 "mode id:",     modeNum);
#endif
  // Print resolution of the active clock:
  struct timespec res;
  if (clock_getres(clockId, &res) == 0) {
    DDM_LOG_DEBUG("TimestampClockPosix::Calibrate",
                   "resolution:", res.tv_nsec);
  }
}

TimestampClockPosix::ClockMode
TimestampClockPosix::clockMode = TimestampClockPosix::UNDEFINED_CLOCK;

clockid_t
TimestampClockPosix::clockId = static_cast<clockid_t>(-1);

Timestamp::counter_t
TimestampClockPosix::frequencyScaling = 1;

} // namespace internal
} // namespace util
} // namespace ddm

#endif // POSIX || OSX || UX
