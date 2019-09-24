
#include "../ddm/util/Timer.h"
#include "../ddm/util/Timestamp.h"

#if defined(DDM__UTIL__TIMER_POSIX) || \
    defined(DDM__UTIL__TIMER_UX)

#include "../ddm/util/internal/TimestampCounterPosix.h"
#include "../ddm/internal/Logging.h"

namespace ddm {
namespace util {
namespace internal {

void TimestampCounterPosix::Calibrate(unsigned int freq) {
  frequencyScaling = freq == 0
                     ? 1900.0f
                     : static_cast<double>(freq);
  DDM_LOG_DEBUG("TimestampCounterPosix::Calibrate(freq)", freq);
  DDM_LOG_DEBUG("TimestampCounterPosix::Calibrate",
                 "timer:", TimerName());
  DDM_LOG_DEBUG("TimestampCounterPosix::Calibrate",
                 "fscale:", frequencyScaling);
}

Timestamp::counter_t TimestampCounterPosix::frequencyScaling = 1;

} // namespace internal
} // namespace util
} // namespace ddm

#endif // DDM__UTIL__TIMER_POSIX || DDM__UTIL__TIMER_UX
