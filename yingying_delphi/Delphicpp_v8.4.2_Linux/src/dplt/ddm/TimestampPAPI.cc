#if defined(DDM_ENABLE_PAPI)

#include "../ddm/util/internal/TimestampPAPI.h"

namespace ddm {
namespace util {
namespace internal {

int TimestampPAPI<TimeMeasure::Clock>::timer_mode   = 0;
int TimestampPAPI<TimeMeasure::Counter>::timer_mode = 0;

Timestamp::counter_t
TimestampPAPI<TimeMeasure::Counter>::frequencyScaling = 1;

} // namespace internal
} // namespace util
} // namespace ddm

#endif // DDM_ENABLE_PAPI
