#ifndef DDM__LIBDDM_H_
#define DDM__LIBDDM_H_

/**
 * The DDM C++ Library Interface
 *
 */
namespace ddm {

}

/**
 * \defgroup DDMConcept DDM C++ Concepts
 * Concepts for C++ components in DDM
 */

#include "internal/Config.h"

#include "Types.h"
#include "Init.h"
#include "Team.h"
#include "Cartesian.h"
#include "TeamSpec.h"

#include "Iterator.h"
#include "View.h"
#include "Range.h"

#include "GlobMem.h"
#include "GlobPtr.h"
#include "GlobRef.h"
#include "GlobAsyncRef.h"

#include "Onesided.h"

#include "LaunchPolicy.h"

#include "Container.h"
#include "Shared.h"
#include "SharedCounter.h"
#include "Exception.h"
#include "Algorithm.h"
#include "Allocator.h"
#include "Atomic.h"

#include "Pattern.h"

#include "util/BenchmarkParams.h"
#include "util/Config.h"
#include "util/Trace.h"
#include "util/PatternMetrics.h"
#include "util/Timer.h"

#include "util/Locality.h"
#include "util/LocalityDomain.h"
#include "util/TeamLocality.h"
#include "util/UnitLocality.h"
#include "util/LocalityJSONPrinter.h"

#include "IO.h"
#include "io/HDF5.h"

#include "internal/Math.h"
#include "internal/Logging.h"

#include "tools/PatternVisualizer.h"

#endif // DDM__LIBDDM_H_
