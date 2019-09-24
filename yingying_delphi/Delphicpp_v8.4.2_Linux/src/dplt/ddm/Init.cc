
#include "../ddm/Init.h"
#include "../ddm/Team.h"
#include "../ddm/Types.h"
#include "../ddm/Shared.h"

#include "../ddm/util/Locality.h"
#include "../ddm/util/Config.h"
#include "../ddm/internal/Logging.h"

#include "../ddm/internal/Annotation.h"


namespace ddm {
  static bool _initialized   = false;
  static bool _multithreaded = false;
}

namespace ddm {
namespace internal {

void wait_breakpoint()
{
  sleep(1);
}

} // namespace internal
} // namespace ddm

void ddm::init(int * argc, char ** *argv)
{
  DDM_LOG_DEBUG("ddm::init()");

  DDM_LOG_DEBUG("ddm::init", "ddm::util::Config::init()");
  ddm::util::Config::init();

#if defined(DDM_ENABLE_THREADSUPPORT)
  DDM_LOG_DEBUG("ddm::init", "dart_init_thread()");
  dart_thread_support_level_t provided_mt;
  dart_init_thread(argc, argv, &provided_mt);
  ddm::_multithreaded = (provided_mt == DART_THREAD_MULTIPLE);
  if (!ddm::_multithreaded) {
    DDM_LOG_WARN("ddm::init",
                  "Support for multi-threading requested at compile time but "
                  "DART does not support multi-threaded access.");
  }
#else
  DDM_LOG_DEBUG("ddm::init", "dart_init()");
  dart_init(argc, argv);
#endif

  ddm::_initialized = true;

  if (ddm::util::Config::get<bool>("DDM_INIT_BREAKPOINT")) {
    if (ddm::myid() == 0) {
      int blockvar = 1;
      ddm::prevent_opt_elimination(blockvar);
      while (blockvar) {
        ddm::internal::wait_breakpoint();
      }
    }
    ddm::barrier();
  }

  DDM_LOG_DEBUG("ddm::init", "ddm::util::Locality::init()");
  ddm::util::Locality::init();
  DDM_LOG_DEBUG("ddm::init >");
}

void ddm::finalize()
{
  DDM_LOG_DEBUG("ddm::finalize()");
  // Check init status of DDM
  if (!ddm::is_initialized()) {
    // DDM has not been initalized or multiple calls of finalize, ignore:
    DDM_LOG_DEBUG("ddm::finalize", "not initialized, ignore");
    return;
  }

  // Wait for all units:
  ddm::barrier();

  // Deallocate global memory allocated in teams:
  DDM_LOG_DEBUG("ddm::finalize", "free team global memory");
  ddm::Team::finalize();

  // Wait for all units:
  ddm::barrier();

  // Finalize DDM runtime:
  DDM_LOG_DEBUG("ddm::finalize", "finalize DDM runtime");
  dart_exit();

  // Mark DDM as finalized (allow subsequent ddm::init):
  ddm::_initialized = false;

  DDM_LOG_DEBUG("ddm::finalize >");
}

bool ddm::is_initialized()
{
  return ddm::_initialized;
}

bool ddm::is_multithreaded()
{
  return ddm::_multithreaded;
}

void ddm::barrier()
{
    //static int count = 0;
    //count++;
    //std::cout << "node " << ddm::myid() << " barrier " << count << std::endl;
  ddm::Team::All().barrier();
  //std::cout << "node " << ddm::myid() << " barrier " << count <<" Done"<<std::endl;
}

ddm::global_unit_t ddm::myid()
{
  return ddm::Team::GlobalUnitID();
}

ssize_t ddm::size()
{
  return ddm::Team::All().size();
}

