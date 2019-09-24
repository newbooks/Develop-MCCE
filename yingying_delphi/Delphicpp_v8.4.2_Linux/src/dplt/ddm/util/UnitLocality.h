#ifndef DDM__UTIL__UNIT_LOCALITY_H__INCLUDED
#define DDM__UTIL__UNIT_LOCALITY_H__INCLUDED

#include "../../ddm/util/Locality.h"
#include "../../ddm/util/LocalityDomain.h"
#include "../../ddm/util/Config.h"

#include "../../ddm/algorithm/internal/String.h"

#include "../dart-impl/dart_types.h"
#include "../dart-impl/dart_locality.h"

#include "../../ddm/Exception.h"
#include "../../ddm/Team.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <iterator>
#include <algorithm>


namespace ddm {
namespace util {

/**
 * Wrapper of a single \c dart_unit_locality_t object.
 */
class UnitLocality
{
private:
  typedef UnitLocality  self_t;

public:

  UnitLocality(
    ddm::Team   & team,
    team_unit_t    unit)
  : _team(&team)
  {
    DDM_ASSERT_RETURNS(
      dart_unit_locality(
        _team->dart_id(), unit, &_unit_locality),
      DART_OK);

    dart_domain_locality_t * team_domain;
    DDM_ASSERT_RETURNS(
      dart_domain_team_locality(
        team.dart_id(), ".", &team_domain),
      DART_OK);

    DDM_ASSERT_RETURNS(
      dart_domain_find(
        team_domain, _unit_locality->domain_tag, &_unit_domain),
      DART_OK);

    dart_domain_locality_t * node_locality = _unit_domain;
    while (node_locality->scope > DART_LOCALITY_SCOPE_NODE) {
      node_locality = node_locality->parent;
    }
    _node_domain = ddm::util::LocalityDomain(node_locality);
  }

  UnitLocality(
    global_unit_t   unit)
  : UnitLocality(ddm::Team::All(), team_unit_t(unit))
  { }

  UnitLocality()
  : UnitLocality(ddm::Team::All(), ddm::Team::All().myid())
  { }

  UnitLocality(const UnitLocality &)             = default;
  UnitLocality & operator=(const UnitLocality &) = default;

  inline const dart_hwinfo_t & hwinfo() const
  {
    DDM_ASSERT(nullptr != _unit_locality);
    return _unit_locality->hwinfo;
  }

  inline dart_hwinfo_t & hwinfo()
  {
    DDM_ASSERT(nullptr != _unit_locality);
    return _unit_locality->hwinfo;
  }

  inline dart_domain_locality_t & domain()
  {
    DDM_ASSERT(nullptr != _unit_locality);
    return *_unit_domain;
  }

  inline const dart_domain_locality_t & domain() const
  {
    DDM_ASSERT(nullptr != _unit_locality);
    return *_unit_domain;
  }

  inline ddm::Team & team()
  {
    if (nullptr == _team) {
      return ddm::Team::Null();
    }
    return *_team;
  }

  inline team_unit_t unit_id() const
  {
    return nullptr == _unit_locality
           ? UNDEFINED_TEAM_UNIT_ID
           : team_unit_t(_unit_locality->unit);
  }

  inline ddm::util::LocalityDomain & node_domain()
  {
    return _node_domain;
  }

  inline ddm::util::LocalityDomain parent()
  {
    return ddm::util::LocalityDomain(*_unit_domain->parent);
  }

  inline ddm::util::LocalityDomain parent_in_scope(
    ddm::util::Locality::Scope scope)
  {
    if (scope == ddm::util::Locality::Scope::Node) {
      return node_domain();
    }

    dart_domain_locality_t * parent_domain = _unit_domain;
    for (int rlevel = _unit_locality->hwinfo.num_scopes;
         rlevel >= 0;
         rlevel--) {
      if (parent_domain == nullptr) {
        DDM_THROW(
          ddm::exception::InvalidArgument,
          "Unit domain is undefined");
      }
      if (static_cast<int>(_unit_locality->hwinfo.scopes[rlevel].scope) <=
          static_cast<int>(scope)) {
        return ddm::util::LocalityDomain(*parent_domain);
      }
      parent_domain = parent_domain->parent;
    }
    DDM_THROW(
      ddm::exception::InvalidArgument,
      "Could not find parent domain of unit in scope " << scope);
  }

  inline std::string domain_tag() const
  {
    DDM_ASSERT(nullptr != _unit_locality);
    return _unit_domain->domain_tag;
  }

  inline std::string host() const
  {
    DDM_ASSERT(nullptr != _unit_locality);
    return _unit_locality->hwinfo.host;
  }

  inline void set_domain_tag(
    const std::string & tag)
  {
    strcpy(_unit_domain->domain_tag, tag.c_str());
  }

  inline void set_host(
    const std::string & hostname)
  {
    strcpy(_unit_locality->hwinfo.host, hostname.c_str());
  }

  inline int num_cores() const
  {
    DDM_ASSERT(nullptr != _unit_locality);
    return (_unit_locality->hwinfo.num_cores);
  }

  inline int min_threads()
  {
    return (_unit_locality == nullptr)
           ? -1 : std::max<int>(_unit_locality->hwinfo.min_threads, 1);
  }

  inline int max_threads()
  {
    return (_unit_locality == nullptr)
           ? -1 : std::max<int>(_unit_locality->hwinfo.max_threads, 1);
  }

  inline int num_threads() const
  {
    DDM_ASSERT(nullptr != _unit_locality);
    return (ddm::util::Config::get<bool>("DDM_MAX_SMT")
               ? _unit_locality->hwinfo.max_threads
               : _unit_locality->hwinfo.min_threads);
  }

  inline int num_numa() const
  {
    dart_domain_locality_t * dom = _unit_domain;
    while (dom->scope >= DART_LOCALITY_SCOPE_NUMA) {
      dom = dom->parent;
    }
    return dom->num_domains;
  }

  inline int numa_id() const
  {
    return (nullptr == _unit_locality ? -1 : _unit_locality->hwinfo.numa_id);
  }

  inline int cpu_id() const
  {
    return (nullptr == _unit_locality ? -1 : _unit_locality->hwinfo.cpu_id);
  }

  inline int cpu_mhz() const
  {
    DDM_ASSERT(nullptr != _unit_locality);
    return (_unit_locality->hwinfo.max_cpu_mhz);
  }

  inline int max_shmem_mbps() const
  {
    DDM_ASSERT(nullptr != _unit_locality);
    return (_unit_locality->hwinfo.max_shmem_mbps);
  }

  inline int max_cpu_mhz()
  {
    return (_unit_locality == nullptr)
           ? -1 : std::max<int>(_unit_locality->hwinfo.max_cpu_mhz, 1);
  }

  inline int min_cpu_mhz()
  {
    return (_unit_locality == nullptr)
           ? -1 : std::max<int>(_unit_locality->hwinfo.min_cpu_mhz, 1);
  }

  inline int cache_line_size(int cache_level)
  {
    return (_unit_locality == nullptr)
           ? 64 : std::max<int>(
                    _unit_locality->hwinfo.cache_line_sizes[cache_level],
                    64);
  }

  inline std::string hostname()
  {
    return (_unit_locality == nullptr) ? "" : _unit_locality->hwinfo.host;
  }

  /**
   * Number of threads currently available to the active unit.
   *
   * The returned value is calculated from unit locality data and hardware
   * specifications and can, for example, be used to set the \c num_threads
   * parameter of OpenMP sections:
   *
   * \code
   * #ifdef DDM_ENABLE_OPENMP
   *   auto n_threads = ddm::util::Locality::NumUnitDomainThreads();
   *   if (n_threads > 1) {
   *     #pragma omp parallel num_threads(n_threads) private(t_id)
   *     {
   *        // ...
   *     }
   * #endif
   * \endcode
   *
   * The following configuration keys affect the number of available
   * threads:
   *
   * - <tt>DDM_DISABLE_THREADS</tt>:
   *   If set, disables multi-threading at unit scope and this method
   *   returns 1.
   * - <tt>DDM_MAX_SMT</tt>:
   *   If set, virtual SMT CPUs (hyperthreads) instead of physical cores
   *   are used to determine availble threads.
   * - <tt>DDM_MAX_UNIT_THREADS</tt>:
   *   Specifies the maximum number of threads available to a single
   *   unit.
   *
   * Note that these settings may differ between hosts.
   *
   * Example for MPI:
   *
   * <tt>
   * mpirun -host node.0 -env DDM_MAX_UNIT_THREADS 4 -n 16 myprogram
   *      : -host node.1 -env DDM_MAX_UNIT_THREADS 2 -n 32 myprogram
   * </tt>
   *
   * The DDM configuration can also be changed at run time with the
   * \c ddm::util::Config interface.
   *
   * \see ddm::util::Config
   * \see ddm::util::TeamLocality
   *
   */
  inline int num_domain_threads()
  {
    auto n_threads = num_cores();
    if (ddm::util::Config::get<bool>("DDM_DISABLE_THREADS")) {
      // Threads disabled in unit scope:
      n_threads  = 1;
    } else if (ddm::util::Config::get<bool>("DDM_MAX_SMT")) {
      // Configured to use SMT (hyperthreads):
      n_threads *= max_threads();
    } else {
      // Start one thread on every physical core assigned to this unit:
      n_threads *= min_threads();
    }
    if (ddm::util::Config::is_set("DDM_MAX_UNIT_THREADS")) {
      n_threads  = std::min(ddm::util::Config::get<int>(
                              "DDM_MAX_UNIT_THREADS"),
                            n_threads);
    }
    return n_threads;
  }

private:

  ddm::Team                 * _team          = nullptr;
  dart_unit_locality_t       * _unit_locality = nullptr;
  dart_domain_locality_t     * _unit_domain   = nullptr;
  ddm::util::LocalityDomain   _node_domain;

}; // class UnitLocality

} // namespace util
} // namespace ddm

#endif // DDM__UTIL__UNIT_LOCALITY_H__INCLUDED
