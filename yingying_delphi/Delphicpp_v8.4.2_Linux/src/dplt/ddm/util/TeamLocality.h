#ifndef DDM__UTIL__TEAM_LOCALITY_H__INCLUDED
#define DDM__UTIL__TEAM_LOCALITY_H__INCLUDED

#include "../../ddm/util/Locality.h"
#include "../../ddm/util/LocalityDomain.h"
#include "../../ddm/util/UnitLocality.h"

#include "../dart-impl/dart_types.h"
#include "../dart-impl/dart_locality.h"

#include "../../ddm/algorithm/internal/String.h"

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
 * Hierarchical locality domains of a specified team.
 *
 * Usage examples:
 *
 * \code
 * ddm::Team & team = ddm::Team::All();
 *
 * ddm::util::TeamLocality tloc(team);
 *
 * // Team locality at first node, split at module scope:
 * tloc.select(".0")
 *     .split(ddm::util::Locality::Scope::Module);
 *
 * size_t num_module_parts = tloc.parts().size();
 * for (ddm::util::LocalityDomain domain : tloc.parts()) {
 *   int module_index            = domain.relative_index();
 *   int domain_max_core_mhz     = domain.hwinfo().max_cpu_mhz;
 *   int domain_min_core_threads = domain.hwinfo().min_threads;
 *   int domain_core_perf        = domain_max_core_mhz *
 *                                 domain_min_core_threads;
 *
 *   size_t num_module_units = domain.units().size();
 *   for (global_unit_t module_unit_id : domain.unit_ids()) {
 *     ddm::util::UnitLocality uloc(team, module_unit_id);
 *
 *     std::string unit_host = uloc.host();
 *     int unit_numa_id      = uloc.hwinfo().numa_id;
 *     int unit_num_cores    = uloc.hwinfo().num_cores;
 *     int unit_num_threads  = uloc.hwinfo().max_threads * unit_num_cores;
 *   }
 * }
 * \endcode
 */
class TeamLocality
{
private:
  typedef TeamLocality                    self_t;
  typedef ddm::util::Locality::Scope     Scope_t;
  typedef ddm::util::LocalityDomain      LocalityDomain_t;

public:
  /**
   * Constructor.
   * Creates new instance of \c ddm::util::TeamLocality by loading the
   * locality domain of a specified team.
   */
  TeamLocality(
    ddm::Team       & team,
    Scope_t            scope      = Scope_t::Global,
    std::string        domain_tag = ".");

  /**
   * Constructor. Creates new instance of \c ddm::util::TeamLocality for
   * a specified team and locality domain.
   */
  TeamLocality(
    ddm::Team       & team,
    LocalityDomain_t & domain);

  /**
   * Default constructor.
   */
  TeamLocality()                           = default;

  /**
   * Copy constructor.
   */
  TeamLocality(const self_t & other)       = default;

  /**
   * Assignment operator.
   */
  self_t & operator=(const self_t & other) = default;

  /**
   * The team locality domain descriptor.
   */
  inline const LocalityDomain_t & domain() const
  {
    return _domain;
  }

  /**
   * The team locality domain descriptor.
   */
  inline LocalityDomain_t & domain()
  {
    return _domain;
  }

  //inline const ddm::util::UnitLocality & unit_locality(
  //  team_unit_t unit) const
  //{
  //  return _domain.unit_locality(unit);
  //}
  //
  //inline ddm::util::UnitLocality & unit_locality(
  //  team_unit_t unit)
  //{
  //  return _domain.unit_locality(unit);
  //}

  /**
   * Split the team locality domain into the given number of parts on the
   * specified locality scope.
   * Team locality domains resulting from the split can be accessed using
   * method \c parts().
   */
  inline self_t & split(Scope_t scope, int num_split_parts = 0)
  {
    _domain.split(scope, num_split_parts);
    return *this;
  }

  /**
   * Split groups in locality domain into separate parts.
   */
  inline self_t & split_groups()
  {
    _domain.split_groups();
    return *this;
  }

  /**
   * Parts of the team locality that resulted from a previous split.
   */
  inline std::vector<LocalityDomain_t::iterator> & groups()
  {
    return _domain.groups();
  }

  /**
   * Parts of the team locality that have been created in a previous split.
   */
  inline const std::vector<LocalityDomain_t::iterator> & groups() const
  {
    return _domain.groups();
  }

  /**
   * Parts of the team locality that resulted from a previous split.
   */
  inline std::vector<LocalityDomain_t> & parts()
  {
    return _domain.parts();
  }

  /**
   * Parts of the team locality that have been created in a previous split.
   */
  inline const std::vector<LocalityDomain_t> & parts() const
  {
    return _domain.parts();
  }

  inline size_t num_nodes() const
  {
    return _domain.size();
  }

  inline size_t num_cores() const
  {
    return _domain.num_cores();
  }

  inline ddm::Team & team()
  {
    return (nullptr == _team) ? ddm::Team::Null() : *_team;
  }

  inline const ddm::Team & team() const
  {
    return (nullptr == _team) ? ddm::Team::Null() : *_team;
  }

  inline const std::vector<global_unit_t> & global_units() const
  {
    return _domain.units();
  }

  inline ddm::util::UnitLocality unit_locality(
    team_unit_t unit_id) const
  {
    return ddm::util::UnitLocality(*_team, unit_id);
  }

  inline ddm::util::UnitLocality unit_locality(
    global_unit_t unit_id) const
  {
    team_unit_t l_unit_id;
    dart_team_unit_g2l(_team->dart_id(), unit_id, &l_unit_id);
    return ddm::util::UnitLocality(*_team, l_unit_id);
  }

  inline LocalityDomain_t & group(
    const std::vector<std::string> & group_subdomain_tags)
  {
    return _domain.group(group_subdomain_tags);
  }

  inline self_t & select(
    const std::vector<std::string> & domain_tags)
  {
    _domain.select(domain_tags);
    return *this;
  }

  inline self_t & exclude(
    const std::vector<std::string> & domain_tags)
  {
    _domain.exclude(domain_tags);
    return *this;
  }

private:
  ddm::Team                        * _team          = nullptr;
  /// Parent scope of the team locality domain hierarchy.
  Scope_t                             _scope         = Scope_t::Undefined;
  /// Locality domain of the team.
  LocalityDomain_t                    _domain;

}; // class TeamLocality

}  // namespace util
}  // namespace ddm

#endif // DDM__UTIL__TEAM_LOCALITY_H__INCLUDED
