#ifndef DDM__UTIL__LOCALITY_H__
#define DDM__UTIL__LOCALITY_H__

#include "../../ddm/Init.h"

#include "../../ddm/util/Config.h"

#include "../dart-impl/dart_types.h"
#include "../dart-impl/dart_locality.h"

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cstring>


std::ostream & operator<<(
  std::ostream                 & os,
  const dart_domain_locality_t & domain_loc);

std::ostream & operator<<(
  std::ostream                 & os,
  const dart_unit_locality_t   & unit_loc);

namespace ddm {

namespace util {

class Locality
{
public:
  friend void ddm::init(int *argc, char ***argv);

public:

  typedef enum
  {
    Undefined = DART_LOCALITY_SCOPE_UNDEFINED,
    Global    = DART_LOCALITY_SCOPE_GLOBAL,
    Group     = DART_LOCALITY_SCOPE_GROUP,
    Network   = DART_LOCALITY_SCOPE_NETWORK,
    Node      = DART_LOCALITY_SCOPE_NODE,
    Module    = DART_LOCALITY_SCOPE_MODULE,
    NUMA      = DART_LOCALITY_SCOPE_NUMA,
    Unit      = DART_LOCALITY_SCOPE_UNIT,
    Package   = DART_LOCALITY_SCOPE_PACKAGE,
    Uncore    = DART_LOCALITY_SCOPE_UNCORE,
    Cache     = DART_LOCALITY_SCOPE_CACHE,
    Core      = DART_LOCALITY_SCOPE_CORE,
    CPU       = DART_LOCALITY_SCOPE_CPU
  }
  Scope;

public:

  static inline int NumNodes()
  {
    return (_team_loc == nullptr)
//         ? -1 : std::max<int>(_team_loc->num_nodes, 1);
           ? -1 : std::max<int>(_team_loc->num_domains, 1);
  }


private:
  static void init();

private:
  static dart_unit_locality_t     * _unit_loc;
  static dart_domain_locality_t   * _team_loc;

};

} // namespace util
} // namespace ddm

std::ostream & operator<<(
  std::ostream                 & os,
  ddm::util::Locality::Scope    scope);

std::ostream & operator<<(
  std::ostream                 & os,
  dart_locality_scope_t          scope);

#endif // DDM__UTIL__LOCALITY_H__
