
#include "../ddm/util/TeamLocality.h"
#include "../ddm/util/Locality.h"

#include "../ddm/Exception.h"

#include "../ddm/internal/Logging.h"

#include "dart-impl/dart_types.h"
#include "dart-impl/dart_locality.h"

#include <vector>
#include <string>


ddm::util::TeamLocality::TeamLocality(
  ddm::Team                  & team,
  ddm::util::Locality::Scope   scope,
  std::string                   domain_tag)
: _team(&team),
  _scope(scope)
{
  DDM_LOG_DEBUG("TeamLocality(t,s,dt)()",
                 "team:",       team.dart_id(),
                 "scope:",      scope,
                 "domain tag:", domain_tag);

  dart_domain_locality_t * domain;
  DDM_ASSERT_RETURNS(
    dart_domain_team_locality(
      _team->dart_id(),
      domain_tag.c_str(),
      &domain),
    DART_OK);

  _domain = ddm::util::LocalityDomain(domain);

  if (_scope != Scope_t::Global) {
    DDM_LOG_DEBUG("TeamLocality(t,s,dt)", "split team domain");
    split(_scope);
  }

  DDM_LOG_DEBUG("TeamLocality(t,s,dt) >");
}

ddm::util::TeamLocality::TeamLocality(
  ddm::Team                  & team,
  ddm::util::LocalityDomain  & domain)
: _team(&team),
  _scope(domain.scope()),
  _domain(domain)
{
  DDM_LOG_DEBUG("TeamLocality(t,d)()",
                 "team:",   team.dart_id(),
                 "domain:", domain.domain_tag());

  DDM_LOG_DEBUG("TeamLocality(t,d) >");
}

