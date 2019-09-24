
#include "../ddm/internal/StreamConversion.h"

#include <ostream>


namespace ddm {

std::ostream & operator<<(
  std::ostream & o,
  dart_global_unit_t uid)
{
  o << uid.id;
  return o;
}

std::ostream & operator<<(
  std::ostream & o,
  dart_team_unit_t uid)
{
  o << uid.id;
  return o;
}

} // namespace ddm
