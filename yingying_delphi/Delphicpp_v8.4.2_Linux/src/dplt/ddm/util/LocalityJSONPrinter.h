#ifndef DDM__UTIL__LOCALITY_JSON_PRINTER_H__INCLUDED
#define DDM__UTIL__LOCALITY_JSON_PRINTER_H__INCLUDED

#include "../../ddm/util/Locality.h"
#include "../../ddm/util/LocalityDomain.h"
#include "../../ddm/util/UnitLocality.h"
#include "../../ddm/util/TeamLocality.h"

#include <string>
#include <iostream>
#include <sstream>



namespace ddm {
namespace util {

class LocalityJSONPrinter
{
private:

  typedef LocalityJSONPrinter self_t;

public:

  LocalityJSONPrinter()
  { }

  self_t & operator<<(const std::string & str) {
    _os << str;
    return *this;
  }

  self_t & operator<<(const dart_unit_locality_t & unit_loc);

  self_t & operator<<(
    const dart_hwinfo_t & hwinfo);

  self_t & operator<<(
    const dart_domain_locality_t & domain_loc) {
    return print_domain(domain_loc.team, &domain_loc, "");
  }

  self_t & operator<<(
    dart_locality_scope_t scope);

  self_t & operator<<(
    ddm::util::Locality::Scope scope) {
    return *this << (static_cast<dart_locality_scope_t>(scope));
  }

  std::string str() const {
    return _os.str();
  }

private:

  self_t & print_domain(
    dart_team_t                    team,
    const dart_domain_locality_t * domain,
    std::string                    indent);

private:

  std::ostringstream _os;
};

template<typename T>
typename std::enable_if<
  std::is_integral<T>::value,
  LocalityJSONPrinter & >::type
operator<<(
  LocalityJSONPrinter & ljp,
  const T & v) {
  std::ostringstream os;
  os << v;
  ljp << os.str();
  return ljp;
}

} // namespace ddm
} // namespace util


#endif // DDM__UTIL__LOCALITY_JSON_PRINTER_H__INCLUDED
