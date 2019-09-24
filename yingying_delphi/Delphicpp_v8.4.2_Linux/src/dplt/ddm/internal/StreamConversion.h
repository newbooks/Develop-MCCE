#ifndef DDM__INTERNAL__STREAM_CONVERSION_H_
#define DDM__INTERNAL__STREAM_CONVERSION_H_

#include "../../ddm/internal/Macro.h"
#include "../../ddm/internal/TypeInfo.h"

#include "../dart-impl/dart_types.h"

#include "../../ddm/Range.h"

#include <array>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <sstream>
#include <iterator>
#include <cstring>


namespace ddm {

std::ostream & operator<<(
  std::ostream & o,
  dart_global_unit_t uid);

std::ostream & operator<<(
  std::ostream & o,
  dart_team_unit_t uid);


/**
 * Write \c std::pair to output stream.
 */
template <class T1, class T2>
std::ostream & operator<<(
  std::ostream            & os,
  const std::pair<T1, T2> & p) {
  os << "(" << p.first << "," << p.second << ")";
  return os;
}


/**
 * Write elements in std::map to output stream.
 */
template <typename T1, typename T2>
std::ostream & operator<<(
  std::ostream & o,
  const std::map<T1,T2> & map)
{
  std::ostringstream ss;
  auto nelem = map.size();
  ss << "{ ";
  int i = 1;
  for (auto kv : map) {
    ss << "(" << kv.first
       << ":" << kv.second << ")"
       << (i++ < nelem ? ", " : "");
  }
  ss << " }";
  operator<<(o, ss.str());
  return o;
}

/**
 * Write elements in std::set to output stream.
 */
template <typename T>
std::ostream & operator<<(
  std::ostream & o,
  const std::set<T> & set)
{
  std::ostringstream ss;
  auto nelem = set.size();
  ss << "{ ";
  int i = 1;
  for (auto e : set) {
    ss << e
       << (i++ < nelem ? ", " : "");
  }
  ss << " }";
  operator<<(o, ss.str());
  return o;
}

/**
 * Write range of random access iterators to output stream.
 */
template <typename Range>
auto operator<<(
  std::ostream & o,
  const Range  & range)
  -> typename std::enable_if<
       (
      // type is range:
         ddm::is_range<Range>::value &&
      // type is not std::string or derivative:
         !std::is_same<Range, std::string>::value &&
         !std::is_base_of<std::string, Range>::value &&
      // range iterator type is random access:
         std::is_same<
           std::random_access_iterator_tag,
           typename std::iterator_traits<
             decltype(ddm::begin(range))>::iterator_category
         >::value
       ),
       std::ostream &
    >::type
{
  typedef typename Range::value_type value_t;

  std::ostringstream ss;
  int pos = 0;
  ss << ddm::internal::typestr(*ddm::begin(range))
     << " { ";
  for (auto it = ddm::begin(range); it != ddm::end(range); ++it, ++pos) {
    ss << static_cast<const value_t>(*it) << " ";
  }
  ss << "}";
  return operator<<(o, ss.str());
}

/**
 * Write range of non-random access iterators to output stream.
 */
template <typename Range>
auto operator<<(
  std::ostream & o,
  const Range  & range)
  -> typename std::enable_if<
       (
      // type is range:
         ddm::is_range<Range>::value &&
      // type is not std::string or derivative:
         !std::is_same<Range, std::string>::value &&
         !std::is_base_of<std::string, Range>::value &&
      // range iterator type is not random access:
         !std::is_same<
           std::random_access_iterator_tag,
           typename std::iterator_traits<
             decltype(ddm::begin(range))>::iterator_category
         >::value
       ),
       std::ostream &
    >::type
{
  typedef typename Range::value_type value_t;

  std::ostringstream ss;
  ss << ddm::internal::typestr(*ddm::begin(range))
     << " { ";
  for (auto it = ddm::begin(range); it != ddm::end(range); ++it) {
    ss << static_cast<const value_t>(*it) << " ";
  }
  ss << "}";
  return operator<<(o, ss.str());
}

} // namespace ddm

#endif // DDM__INTERNAL__STREAM_CONVERSION_H_
