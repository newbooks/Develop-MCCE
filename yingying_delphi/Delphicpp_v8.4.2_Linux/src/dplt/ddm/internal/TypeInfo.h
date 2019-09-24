#ifndef DDM__INTERNAL__TYPE_INFO_H__INCLUDED
#define DDM__INTERNAL__TYPE_INFO_H__INCLUDED

#include <string>
#include <typeinfo>

namespace ddm {
namespace internal {

std::string demangle(const char * typeid_name);

template <class T>
std::string typestr(const T & obj) {
  return ddm::internal::demangle(
           typeid(obj).name()
         );
}

} // namespace internal
} // namespace ddm

#endif // DDM__INTERNAL__TYPE_INFO_H__INCLUDED
