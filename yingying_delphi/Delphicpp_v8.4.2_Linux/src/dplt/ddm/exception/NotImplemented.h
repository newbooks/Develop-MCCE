#ifndef DDM__EXCEPTION__NOT_IMPLEMENTED_H_
#define DDM__EXCEPTION__NOT_IMPLEMENTED_H_

#include <stdexcept>
#include <string>

namespace ddm {
namespace exception {

class NotImplemented : public ::std::runtime_error {
public:
  NotImplemented(const ::std::string & msg)
  : ::std::runtime_error(msg) {
  }
};

} // namespace exception
} // namespace ddm

#endif // DDM__EXCEPTION__NOT_IMPLEMENTED_H_
