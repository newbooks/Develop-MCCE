#ifndef DDM__EXCEPTION__RUNTIME_ERROR_H_
#define DDM__EXCEPTION__RUNTIME_ERROR_H_

#include <stdexcept>
#include <string>

namespace ddm {
namespace exception {

class RuntimeError : public ::std::runtime_error {
public:
  RuntimeError(const ::std::string & message)
  : ::std::runtime_error(message) {
  }
};

} // namespace exception
} // namespace ddm

#endif // DDM__EXCEPTION__RUNTIME_ERROR_H_
