#ifndef DDM__EXCEPTION__INVALID_ARGUMENT_H_
#define DDM__EXCEPTION__INVALID_ARGUMENT_H_

#include "../../ddm/exception/RuntimeError.h"
#include <string>

namespace ddm {
namespace exception {

class InvalidArgument : public RuntimeError {
public:
  InvalidArgument(const ::std::string & msg)
  : RuntimeError(msg) {
  }
};

} // namespace exception
} // namespace ddm

#endif // DDM__EXCEPTION__INVALID_ARGUMENT_H_
