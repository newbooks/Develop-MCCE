#ifndef DDM__EXCEPTION__ASSERTION_FAILED_H_
#define DDM__EXCEPTION__ASSERTION_FAILED_H_

#include <stdexcept>
#include <string>

namespace ddm {
namespace exception {

class AssertionFailed : public ::std::runtime_error {
public:
  AssertionFailed(const ::std::string & msg)
  : ::std::runtime_error(msg) {
  }
};

} // namespace exception
} // namespace ddm

#endif // DDM__EXCEPTION__ASSERTION_FAILED_H_
