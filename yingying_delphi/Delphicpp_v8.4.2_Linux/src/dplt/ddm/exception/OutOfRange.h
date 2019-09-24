#ifndef DDM__EXCEPTION__OUT_OF_RANGE_H_
#define DDM__EXCEPTION__OUT_OF_RANGE_H_

#include "../../ddm/exception/InvalidArgument.h"
#include <stdexcept>
#include <string>

namespace ddm {
namespace exception {

class OutOfRange : public ::std::out_of_range {
public:
  OutOfRange(const ::std::string & msg)
  : ::std::out_of_range(msg) {
  }
};

} // namespace exception
} // namespace ddm

#endif // DDM__EXCEPTION__OUT_OF_RANGE_H_
