#ifndef DDM__LAUNCH__H__INCLUDED
#define DDM__LAUNCH__H__INCLUDED

#include <cstdint>

namespace ddm {

enum class launch : uint16_t {
/// synchronous launch policy
sync     = 0x1,
/// async launch policy
async    = 0x2
};

}


#endif // DDM__LAUNCH__H__INCLUDED
