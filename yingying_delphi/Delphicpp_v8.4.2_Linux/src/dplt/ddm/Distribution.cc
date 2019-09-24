#include "../ddm/Types.h"
#include "../ddm/Distribution.h"

const ddm::Distribution ddm::BLOCKED =
  ddm::Distribution(ddm::internal::DIST_BLOCKED, -1);

const ddm::Distribution ddm::CYCLIC =
  ddm::Distribution(ddm::internal::DIST_CYCLIC, 1);

const ddm::Distribution ddm::NONE =
  ddm::Distribution(ddm::internal::DIST_NONE, -1);

ddm::Distribution ddm::TILE(int blockSize) {
  return Distribution(ddm::internal::DIST_TILE, blockSize);
}

ddm::Distribution ddm::BLOCKCYCLIC(int blockSize) {
  return Distribution(ddm::internal::DIST_BLOCKCYCLIC, blockSize);
}

namespace ddm {

std::ostream & operator<<(
  std::ostream & os,
  const ddm::Distribution & distribution)
{
  os << "Distribution(";
  if (distribution.type == ddm::internal::DIST_TILE) {
    os << "TILE(" << distribution.blocksz << ")";
  }
  else if (distribution.type == ddm::internal::DIST_BLOCKCYCLIC) {
    os << "BLOCKCYCLIC(" << distribution.blocksz << ")";
  }
  else if (distribution.type == ddm::internal::DIST_CYCLIC) {
    os << "CYCLIC";
  }
  else if (distribution.type == ddm::internal::DIST_BLOCKED) {
    os << "BLOCKED";
  }
  else if (distribution.type == ddm::internal::DIST_NONE) {
    os << "NONE";
  }
  os << ")";
  return os;
}

}
