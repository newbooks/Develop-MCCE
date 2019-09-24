/**
 * \file ddm/dart/mpi/dart_locality_priv.c
 *
 */

#include "dart_locality_priv.h"

#include "dart_types.h"
#include "dart_communication.h"
#include "dart_locality.h"
#include "dart_team_group.h"

#include "logging.h"
#include "locality.h"


dart_ret_t dart__mpi__locality_init()
{
  DART_LOG_DEBUG("dart__mpi__locality_init()");
  dart_ret_t ret;

  ret = dart__base__locality__init();
  if (ret != DART_OK) {
    DART_LOG_ERROR("dart__mpi__locality_init ! "
                   "dart__base__locality__init failed: %d", ret);
    return ret;
  }
  DART_LOG_DEBUG("dart__mpi__locality_init >");
  return DART_OK;
}

dart_ret_t dart__mpi__locality_finalize()
{
  DART_LOG_DEBUG("dart__mpi__locality_finalize()");
  dart_ret_t ret;

  ret = dart__base__locality__finalize();

  dart_barrier(DART_TEAM_ALL);

  if (ret != DART_OK) {
    DART_LOG_ERROR("dart__mpi__locality_finalize ! "
                   "dart__base__locality__finalize failed: %d", ret);
    return ret;
  }
  DART_LOG_DEBUG("dart__mpi__locality_finalize >");
  return DART_OK;
}

