#ifndef DDM__INIT_H_
#define DDM__INIT_H_

#include "dart-impl/dart.h"
#include "../ddm/Types.h"

/**
 * \defgroup  DDMLib  High-level management functionality of the DDM library
 *
 * Functions controlling the initialization and finalization of the DDM library.
 * The library has to be initialized using \ref ddm::init before any other
 * DDM functionality can be used and should be finalized using \ref ddm::finalize
 * before program exit.
 */

namespace ddm
{
  /**
   * Initialize the DDM library and the underlying runtime system.
   *
   *
   * \param argc Pointer to the \c argc argument to \c main.
   * \param argv Pointer to the \c argv argument to \c main.
   *
   * \ingroup DDMLib
   */
  void   init(int *argc, char ***argv);

  /**
   * Finalize the DDM library and the underlying runtime system.
   *
   * \ingroup DDMLib
   */
  void   finalize();

  /**
   * Check whether DDM has been initialized already.
   *
   * \return True if DDM has been initialized successfully.
   *         False if DDM is not initialized properly or has been finalized.
   *
   * \ingroup DDMLib
   */
  bool   is_initialized();

  /**
   * Check whether DDM has been initialized with support for
   * multi-threaded access.
   *
   * \return True if DDM and the underlying runtime has been compiled
   *         with support for thread-concurrent access. False otherwise.
   *
   * \ingroup DDMLib
   */
  bool   is_multithreaded();

  /**
   * Shortcut to query the global unit ID of the calling unit.
   *
   * \return The unit ID of the calling unit relative to \ref ddm::Team::All
   *
   * \sa ddm::Team::GlobalUnitID
   * \sa ddm::Team::All
   *
   * \ingroup DDMLib
   */
  global_unit_t    myid();

  /**
   * Return the number of units in the global team.
   *
   * \return The number of units available.
   *         -1 if DDM is not initialized (anymore).
   *
   * \sa ddm::Team::size
   * \sa ddm::Team::All
   *
   * \ingroup DDMLib
   */
  ssize_t size();

  /**
   * A global barrier involving all units.
   *
   * \sa ddm::Team::barrier
   * \sa ddm::Team::All
   *
   * \ingroup DDMLib
   */
  void   barrier();
}

#endif // DDM__INIT_H_
