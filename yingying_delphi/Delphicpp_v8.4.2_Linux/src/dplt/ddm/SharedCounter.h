#ifndef DDM__SHARED_COUNTER_H_
#define DDM__SHARED_COUNTER_H_

#include "../ddm/Array.h"

namespace ddm {

/**
 * A simple shared counter that allows atomic increment-
 * and decrement operations.
 *
 * TODO: Should probably be based on MPI_ACCUMULATE.
 */
template<typename ValueType = int>
class SharedCounter {
private:
  typedef SharedCounter<ValueType> self_t;

public:
  /**
   * Constructor.
   */
  SharedCounter()
  : _num_units(ddm::Team::All().size()),
    _myid(ddm::Team::All().myid()),
    _local_counts(_num_units)
  {
    _local_counts.local[0] = 0;
    _local_counts.barrier();
  }

  SharedCounter(ddm::Team& team)
  : _num_units(team.size()),
    _myid(team.myid()),
    _local_counts(_num_units, team)
  {
    _local_counts.local[0] = 0;
    _local_counts.barrier();
  }

  /**
   * Increment the shared counter value, atomic operation.
   */
  void inc(
    /// Increment value
    ValueType increment)
  {
    _local_counts[_myid] += increment;
  }

  /**
   * Decrement the shared counter value, atomic operation.
   */
  void dec(
    /// Decrement value
    ValueType increment)
  {
    _local_counts[_myid] -= increment;
  }

  /**
   * Read the current value of the shared counter.
   * Accumulates increment/decrement values of every unit.
   * Reading a shared is not atomic, use Team::barrier() to synchronize.
   *
   * \complexity  O(u) for \c u units in the associated team
   */
  ValueType get() const
  {
    ValueType acc = 0;
    for (team_unit_t i{0}; i < _num_units; ++i) {
      // use local access on own counter value:
      acc += (i == _myid
                ? _local_counts.local[0]
                : _local_counts[i]);
    }
    return acc;
  }

private:
  /// The number of units interacting with the counter
  size_t                 _num_units;
  /// The DART id of the unit that created this local counter intance
  team_unit_t            _myid;
  /// Buffer containing counter increments/decrements of every unit
  ddm::Array<ValueType> _local_counts;
};

} // namespace ddm

#endif // DDM__SHARED_COUNTER_H_
