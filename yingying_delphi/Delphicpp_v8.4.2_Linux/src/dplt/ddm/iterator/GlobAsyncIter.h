#ifndef DDM__GLOB_ASYNC_ITER_H__
#define DDM__GLOB_ASYNC_ITER_H__

#include "../../ddm/GlobAsyncRef.h"
#include "../../ddm/GlobAsyncPtr.h"
#include "../../ddm/Pattern.h"

#include "../../ddm/iterator/GlobIter.h"

#include "../dart-impl/dart_communication.h"

#include <iostream>

namespace ddm {

template <
  typename ElementType,
  class    PatternType = Pattern<1> >
class GlobAsyncIter
: public GlobIter<
           ElementType,
           PatternType,
           GlobAsyncPtr<ElementType, PatternType>,
           GlobAsyncRef<ElementType> > {
private:
  typedef GlobAsyncIter<
    ElementType,
    PatternType,
    PointerType,
    ReferenceType > self_t;

public:
  /**
   * Default constructor.
   */
  GlobAsyncIter()
  : GlobIter() {
    DDM_LOG_TRACE_VAR("GlobAsyncIter()", this->_idx);
  }

  GlobAsyncIter(
      const self_t & other) = default;
  GlobAsyncIter<ElementType, PatternType> & operator=(
      const self_t & other) = default;

  /**
   * Wait for completion of non-blocking read- and write operations
   * that have been executed on this global iterator since the last call
   * of \c wait.
   */
  void wait()
  {
    dart_flush_all(this->dart_gptr());
  }

  /**
   * Wait for local completion of non-blocking read operations that have 
   * been executed on this global iterator since the last call of \c wait.
   */
  void get()
  {
    dart_flush_all(this->dart_gptr());
  }

  /**
   * Block until all non-blocking write operations that have been executed
   * on this global iterator since the last call of \c wait have been
   * published.
   * Does not guarantee remote completion.
   */
  void push()
  {
    dart_flush_local_all(this->dart_gptr());
  }

}; // class GlobAsyncIter

}  // namespace ddm

#endif // DDM__GLOB_ASYNC_ITER_H__
