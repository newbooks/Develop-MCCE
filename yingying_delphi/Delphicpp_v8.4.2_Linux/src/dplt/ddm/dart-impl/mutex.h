#ifndef DDM_DART_BASE_MUTEX__H_
#define DDM_DART_BASE_MUTEX__H_

#include "dart_types.h"

#if defined(DART_ENABLE_THREADSUPPORT) && !defined(DART_HAVE_PTHREADS)
#error "Thread support has been enabled but PTHREADS support is not available!"
#endif

#if !defined(DART_ENABLE_THREADSUPPORT) && defined(DART_HAVE_PTHREADS)
#undef DART_HAVE_PTHREADS
#endif

#ifdef DART_HAVE_PTHREADS
#include <pthread.h>
#endif


#ifdef DART_HAVE_PTHREADS
#define DART_MUTEX_INITIALIZER { PTHREAD_MUTEX_INITIALIZER }
#else
#define DART_MUTEX_INITIALIZER { 0 }
#endif


typedef struct dart_mutex {
#ifdef DART_HAVE_PTHREADS
pthread_mutex_t mutex;
#else 
// required since C99 does not allow empty structs
// TODO: this could be used for correctness checking
char __dummy;
#endif
} dart_mutex_t;

static inline
dart_ret_t
dart__base__mutex_init(dart_mutex_t *mutex)
{
#ifdef DART_HAVE_PTHREADS
  pthread_mutex_init(&mutex->mutex, NULL);
  return DART_OK;
#else
  return DART_ERR_INVAL;
#endif
}

static inline
dart_ret_t
dart__base__mutex_lock(dart_mutex_t *mutex)
{
#ifdef DART_HAVE_PTHREADS
  pthread_mutex_lock(&mutex->mutex);
  return DART_OK;
#else
  return DART_ERR_INVAL;
#endif
}

static inline
dart_ret_t
dart__base__mutex_unlock(dart_mutex_t *mutex)
{
#ifdef DART_HAVE_PTHREADS
  pthread_mutex_unlock(&mutex->mutex);
  return DART_OK;
#else
  return DART_ERR_INVAL;
#endif
}

static inline
dart_ret_t
dart__base__mutex_trylock(dart_mutex_t *mutex)
{
#ifdef DART_HAVE_PTHREADS
  pthread_mutex_trylock(&mutex->mutex);
  return DART_OK;
#else
  return DART_ERR_INVAL;
#endif
}


static inline
dart_ret_t
dart__base__mutex_destroy(dart_mutex_t *mutex)
{
#ifdef DART_HAVE_PTHREADS
  pthread_mutex_destroy(&mutex->mutex);
  return DART_OK;
#else
  return DART_ERR_INVAL;
#endif
}

#endif /* DDM_DART_BASE_MUTEX__H_ */
