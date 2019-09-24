// stacktrace.h (c)
// 2008, Timo Bingmann from http://idlebox.net/
// published under the WTFPL v2.0
//
// Source: https://panthema.net/2008/0901-stacktrace-demangled/

#ifndef DDM__EXCEPTION__STACK_TRACE_H_
#define DDM__EXCEPTION__STACK_TRACE_H_

#include "../../ddm/internal/Macro.h"

#ifdef __GNUC__

#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>
#include <cxxabi.h>

/**
 * Print a demangled stack backtrace of the caller function
 * to FILE* out.
 */
void ddm__print_stacktrace(
  FILE         * out        = stderr,
  unsigned int   max_frames = 63);

#else  // ifdef __GNUC__

static inline void ddm__print_stacktrace(
  FILE         * out        = stderr,
  unsigned int)
{
  fprintf(out, "No stack trace\n");
}

#endif // ifdef __GNUC__

#endif // DDM__EXCEPTION__STACK_TRACE_H_
