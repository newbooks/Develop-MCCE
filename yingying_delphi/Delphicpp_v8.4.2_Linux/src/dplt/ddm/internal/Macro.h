#ifndef DDM__INTERNAL__MACRO_H_
#define DDM__INTERNAL__MACRO_H_

/**
 * Macro parameter value string expansion.
 */
#define ddm__tostr(s) #s
/**
 * Macro parameter string expansion.
 */
#define ddm__toxstr(s) ddm__tostr(s)
/**
 * Mark variable as intentionally or potentially unused to avoid compiler
 * warning from unused variable.
 */
#define ddm__unused(x) (void)(x)

/**
 * Workaround for GCC versions that do not support the noinline attribute.
 */
#ifndef NOINLINE_ATTRIBUTE
  #ifdef __ICC
    #define NOINLINE_ATTRIBUTE __attribute__(( noinline ))
  #else
    #define NOINLINE_ATTRIBUTE
  #endif // __ICC
#endif // NOINLINE_ATTRIBUTE

#endif // DDM__INTERNAL__MACRO_H_
