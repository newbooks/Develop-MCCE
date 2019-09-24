#ifndef DDM__INTERNAL__CONFIG_H_
#define DDM__INTERNAL__CONFIG_H_

/**
 * Input for configuration file generated during build.
 * Provides platform-specific definitions.
 */

#ifdef DOXYGEN

/**
 * \defgroup{Config}
 *
 * \ingroup{Config}
 *
 * \par{Architecture-specific Definitions}
 *
 * Definition                             | Defined for                                |
 * -------------------------------------- | ------------------------------------------ |
 * <tt>DDM__ARCH__ARCH_32</tt>           | Any 32-bit architecture.                   |
 * <tt>DDM__ARCH__ARCH_64</tt>           | Any 64-bit architecture.                   |
 * <tt>DDM__ARCH__ARCH_X86_32</tt>       | Intel x86 compatible 32-bit architecture.  |
 * <tt>DDM__ARCH__ARCH_X86_64</tt>       | Intel x86 compatible 64-bit architecture.  |
 * <tt>DDM__ARCH__ARCH_ARM</tt>          | Any ARM architecture.                      |
 * <tt>DDM__ARCH__ARCH_ARMV<i>X</i></tt> | ARM architecture version <i>X</i>          |
 * &nbsp;                                 | e.g. <tt>DDM__ARCH__ARMV7</tt> for ARMv7. |
 * <tt>DDM__ARCH__CACHE_LINE_SIZE</tt>   | Width of a single cache line, in bytes.    |
 * <tt>DDM__ARCH__PAGE_SIZE</tt>         | Width of a single memory page, in bytes.   |
 * <tt>DDM__ARCH__HAS_CAS</tt>           | Atomic Compare-And-Swap supported.         |
 * <tt>DDM__ARCH__HAS_CAS_64</tt>        | CAS on 64-bit wide values supported.       |
 * <tt>DDM__ARCH__HAS_CAS_32</tt>        | CAS on 32-bit wide values supported.       |
 * <tt>DDM__ARCH__HAS_LLSC</tt>          | Load-Linked/Store-Conditional supported.   |
 * <tt>DDM__ARCH__HAS_LLSC_32</tt>       | LL/SC on 32-bit wide values supported.     |
 * <tt>DDM__ARCH__HAS_LLSC_64</tt>       | LL/SC on 64-bit wide values supported.     |
 *
 * \par{OS-specific Definitions}
 *
 * Definition                             | Defined for                                |
 * -------------------------------------- | ------------------------------------------ |
 * <tt>DDM__PLATFORM__POSIX</tt>         | POSIX-compatible platform.                 |
 * <tt>DDM__PLATFORM__LINUX</tt>         | Linux platform.                            |
 * <tt>DDM__PLATFORM__FREEBSD</tt>       | FreeBSD platform.                          |
 * <tt>DDM__PLATFORM__OSX</tt>           | Apple OSX platform.                        |
 * <tt>DDM__PLATFORM__UX</tt>            | HP-UX/Sun platform.                        |
 *
 */

#else // !DOXYGEN

// Architecture defines

#if defined(__x86_64__)
#  define DDM__ARCH__ARCH_X86_64
#  define DDM__ARCH__ARCH_X86
#  define DDM__ARCH__ARCH_64
#  define DDM__ARCH__HAS_CAS_64
#elif defined(__i386)
#  define DDM__ARCH__ARCH_X86_32
#  define DDM__ARCH__ARCH_X86
#  define DDM__ARCH__ARCH_32
#  define DDM__ARCH__HAS_CAS_32
#elif defined(__arm__)
#  define DDM__ARCH__ARCH_ARM
// ARM versions consolidated to major architecture version. 
// See: https://wiki.edubuntu.org/ARM/Thumb2PortingHowto
#  if defined(__ARM_ARCH_7__) || \
      defined(__ARM_ARCH_7R__) || \
      defined(__ARM_ARCH_7A__)
#      define DDM__ARCH__ARCH_ARMV7 1
#  endif
#  if defined(DDM__ARCH__ARCH_ARMV7) || \
      defined(__ARM_ARCH_6__) || \
      defined(__ARM_ARCH_6J__) || \
      defined(__ARM_ARCH_6K__) || \
      defined(__ARM_ARCH_6Z__) || \
      defined(__ARM_ARCH_6T2__) || \
      defined(__ARM_ARCH_6ZK__)
#      define DDM__ARCH__ARCH_ARMV6 1
#  endif
#  if defined(DDM__ARCH__ARCH_ARMV6) || \
      defined(__ARM_ARCH_5T__) || \
      defined(__ARM_ARCH_5E__) || \
      defined(__ARM_ARCH_5TE__) || \
      defined(__ARM_ARCH_5TEJ__)
#      define DDM__ARCH__ARCH_ARMV5 1
#  endif
#  if defined(DDM__ARCH__ARCH_ARMV5) || \
      defined(__ARM_ARCH_4__) || \
      defined(__ARM_ARCH_4T__)
#      define DDM__ARCH__ARCH_ARMV4 1
#  endif
#  if defined(DDM__ARCH__ARCH_ARMV4) || \
      defined(__ARM_ARCH_3__) || \
      defined(__ARM_ARCH_3M__)
#      define DDM__ARCH__ARCH_ARMV3 1
#  endif
#  if defined(DDM__ARCH__ARCH_ARMV3) || \
      defined(__ARM_ARCH_2__)
#    define DDM__ARCH__ARCH_ARMV2 1
#    define DDM__ARCH__ARCH_ARM 1
#  endif

#else
#  define DDM__ARCH__ARCH_UNKNOWN
#endif

// Atomic instructions:
//
// LL/SC:
#if defined(__ARM_ARCH_7A__)
#define DDM__ARCH__HAS_LLSC
#define DDM__ARCH__HAS_LLSC_64
#endif
// CAS:
#if defined(DDM__ARCH__HAS_CAS_64) || \
    defined(DDM__ARCH__HAS_CAS_32)
#  define DDM__ARCH__HAS_CAS
#endif
#if defined(DDM__ARCH__HAS_LLSC_64) || \
    defined(DDM__ARCH__HAS_LLSC_32)
#  define DDM__ARCH__HAS_LLSC
#endif

#if defined(DDM__ARCH__ARCH_ARM)
// Assuming 32-bit architecture for ARM:
#  define DDM__ARCH__ARCH_32
#endif

// Cache line and page size, in bytes
#if defined(DDM__ARCH__ARCH_64)
#  define DDM__ARCH__CACHE_LINE_SIZE 64
#  define DDM__ARCH__PAGE_SIZE 0x1000
#else
#  define DDM__ARCH__CACHE_LINE_SIZE 32
#  define DDM__ARCH__PAGE_SIZE 0x1000
#endif

// Platform defines

// OSX
#if defined(__MACH__) && defined(__APPLE__)
#  define DDM__PLATFORM__OSX
#endif
// UX
#if (defined(__hpux) || defined(hpux)) || \
     ((defined(__sun__) || defined(__sun) || defined(sun)) && \
      (defined(__SVR4) || defined(__svr4__)))
#  define DDM__PLATFORM__UX
#endif
// Linux
#if defined(__linux__)
#  define DDM__PLATFORM__LINUX
#  define DDM__PLATFORM__POSIX
#endif
// FreeBSD
#if defined(__FreeBSD__)
#  define DDM__PLATFORM__FREEBSD
#  define DDM__PLATFORM__POSIX
#endif

#endif // DOXYGEN

#endif // DDM__INTERNAL__CONFIG_H_
