#ifndef DDM__UTIL__STATIC_CONFIG_H__INCLUDED
#define DDM__UTIL__STATIC_CONFIG_H__INCLUDED

/*
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !!!!! ----------- AUTO-GENERATED FILE - DO NOT EDIT ----------------!!!!!
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *
 *       Do not modify the auto-generated file `StaticConfig.h`,
 *       ensure to edit the header template `StaticConfig.h.in`.
 */

namespace ddm {
namespace util {

  static struct StaticConfig {
    bool avail_papi            = false;
    bool avail_hwloc           = false;
    bool avail_likwid          = false;
    bool avail_numa            = false;
    bool avail_plasma          = false;
    bool avail_hdf5            = false;
    bool avail_mkl             = false;
    bool avail_blas            = false;
    bool avail_lapack          = false;
    bool avail_scalapack       = false;
    /* Available Algorithms */
    bool avail_algo_summa      = false;
  } DDMConfig;

}
}

#endif // DDM__UTIL__STATIC_CONFIG_H__INCLUDED
