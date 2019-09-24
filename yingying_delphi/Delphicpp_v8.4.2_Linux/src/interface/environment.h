/**
 * @file environment.h
 * @brief pre-processor marco's used for compiling the code
 *
 * @author Chuan Li, chuanli@clemson.edu
 */

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

/**
 * preprocessor marco to choose data types for delphi_real and delphi_integer numbers \n
 * - MX -- double      &      int [mixed double/single precision, mostly used]
 * - SP -- float       &      int [single precision]
 * - DP -- double      & long int [double precision]
 * - LD -- long double & long int [extreme case,rarely used]
 */
#define MX

/**
 * flag for additional outputs at run time
 */
//#define VERBOSE

/**
 * flag to control the precision of output numbers
 */
//#define DEVELOPER

/**
 * flag to embed delphicpp in mcce
 */
// #define MCCE

/**
 * flag to embed delphicpp in PrimePKA
 */
// #define PRIME

/**
 * flag for OMP
 */
//#define PARALLEL_OMP

/**
 * flag for MPI parallel code
 */
//#define PARALLEL_MPI

/**
 * Introduce namespace "delphi" here and used in getFunction and getStatement to avoid ambiguous reference
 * to data type "delphi_real" when compiling the code in Mac system
 */

namespace delphi
{
/*
 * single precision
 */
#ifdef SP
typedef int   delphi_integer;
typedef float delphi_real;

#ifdef PARALLEL_MPI
typedef signed int         boost_delphi_integer;
typedef float              boost_delphi_real;
#define mpi_delphi_integer MPI_INT
#define mpi_delphi_real    MPI_FLOAT
#endif

#endif

/*
 * mixed precision
 */
#ifdef MX
typedef int    delphi_integer;
typedef double delphi_real;

#ifdef PARALLEL_MPI
typedef signed int         boost_delphi_integer;
typedef double             boost_delphi_real;
#define mpi_delphi_integer MPI_INT
#define mpi_delphi_real    MPI_DOUBLE
#endif

#endif

/*
 * double precision
 */
#ifdef DP
typedef long int delphi_integer;
typedef double   delphi_real;

#ifdef PARALLEL_MPI
typedef signed long int    boost_delphi_integer;
typedef double             boost_delphi_real;
#define mpi_delphi_integer MPI_LONG
#define mpi_delphi_real    MPI_DOUBLE
#endif

#endif

/*
 * long precision
 */
#ifdef LD
typedef long int    delphi_integer;
typedef long double delphi_real;

#ifdef PARALLEL_MPI
typedef signed long int    boost_delphi_integer;
typedef long double        boost_delphi_real;
#define mpi_delphi_integer MPI_LONG
#define mpi_delphi_real    MPI_LONG_DOUBLE
#endif

#endif

/**
 * flag to debug the function(s) of reading SIZE file
 */
//#define DEBUG_IO_SIZE

#ifdef DEBUG_IO_SIZE
#define DEBUG_IO_FORCE
#endif

/**
 * flag to debug the function(s) of reading CHARGE file
 */
//#define DEBUG_IO_CHARGE

#ifdef DEBUG_IO_CHARGE
#define DEBUG_IO_FORCE
#endif

/**
 * flag to debug the function(s) of reading PDB file
 */
//#define DEBUG_IO_PDB

/**
 * flag to debug the function(s) of IDataMarshal class
 */
//#define DEBUG_DATAMARSHAL

/**
 * flag to debug the function(s) of CDelphiData class
 */
//#define DEBUG_DELPHI_MAP

/**
 * flag to debug the function(s) of CSpace class
 */
//#define DEBUG_DELPHI_SPACE

/**
 * flag to debug the function(s) of CDelphiFastSOR class
 */
//#define DEBUG_DELPHI_SOLVER

#ifdef DEBUG_DELPHI_SOLVER
#define DEBUG_DELPHI_SOLVER_MKDBSF1
#define DEBUG_DELPHI_SOLVER_MKDBSF
#define DEBUG_DELPHI_SOLVER_ITIT
#define DEBUG_DELPHI_SOLVER_RELFAC
#define DEBUG_DELPHI_SOLVER_SETBC
#define DEBUG_DELPHI_SOLVER_SETCRG
#endif

/**
 * flag to debug the function(s) of CDelphiEnergy class
 */
//#define DEBUG_DELPHI_ENERGY

/**
 * flag to debug the function(s) of CSite class
 */
//#define DEBUG_DELPHI_SITE

/**
 * flag to debug calling delphicpp in mcce
 */
//#define DEBUG_MCCE

/**
 * flag to debug constructing/destroying objects
 */
//#define DEBUG_OBJECT

}

using namespace delphi;

#endif //ENVIRONMENT_H_
