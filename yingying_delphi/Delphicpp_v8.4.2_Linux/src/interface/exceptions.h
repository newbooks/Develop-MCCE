/**
 * @file exceptions.h
 * @brief This file defines the format of warnings and exceptions
 *
 * @author Chuan Li, chuanli@clemson.edu
 */

#ifndef CEXCEPTIONS_H_
#define CEXCEPTIONS_H_

#include <iostream>
#include <sstream>
#include "environment.h"

using namespace std;
extern ostringstream cwarn;

/**
 * class CException defines the format of exceptions shown in the log file.
 *
 * Exception is a severe problem that the program does not know how to correct it. Program terminates
 * immediately with a message when an exception is thrown.
 */
class CException
{
   public:
      CException()
      {
         cerr << "\033[1;31m" << "[FATAL ERROR] ";
      }

      ~CException()
      {
         cerr << "\033[0m";
      }
};

/**
 * class CWarning defines the format of warnings shown in the log file.
 *
 * Warning is a problem that the program can correct (usually reset the variable to its default value)
 * for the user without terminating the current run. Instead, a message is shown to remind the user
 * what has been corrected.
 */

#ifndef PARALLEL_MPI
class CWarning // Program continues its run with warnings
{
   private:
                                                                                  
   public:
      static int iWarningNum;

      CWarning()
      {
         iWarningNum++;
      
         cwarn << "   " << "\033[1;34m" << "|[WARNING #" << iWarningNum << "] ";
      }

      ~CWarning()
      {
         cwarn << "\033[0m";
      }
};
#endif

#ifdef PARALLEL_MPI
class CWarning // Program continues its run with warnings
{
   private:
                                                                                  
   public:
      static int iWarningNum;
      CWarning()
      {
         iWarningNum++;
      
         cwarn << "   " << "\033[1;34m" << "[WARNING #" << iWarningNum << "] ";
      }

      ~CWarning()
      {
         cwarn << "\033[0m";
      }
};
#endif

#endif // CEXCEPTIONS_H_
