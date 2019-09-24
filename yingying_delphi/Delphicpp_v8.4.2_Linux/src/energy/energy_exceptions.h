/*
 * energy_exceptions.h
 *
 *      Author: Lin Wang, lwang3@g.clemson.edu
 *
 * These warnings and exceptions should inherit from the classes CWarning and CException in order to
 * keep track the number of warnings and maintain consistent format.
 */

#ifndef ENERGY_EXCEPTIONS_H_
#define ENERGY_EXCEPTIONS_H_

#include "../interface/exceptions.h"

class CIonicCalcUnderTest : public CWarning
{
   public:
    CIonicCalcUnderTest()
      {
         cwarn << "THIS PART IS STILL UNDER TESTING. \n";
      }
};

class CReacFieldEnergyOverEST : public CWarning
{
   public:
    CReacFieldEnergyOverEST()
      {

         cwarn << "BE CAREFUL!! SELF REACTION FIELD ENERGY MIGHT HAVE BEEN OVERESTIMATED BECAUSE OF RADII CHANGES. \n";
     
      }
};

class CThreadsLessThanProcs : public CWarning
{
    public:
     CThreadsLessThanProcs(const int& iProcs, const int&iThreads)
        {
            cwarn << "PROGRAM DETECTS " << iProcs << " PROCESSORS AVAILABLE. HOWEVER, PROGRAM IS SET " << iThreads << " THREADS FOR CALCULATION. \n";
            cwarn << "FOR BEST PERFORMANCE, SET THREADS NUMBER NO LESS THAN TOTAL NUMBER OF PROCESSORS BY USING COMMAND 'export OMP_NUM_THREADS=Processors_Number'. \n";
        }   
};

class CAtomRadiusZero : public CWarning
{
   public:
    CAtomRadiusZero(const delphi_integer& iAtom)
      {

            cwarn << "CHARGED ATOM NUMBER " << setw(4) << fixed << right << iAtom << "  RADIUS CHANGED FROM ZERO TO 1.0 (AFFECTS ONLY THE SELF-REACTION ENERGY)" << endl;      

      }
};

#endif /* ENERGY_EXCEPTIONS_H_ */
