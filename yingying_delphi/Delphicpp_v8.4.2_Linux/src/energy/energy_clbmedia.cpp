/*
 * energy_clbmedia.cpp
 *
 * Author: Lin Wang, lwang3@g.clemson.edu
 * This file is function of coulombic energy calculation in linear case. 
 * It's called by energy_run function if linear case evokes and media number is greater than one.
 * 
 */

#include "energy.h"

void CDelphiEnergy::energy_clbmedia(delphi_real& fEnergy_Coulombic)
{
    int i, j;
    delphi_real dx, dy, dz;
    delphi_real fEnergy_Temp = 0.0;
    delphi_real fDistance = 0.0;
    delphi_real fEnergy_Coulombic_Temp = 0.0;

    /*
     * Define these variables in this scope for OpenMP
     */
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> prggvAtomicCrg_nValue(iCrgGridNum);
    vector<delphi_real> prgfAtomEps_Val(iCrgGridNum);

    for (j = 0; j < iCrgGridNum; j++)
    {
        prgfgCrgPoseA_nX[j] = prgfgCrgPoseA[j].nX;
        prgfgCrgPoseA_nY[j] = prgfgCrgPoseA[j].nY;
        prgfgCrgPoseA_nZ[j] = prgfgCrgPoseA[j].nZ;
        prggvAtomicCrg_nValue[j] = prggvAtomicCrg[j].nValue;
        prgfAtomEps_Val[j] = prgfAtomEps[j];
    }

    // ++++++++ Starting parallel computing and evoke OpenMP +++++++++++ //
    // Refer OpenMP website for more details about tutorials and manuals //
    // http://openmp.org/wp/openmp-specifications                        //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

#ifdef PARALLEL_OMP

#pragma omp parallel shared(prgfgCrgPoseA_nX,prgfgCrgPoseA_nY,prgfgCrgPoseA_nZ,prggvAtomicCrg_nValue,prgfAtomEps_Val) private(j,dx,dy,dz,fDistance,fEnergy_Temp)
    {
#pragma omp for reduction( + : fEnergy_Coulombic_Temp )

#endif

    for (i = 0; i < iCrgGridNum; i++)
    {
        fEnergy_Temp = 0.0;

        for (j = 0; j < iCrgGridNum; j++)
        {
            if (i != j)
            {
                /*
                 * Calculation dot product of the vector in x,y,z diction
                 */
                dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
            }
        }

        fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp / prgfAtomEps_Val[j];
    }

#ifdef PARALLEL_OMP
}   // end of #pragma omp parallel
#endif

    fEnergy_Coulombic = fEnergy_Coulombic_Temp / 2.0;

    /*
     * Clean memory occupied by these vectors
     */
    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);
    vector<delphi_real>().swap(prgfAtomEps_Val);
}

#ifdef PARALLEL_MPI


/*
 * MPI version of energy_clbmedia function
 */
void CDelphiEnergy::mpi_energy_clbmedia(delphi_real& fEnergy_Coulombic)
{
    int i, j;
    delphi_real dx, dy, dz;
    delphi_real fEnergy_Temp = 0.0;
    delphi_real fDistance = 0.0;
    delphi_real fEnergy_Coulombic_Temp = 0.0;

    /*
     * Define these variables in this scope for OpenMP
     */
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> prggvAtomicCrg_nValue(iCrgGridNum);
    vector<delphi_real> prgfAtomEps_Val(iCrgGridNum);

    if (0 == mpi_rank) /* master process */
    {
        for (j = 0; j < iCrgGridNum; j++)
        {
            prgfgCrgPoseA_nX[j] = prgfgCrgPoseA[j].nX;
            prgfgCrgPoseA_nY[j] = prgfgCrgPoseA[j].nY;
            prgfgCrgPoseA_nZ[j] = prgfgCrgPoseA[j].nZ;
            prggvAtomicCrg_nValue[j] = prggvAtomicCrg[j].nValue;
            prgfAtomEps_Val[j] = prgfAtomEps[j];
        }

        /*
         * -------------------- parallel computing starts --------------------
         */
        delphi_real * mpi_prgfgCrgPoseA_nX      = prgfgCrgPoseA_nX.data();
        delphi_real * mpi_prgfgCrgPoseA_nY      = prgfgCrgPoseA_nY.data();
        delphi_real * mpi_prgfgCrgPoseA_nZ      = prgfgCrgPoseA_nZ.data();
        delphi_real * mpi_prggvAtomicCrg_nValue = prggvAtomicCrg_nValue.data();
        delphi_real * mpi_prgfAtomEps_Val       = prgfAtomEps_Val.data();

        MPI_Bcast(mpi_prgfgCrgPoseA_nX,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nY,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nZ,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prggvAtomicCrg_nValue, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfAtomEps_Val,       iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

        /*
         * the following loop is given to slaves for parallel computing
         *
        for (i = 0; i < iCrgGridNum; i++)
        {
            fEnergy_Temp = 0.0;

            for (j = 0; j < iCrgGridNum; j++)
            {
                if (i != j)
                {

                    //Calculation dot product of the vector in x,y,z diction
                    dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                    dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                    dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                    fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                    fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
                }
            }

            fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp / prgfAtomEps_Val[j];
        }*/

        MPI_Reduce(MPI_IN_PLACE, &fEnergy_Coulombic_Temp, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);

        fEnergy_Coulombic = fEnergy_Coulombic_Temp / 2.0;
    }
    else /* slave processes */
    {
        /*
         * -------------------- parallel computing starts --------------------
         */
        delphi_real * mpi_prgfgCrgPoseA_nX      = prgfgCrgPoseA_nX.data();
        delphi_real * mpi_prgfgCrgPoseA_nY      = prgfgCrgPoseA_nY.data();
        delphi_real * mpi_prgfgCrgPoseA_nZ      = prgfgCrgPoseA_nZ.data();
        delphi_real * mpi_prggvAtomicCrg_nValue = prggvAtomicCrg_nValue.data();
        delphi_real * mpi_prgfAtomEps_Val       = prgfAtomEps_Val.data();

        MPI_Bcast(mpi_prgfgCrgPoseA_nX,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nY,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nZ,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prggvAtomicCrg_nValue, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfAtomEps_Val,       iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

        if ( iCrgGridNum <= mpi_num_workers )
        {
            if (mpi_rank <= iCrgGridNum)
            {
                i = mpi_rank - 1;

                fEnergy_Temp = 0.0;

                for (j = 0; j < iCrgGridNum; j++)
                {
                    if (i != j)
                    {
                        /*
                         * Calculation dot product of the vector in x,y,z diction
                         */
                        dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                        dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                        dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                        fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                        fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
                    }
                }

                fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp / prgfAtomEps_Val[j];
            }
        }
        /*
         * a simple parallelization but does NOT produce high speedups
         */
//        else
//        {
//            int increment = mpi_num_workers;
//
//            for (i = mpi_rank - 1; i < iCrgGridNum; i += increment)
//            {
//                fEnergy_Temp = 0.0;
//
//                for (j = 0; j < iCrgGridNum; j++)
//                {
//                    if (i != j)
//                    {
//                        /*
//                         * Calculation dot product of the vector in x,y,z diction
//                         */
//                        dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
//                        dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
//                        dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
//                        fDistance = sqrt(dx * dx + dy * dy + dz * dz);
//                        fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
//                    }
//                }
//
//                fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp / prgfAtomEps_Val[j];
//            }
//        }
        /*
         * another parallelization to improve performance by using cache more effectively
         */
        else
        {
            int mpi_mystart, mpi_myend, mpi_split = iCrgGridNum % mpi_num_workers;

            if (mpi_rank <= mpi_split)
            {
                mpi_mystart = (mpi_rank - 1) * (iCrgGridNum / mpi_num_workers + 1);
                mpi_myend   = mpi_mystart + (iCrgGridNum / mpi_num_workers + 1) - 1;
            }
            else
            {
                mpi_mystart = mpi_split * (iCrgGridNum / mpi_num_workers + 1) + (mpi_rank - mpi_split - 1) * (iCrgGridNum / mpi_num_workers);
                mpi_myend   = mpi_mystart + (iCrgGridNum / mpi_num_workers) - 1;
            }

            if (mpi_myend >= iCrgGridNum - 1) mpi_myend = iCrgGridNum - 1;

            for (i = mpi_mystart; i <= mpi_myend; i++)
            {
                fEnergy_Temp = 0.0;

                for (j = 0; j < iCrgGridNum; j++)
                {
                    if (i != j)
                    {
                        /*
                         * Calculation dot product of the vector in x,y,z diction
                         */
                        dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                        dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                        dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                        fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                        fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
                    }
                }

                fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp / prgfAtomEps_Val[j];
            }
        }

        MPI_Reduce(&fEnergy_Coulombic_Temp, &fEnergy_Coulombic_Temp, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);

    } /* end of splitting master and slave processes */

    /*
     * Clean memory occupied by these vectors
     */
    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);
    vector<delphi_real>().swap(prgfAtomEps_Val);

}

#endif
