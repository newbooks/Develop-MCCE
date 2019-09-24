/*
 * energy_clb.cpp
 *
 * Author: Lin Wang, lwang3@g.clemson.edu
 *
 * This file is function of coulombic energy calculation in linear case. 
 * It's called by energy_run function if linear case evokes and media number is only 1.
 * 
 */

#include "energy.h"

void CDelphiEnergy::energy_clb(delphi_real& fEnergy_Coulombic)
{
    int i, j;
    delphi_real dx, dy, dz;
    delphi_real fEnergy_Temp = 0.0;
    delphi_real fEnergy_Coulombic_Temp = 0.0;
    delphi_real fDistance = 0.0;

    /*
     * Define these variables in this scope for OpenMP
     */
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> prggvAtomicCrg_nValue(iCrgGridNum);

    for (j = 0; j < iCrgGridNum; j++)
    {
        prgfgCrgPoseA_nX[j]      = prgfgCrgPoseA[j].nX;
        prgfgCrgPoseA_nY[j]      = prgfgCrgPoseA[j].nY;
        prgfgCrgPoseA_nZ[j]      = prgfgCrgPoseA[j].nZ;
        prggvAtomicCrg_nValue[j] = prggvAtomicCrg[j].nValue;
    }

    // ++++++++ Starting parallel computing and evoke OpenMP +++++++++++ //
    // Refer OpenMP website for more details about tutorials and manuals //
    // http://openmp.org/wp/openmp-specifications                        //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

#ifdef PARALLEL_OMP

#pragma omp parallel shared(prgfgCrgPoseA_nX,prgfgCrgPoseA_nY,prgfgCrgPoseA_nZ,prggvAtomicCrg_nValue) private(j,dx,dy,dz,fDistance,fEnergy_Temp)
    {
#pragma omp for reduction( + : fEnergy_Coulombic_Temp )

#endif

    for (i = 0; i < iCrgGridNum - 1; i++)
    {
        fEnergy_Temp = 0.0;

        for (j = i + 1; j < iCrgGridNum; j++)
        {
            /*
             * Calculation dot product of the vector in x,y,z diction
             */
            dx            = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
            dy            = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
            dz            = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
            fDistance     = sqrt(dx * dx + dy * dy + dz * dz);
            fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
        }

        fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp;
    }

#ifdef PARALLEL_OMP

}   // end of #pragma omp parallel

#endif

    fEnergy_Coulombic = fEnergy_Coulombic_Temp;

    /*
     * Clean memory occupied by these vectors
     */
    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);
}


#ifdef PARALLEL_MPI

/*
 * MPI version of energy_clb function
 */
void CDelphiEnergy::mpi_energy_clb(delphi_real& fEnergy_Coulombic)
{
    int i, j;
    delphi_real dx, dy, dz;
    delphi_real fEnergy_Temp = 0.0;
    delphi_real fEnergy_Coulombic_Temp = 0.0;
    delphi_real fDistance = 0.0;

    /*
     * Define these variables in this scope for OpenMP
     */
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> prggvAtomicCrg_nValue(iCrgGridNum);

    if (0 == mpi_rank) /* master process */
    {
        for (j = 0; j < iCrgGridNum; j++)
        {
            prgfgCrgPoseA_nX[j]      = prgfgCrgPoseA[j].nX;
            prgfgCrgPoseA_nY[j]      = prgfgCrgPoseA[j].nY;
            prgfgCrgPoseA_nZ[j]      = prgfgCrgPoseA[j].nZ;
            prggvAtomicCrg_nValue[j] = prggvAtomicCrg[j].nValue;
        }

        /*
         * -------------------- parallel computing starts --------------------
         */
        delphi_real * mpi_prgfgCrgPoseA_nX      = prgfgCrgPoseA_nX.data();
        delphi_real * mpi_prgfgCrgPoseA_nY      = prgfgCrgPoseA_nY.data();
        delphi_real * mpi_prgfgCrgPoseA_nZ      = prgfgCrgPoseA_nZ.data();
        delphi_real * mpi_prggvAtomicCrg_nValue = prggvAtomicCrg_nValue.data();

        MPI_Bcast(mpi_prgfgCrgPoseA_nX,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nY,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nZ,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prggvAtomicCrg_nValue, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

        /*
         * the following loop is given to slaves for parallel computing
         *
        for (i = 0; i < iCrgGridNum - 1; i++)
        {
            fEnergy_Temp = 0.0;

            for (j = i + 1; j < iCrgGridNum; j++)
            {

                //Calculation dot product of the vector in x,y,z diction
                dx            = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                dy            = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                dz            = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                fDistance     = sqrt(dx * dx + dy * dy + dz * dz);
                fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
            }

            fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp;
        }*/

        MPI_Reduce(MPI_IN_PLACE, &fEnergy_Coulombic_Temp, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);

        fEnergy_Coulombic = fEnergy_Coulombic_Temp;
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

        MPI_Bcast(mpi_prgfgCrgPoseA_nX,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nY,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nZ,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prggvAtomicCrg_nValue, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

        if ( iCrgGridNum <= mpi_num_workers )
        {
            if (mpi_rank <= iCrgGridNum)
            {
                i = mpi_rank - 1;

                fEnergy_Temp = 0.0;

                for (j = i + 1; j < iCrgGridNum; j++)
                {
                    /*
                     * Calculation dot product of the vector in x,y,z diction
                     */
                    dx            = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                    dy            = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                    dz            = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                    fDistance     = sqrt(dx * dx + dy * dy + dz * dz);
                    fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
                }

                fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp;
            }
        }
        /*
         * a simple parallelization but does NOT produce high speedups
         */
//        else
//        {
//            int increment = mpi_num_workers;
//
//            for (i = mpi_rank - 1; i < iCrgGridNum - 1; i += increment)
//            {
//                fEnergy_Temp = 0.0;
//
//                for (j = i + 1; j < iCrgGridNum; j++)
//                {
//                    /*
//                     * Calculation dot product of the vector in x,y,z diction
//                     */
//                    dx            = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
//                    dy            = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
//                    dz            = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
//                    fDistance     = sqrt(dx * dx + dy * dy + dz * dz);
//                    fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
//                }
//
//                fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp;
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

            if (mpi_myend >= iCrgGridNum - 2) mpi_myend = iCrgGridNum - 2;

            for (i = mpi_mystart; i <= mpi_myend; i++)
            {
                fEnergy_Temp = 0.0;

                for (j = i + 1; j < iCrgGridNum; j++)
                {
                    /*
                     * Calculation dot product of the vector in x,y,z diction
                     */
                    dx            = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                    dy            = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                    dz            = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                    fDistance     = sqrt(dx * dx + dy * dy + dz * dz);
                    fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
                }

                fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp;
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
}

#endif
