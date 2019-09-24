/*
 * energy_clbnonl.cpp
 *
 * Author: Lin Wang, lwang3@g.clemson.edu
 *
 * This file is function of coulombic energy calculation in non-linear case. 
 * It's called by energy_run function if non-linear case evokes.
 * 
 */

#include "energy.h"

void CDelphiEnergy::energy_clbnonl(delphi_real& fEnergy_Coulombic, delphi_real& fEnergy_SolvToChgIn, int& iGridOutput)
{
    int i, j, n;
    delphi_real dx, dy, dz;
    delphi_real fEnergy_Temp1, fEnergy_Temp2, fDistance = 0.0;
    delphi_real fEnergy_Temp = 0.0, fEnergy_SolvToChgIn_Temp = 0.0;
    delphi_real c = 0.0006023;

    /*
     * Define these variables in this scope for OpenMP
     */
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> sout_nX(iGridOutput);
    vector<delphi_real> sout_nY(iGridOutput);
    vector<delphi_real> sout_nZ(iGridOutput);
    vector<delphi_real> sout_nValue(iGridOutput);
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

    for (n = 0; n < iGridOutput; n++)
    {
        sout_nX[n] = sout[n].nGrid.nX;
        sout_nY[n] = sout[n].nGrid.nY;
        sout_nZ[n] = sout[n].nGrid.nZ;
        sout_nValue[n] = sout[n].nValue;
    }

    // ++++++++ Starting parallel computing and evoke OpenMP +++++++++++ //
    // Refer OpenMP website for more details about tutorials and manuals //
    // http://openmp.org/wp/openmp-specifications                       //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

#ifdef PARALLEL_OMP

#pragma omp parallel shared(prgfgCrgPoseA_nX,prgfgCrgPoseA_nY,prgfgCrgPoseA_nZ,sout_nX,sout_nY,sout_nZ,sout_nValue,prggvAtomicCrg_nValue,prgfAtomEps_Val) private(j,n,dx,dy,dz,fDistance,fEnergy_Temp1,fEnergy_Temp2)
    {
#pragma omp for reduction( + : fEnergy_Temp, fEnergy_SolvToChgIn_Temp)

#endif

    for (i = 0; i < iCrgGridNum; i++)
    {
        fEnergy_Temp1 = 0.0;
        fEnergy_Temp2 = 0.0;

        for (j = 0; j < iCrgGridNum; j++)
        {
            if (i != j)
            {
                dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                fEnergy_Temp1 += prggvAtomicCrg_nValue[j] / fDistance;
            }
        }

        fEnergy_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp1 / prgfAtomEps_Val[i];

        // Calculation for the solvent contribution
        // Array sout is declared in pointers module and allocated in nlener subroutine

        if (fIonStrength > 1.0e-6)
        {
            for (n = 0; n < iGridOutput; n++)
            {
                dx = prgfgCrgPoseA_nX[i] - sout_nX[n];
                dy = prgfgCrgPoseA_nY[i] - sout_nY[n];
                dz = prgfgCrgPoseA_nZ[i] - sout_nZ[n];
                fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                fEnergy_Temp2 += sout_nValue[n] / fDistance;
            }

            fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
        }
    }

#ifdef PARALLEL_OMP
}   // end of #pragma omp parallel
#endif

    fEnergy_Coulombic = fEnergy_Temp / 2.0;
    fEnergy_SolvToChgIn = fEnergy_SolvToChgIn_Temp * c / (2.0 * fEpsOut);

    /*
     * Clean memory occupied by these vectors
     */
    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);
    vector<delphi_real>().swap(sout_nX);
    vector<delphi_real>().swap(sout_nY);
    vector<delphi_real>().swap(sout_nZ);
    vector<delphi_real>().swap(sout_nValue);
    vector<delphi_real>().swap(prgfAtomEps_Val);
}


#ifdef PARALLEL_MPI

/*
 * MPI version of energy_clbnonl function
 */
void CDelphiEnergy::mpi_energy_clbnonl(delphi_real& fEnergy_Coulombic, delphi_real& fEnergy_SolvToChgIn, int& iGridOutput)
{
    int i, j, n;
    delphi_real dx, dy, dz;
    delphi_real fEnergy_Temp1, fEnergy_Temp2, fDistance = 0.0;
    delphi_real fEnergy_Temp = 0.0, fEnergy_SolvToChgIn_Temp = 0.0;
    delphi_real c = 0.0006023;

    /*
     * Define these variables in this scope for OpenMP
     */
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> sout_nX(iGridOutput);
    vector<delphi_real> sout_nY(iGridOutput);
    vector<delphi_real> sout_nZ(iGridOutput);
    vector<delphi_real> sout_nValue(iGridOutput);
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

        for (n = 0; n < iGridOutput; n++)
        {
            sout_nX[n] = sout[n].nGrid.nX;
            sout_nY[n] = sout[n].nGrid.nY;
            sout_nZ[n] = sout[n].nGrid.nZ;
            sout_nValue[n] = sout[n].nValue;
        }

        /*
         * -------------------- parallel computing starts --------------------
         */
        delphi_real * mpi_prgfgCrgPoseA_nX      = prgfgCrgPoseA_nX.data();
        delphi_real * mpi_prgfgCrgPoseA_nY      = prgfgCrgPoseA_nY.data();
        delphi_real * mpi_prgfgCrgPoseA_nZ      = prgfgCrgPoseA_nZ.data();
        delphi_real * mpi_prggvAtomicCrg_nValue = prggvAtomicCrg_nValue.data();
        delphi_real * mpi_prgfAtomEps_Val       = prgfAtomEps_Val.data();
        delphi_real * mpi_sout_nX               = sout_nX.data();
        delphi_real * mpi_sout_nY               = sout_nY.data();
        delphi_real * mpi_sout_nZ               = sout_nZ.data();
        delphi_real * mpi_sout_nValue           = sout_nValue.data();

        MPI_Bcast(mpi_prgfgCrgPoseA_nX,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nY,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nZ,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prggvAtomicCrg_nValue, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfAtomEps_Val,       iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout_nX,               iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout_nY,               iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout_nZ,               iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout_nValue,           iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);

        /*
         * the following loop is given to slaves for parallel computing
         *
        for (i = 0; i < iCrgGridNum; i++)
        {
            fEnergy_Temp1 = 0.0;
            fEnergy_Temp2 = 0.0;

            for (j = 0; j < iCrgGridNum; j++)
            {
                if (i != j)
                {
                    dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                    dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                    dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                    fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                    fEnergy_Temp1 += prggvAtomicCrg_nValue[j] / fDistance;
                }
            }

            fEnergy_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp1 / prgfAtomEps_Val[i];

            // Calculation for the solvent contribution
            // Array sout is declared in pointers module and allocated in nlener subroutine

            if (fIonStrength > 1.0e-6)
            {
                for (n = 0; n < iGridOutput; n++)
                {
                    dx = prgfgCrgPoseA_nX[i] - sout_nX[n];
                    dy = prgfgCrgPoseA_nY[i] - sout_nY[n];
                    dz = prgfgCrgPoseA_nZ[i] - sout_nZ[n];
                    fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                    fEnergy_Temp2 += sout_nValue[n] / fDistance;
                }

                fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
            }
        }*/

        MPI_Reduce(MPI_IN_PLACE, &fEnergy_Temp,             1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, &fEnergy_SolvToChgIn_Temp, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);

        fEnergy_Coulombic   = fEnergy_Temp / 2.0;
        fEnergy_SolvToChgIn = fEnergy_SolvToChgIn_Temp * c / (2.0 * fEpsOut);
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
        delphi_real * mpi_sout_nX               = sout_nX.data();
        delphi_real * mpi_sout_nY               = sout_nY.data();
        delphi_real * mpi_sout_nZ               = sout_nZ.data();
        delphi_real * mpi_sout_nValue           = sout_nValue.data();

        MPI_Bcast(mpi_prgfgCrgPoseA_nX,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nY,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfgCrgPoseA_nZ,      iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prggvAtomicCrg_nValue, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_prgfAtomEps_Val,       iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout_nX,               iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout_nY,               iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout_nZ,               iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout_nValue,           iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);

        if ( iCrgGridNum <= mpi_num_workers )
        {
            if (mpi_rank <= iCrgGridNum)
            {
                i = mpi_rank - 1;

                fEnergy_Temp1 = 0.0;
                fEnergy_Temp2 = 0.0;

                for (j = 0; j < iCrgGridNum; j++)
                {
                    if (i != j)
                    {
                        dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                        dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                        dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                        fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                        fEnergy_Temp1 += prggvAtomicCrg_nValue[j] / fDistance;
                    }
                }

                fEnergy_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp1 / prgfAtomEps_Val[i];

                // Calculation for the solvent contribution
                // Array sout is declared in pointers module and allocated in nlener subroutine

                if (fIonStrength > 1.0e-6)
                {
                    for (n = 0; n < iGridOutput; n++)
                    {
                        dx = prgfgCrgPoseA_nX[i] - sout_nX[n];
                        dy = prgfgCrgPoseA_nY[i] - sout_nY[n];
                        dz = prgfgCrgPoseA_nZ[i] - sout_nZ[n];
                        fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                        fEnergy_Temp2 += sout_nValue[n] / fDistance;
                    }

                    fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
                }
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
//                fEnergy_Temp1 = 0.0;
//                fEnergy_Temp2 = 0.0;
//
//                for (j = 0; j < iCrgGridNum; j++)
//                {
//                    if (i != j)
//                    {
//                        dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
//                        dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
//                        dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
//                        fDistance = sqrt(dx * dx + dy * dy + dz * dz);
//                        fEnergy_Temp1 += prggvAtomicCrg_nValue[j] / fDistance;
//                    }
//                }
//
//                fEnergy_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp1 / prgfAtomEps_Val[i];
//
//                // Calculation for the solvent contribution
//                // Array sout is declared in pointers module and allocated in nlener subroutine
//
//                if (fIonStrength > 1.0e-6)
//                {
//                    for (n = 0; n < iGridOutput; n++)
//                    {
//                        dx = prgfgCrgPoseA_nX[i] - sout_nX[n];
//                        dy = prgfgCrgPoseA_nY[i] - sout_nY[n];
//                        dz = prgfgCrgPoseA_nZ[i] - sout_nZ[n];
//                        fDistance = sqrt(dx * dx + dy * dy + dz * dz);
//                        fEnergy_Temp2 += sout_nValue[n] / fDistance;
//                    }
//
//                    fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
//                }
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
                fEnergy_Temp1 = 0.0;
                fEnergy_Temp2 = 0.0;

                for (j = 0; j < iCrgGridNum; j++)
                {
                    if (i != j)
                    {
                        dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                        dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                        dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                        fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                        fEnergy_Temp1 += prggvAtomicCrg_nValue[j] / fDistance;
                    }
                }

                fEnergy_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp1 / prgfAtomEps_Val[i];

                // Calculation for the solvent contribution
                // Array sout is declared in pointers module and allocated in nlener subroutine

                if (fIonStrength > 1.0e-6)
                {
                    for (n = 0; n < iGridOutput; n++)
                    {
                        dx = prgfgCrgPoseA_nX[i] - sout_nX[n];
                        dy = prgfgCrgPoseA_nY[i] - sout_nY[n];
                        dz = prgfgCrgPoseA_nZ[i] - sout_nZ[n];
                        fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                        fEnergy_Temp2 += sout_nValue[n] / fDistance;
                    }

                    fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
                }
            }
        }

        MPI_Reduce(&fEnergy_Temp,             &fEnergy_Temp,             1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&fEnergy_SolvToChgIn_Temp, &fEnergy_SolvToChgIn_Temp, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);

    } /* end of splitting master and slave processes */

    /*
     * Clean memory occupied by these vectors
     */
    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);
    vector<delphi_real>().swap(sout_nX);
    vector<delphi_real>().swap(sout_nY);
    vector<delphi_real>().swap(sout_nZ);
    vector<delphi_real>().swap(sout_nValue);
    vector<delphi_real>().swap(prgfAtomEps_Val);
}

#endif
