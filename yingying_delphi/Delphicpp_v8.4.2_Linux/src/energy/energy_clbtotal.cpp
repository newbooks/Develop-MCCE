/*
 * energy_clbtotal.cpp
 *
 *  Created on: Feb 23, 2014
 *      Author: Lin Wang
 */

#include "energy.h"

void CDelphiEnergy::energy_clbtotal(delphi_real& fEnergy_SolvToChgIn, delphi_real& fEnergy_Coulombic)
{
    int i, j, k, n, iw, iGridOutput;
    delphi_real fCutedges_i, fCutedges_j, fCutedges_k;
    delphi_real fEnergy_Temp = 0.0, fEnergy_Temp1, fEnergy_Temp2, fEnergy_SolvToChgIn_Temp = 0.0;
    delphi_real fDistance, gX, gY, gZ;
    delphi_real goff, carica, fPhi, fTemp, c;
    vector<SGridValue<delphi_real> > sout2;

    SGridValue<delphi_real> fGridValue_Temp;

    fEnergy_Coulombic = 0.0;
    goff = (iGrid + 1.0) / 2.0;
    gX = (-goff / fScale) + fgBoxCenter.nX;
    gY = (-goff / fScale) + fgBoxCenter.nY;
    gZ = (-goff / fScale) + fgBoxCenter.nZ;

    c = fScale * fScale * fScale;

    fTemp = -2.0 * fIonStrength / c;

    /*
     * Define these variables in this scope for vectorization and omp
     */
    delphi_real dx, dy, dz;
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

    for (k = 0; k < iGrid; k++)
    {
        fCutedges_k = 1.0;
        if (k == 0 || k == iGrid - 1) fCutedges_k = 0.5;

        for (j = 0; j < iGrid; j++)
        {
            fCutedges_j = fCutedges_k;
            if (j == 0 || j == iGrid - 1) fCutedges_j = fCutedges_k * 0.5;

            for (i = 0; i < iGrid; i++)
            {
                fCutedges_i = fCutedges_j;
                if (i == 0 || i == iGrid - 1) fCutedges_i = fCutedges_j * 0.5;

                iw = k * iGrid * iGrid + j * iGrid + i;
                if (prgbDielecMap[iw])
                {
                    fPhi = prgfPhimap[iw];
                    carica = fPhi * fTemp;

                    fGridValue_Temp.nGrid.nX = i / fScale + gX + 1;
                    fGridValue_Temp.nGrid.nY = j / fScale + gY + 1;
                    fGridValue_Temp.nGrid.nZ = k / fScale + gZ + 1;
                    fGridValue_Temp.nValue = carica * fCutedges_i;

                    sout2.push_back(fGridValue_Temp);
                }
            }
        }
    }

    iGridOutput = sout2.size();

    /*
     * Define these variables in this scope for vectorization and omp
     */
    vector<delphi_real> sout2_nX(iGridOutput);
    vector<delphi_real> sout2_nY(iGridOutput);
    vector<delphi_real> sout2_nZ(iGridOutput);
    vector<delphi_real> sout2_nValue(iGridOutput);

    for (n = 0; n < iGridOutput; n++)
    {
        sout2_nX[n] = sout2[n].nGrid.nX;
        sout2_nY[n] = sout2[n].nGrid.nY;
        sout2_nZ[n] = sout2[n].nGrid.nZ;
        sout2_nValue[n] = sout2[n].nValue;
    }

    #ifdef VERBOSE
    cout << " number of g.p. in solution contributing to the energy :  " << iGridOutput << endl;
    #endif

    // ++++++++ Starting parallel computing and evoke OpenMP +++++++++++ //
    // Refer OpenMP website for more details about tutorials and manuals //
    // http://openmp.org/wp/openmp-specifications                       //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

#ifdef PARALLEL_OMP
#pragma omp parallel shared(prgfgCrgPoseA_nX,prgfgCrgPoseA_nY,prgfgCrgPoseA_nZ,sout2_nX,sout2_nY,sout2_nZ,sout2_nValue,prggvAtomicCrg_nValue,prgfAtomEps_Val) private(j,n,dx,dy,dz,fDistance,fEnergy_Temp1,fEnergy_Temp2)
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

        /*
         * calculation for the solvent contribution.
         */
        for (n = 0; n < iGridOutput; n++)
        {
            dx = prgfgCrgPoseA_nX[i] - sout2_nX[n];
            dy = prgfgCrgPoseA_nY[i] - sout2_nY[n];
            dz = prgfgCrgPoseA_nZ[i] - sout2_nZ[n];
            fDistance = sqrt(dx * dx + dy * dy + dz * dz);
            fEnergy_Temp2 += sout2_nValue[n] / fDistance;
        }

        fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
    }

#ifdef PARALLEL_OMP
}   // end of #pragma omp parallel
#endif

    fEnergy_Coulombic   = fEnergy_Temp / 2.0;
    fEnergy_SolvToChgIn = fEnergy_SolvToChgIn_Temp * 0.0006023 / (2.0 * fEpsOut);

    /*
     * Clean memory occupied by these vectors
     */
    vector<SGridValue<delphi_real> >().swap(sout2);
    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);
    vector<delphi_real>().swap(prgfAtomEps_Val);
    vector<delphi_real>().swap(sout2_nX);
    vector<delphi_real>().swap(sout2_nY);
    vector<delphi_real>().swap(sout2_nZ);
    vector<delphi_real>().swap(sout2_nValue);
}


#ifdef PARALLEL_MPI

/*
 * MPI version of energy_clbtotal function
 */
void CDelphiEnergy::mpi_energy_clbtotal(delphi_real& fEnergy_SolvToChgIn, delphi_real& fEnergy_Coulombic)
{
    int i, j, k, n, iw, iGridOutput;
    delphi_real fCutedges_i, fCutedges_j, fCutedges_k;
    delphi_real fEnergy_Temp, fEnergy_Temp1, fEnergy_Temp2, fEnergy_SolvToChgIn_Temp;
    delphi_real fDistance, gX, gY, gZ;
    delphi_real goff, carica, fPhi, fTemp, c;
    vector<SGridValue<delphi_real> > sout2;

    SGridValue<delphi_real> fGridValue_Temp;

    if (0 == mpi_rank) /* master process */
    {
        fEnergy_Coulombic        = 0.0;
        goff = (iGrid + 1.0) / 2.0;
        gX = (-goff / fScale) + fgBoxCenter.nX;
        gY = (-goff / fScale) + fgBoxCenter.nY;
        gZ = (-goff / fScale) + fgBoxCenter.nZ;

        c = fScale * fScale * fScale;

        fTemp = -2.0 * fIonStrength / c;

        /*
         * Define these variables in this scope for vectorization and omp
         */
        delphi_real dx, dy, dz;
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

        for (k = 0; k < iGrid; k++)
        {
            fCutedges_k = 1.0;
            if (k == 0 || k == iGrid - 1) fCutedges_k = 0.5;

            for (j = 0; j < iGrid; j++)
            {
                fCutedges_j = fCutedges_k;
                if (j == 0 || j == iGrid - 1) fCutedges_j = fCutedges_k * 0.5;

                for (i = 0; i < iGrid; i++)
                {
                    fCutedges_i = fCutedges_j;
                    if (i == 0 || i == iGrid - 1) fCutedges_i = fCutedges_j * 0.5;

                    iw = k * iGrid * iGrid + j * iGrid + i;
                    if (prgbDielecMap[iw])
                    {
                        fPhi = prgfPhimap[iw];
                        carica = fPhi * fTemp;

                        fGridValue_Temp.nGrid.nX = i / fScale + gX + 1;
                        fGridValue_Temp.nGrid.nY = j / fScale + gY + 1;
                        fGridValue_Temp.nGrid.nZ = k / fScale + gZ + 1;
                        fGridValue_Temp.nValue = carica * fCutedges_i;

                        sout2.push_back(fGridValue_Temp);
                    }
                }
            }
        }

        iGridOutput = sout2.size();

        /*
         * Define these variables in this scope for vectorization and omp
         */
        vector<delphi_real> sout2_nX(iGridOutput);
        vector<delphi_real> sout2_nY(iGridOutput);
        vector<delphi_real> sout2_nZ(iGridOutput);
        vector<delphi_real> sout2_nValue(iGridOutput);

        for (n = 0; n < iGridOutput; n++)
        {
            sout2_nX[n] = sout2[n].nGrid.nX;
            sout2_nY[n] = sout2[n].nGrid.nY;
            sout2_nZ[n] = sout2[n].nGrid.nZ;
            sout2_nValue[n] = sout2[n].nValue;
        }

        #ifdef VERBOSE
        cout << " number of g.p. in solution contributing to the energy :  " << iGridOutput << endl;
        #endif

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

        MPI_Bcast(&iGridOutput,              1,           MPI_INT,         0, MPI_COMM_WORLD);

        delphi_real * mpi_sout2_nX              = sout2_nX.data();
        delphi_real * mpi_sout2_nY              = sout2_nY.data();
        delphi_real * mpi_sout2_nZ              = sout2_nZ.data();
        delphi_real * mpi_sout2_nValue          = sout2_nValue.data();

        MPI_Bcast(mpi_sout2_nX,              iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout2_nY,              iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout2_nZ,              iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout2_nValue,          iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);

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


            //calculation for the solvent contribution.
            for (n = 0; n < iGridOutput; n++)
            {
                dx = prgfgCrgPoseA_nX[i] - sout2_nX[n];
                dy = prgfgCrgPoseA_nY[i] - sout2_nY[n];
                dz = prgfgCrgPoseA_nZ[i] - sout2_nZ[n];
                fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                fEnergy_Temp2 += sout2_nValue[n] / fDistance;
            }

            fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
        }*/

        MPI_Reduce(MPI_IN_PLACE, &fEnergy_Temp,             1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, &fEnergy_SolvToChgIn_Temp, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);

        fEnergy_Coulombic   = fEnergy_Temp / 2.0;
        fEnergy_SolvToChgIn = fEnergy_SolvToChgIn_Temp * 0.0006023 / (2.0 * fEpsOut);

        /*
         * Clean memory occupied by these vectors
         */
        vector<SGridValue<delphi_real> >().swap(sout2);
        vector<delphi_real>().swap(prgfgCrgPoseA_nX);
        vector<delphi_real>().swap(prgfgCrgPoseA_nY);
        vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
        vector<delphi_real>().swap(prggvAtomicCrg_nValue);
        vector<delphi_real>().swap(prgfAtomEps_Val);
        vector<delphi_real>().swap(sout2_nX);
        vector<delphi_real>().swap(sout2_nY);
        vector<delphi_real>().swap(sout2_nZ);
        vector<delphi_real>().swap(sout2_nValue);
    }
    else /* slave processes */
    {
        delphi_real dx, dy, dz;
        vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
        vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
        vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
        vector<delphi_real> prggvAtomicCrg_nValue(iCrgGridNum);
        vector<delphi_real> prgfAtomEps_Val(iCrgGridNum);

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

        MPI_Bcast(&iGridOutput,              1,           MPI_INT,         0, MPI_COMM_WORLD);

        vector<delphi_real> sout2_nX(iGridOutput);
        vector<delphi_real> sout2_nY(iGridOutput);
        vector<delphi_real> sout2_nZ(iGridOutput);
        vector<delphi_real> sout2_nValue(iGridOutput);

        delphi_real * mpi_sout2_nX     = sout2_nX.data();
        delphi_real * mpi_sout2_nY     = sout2_nY.data();
        delphi_real * mpi_sout2_nZ     = sout2_nZ.data();
        delphi_real * mpi_sout2_nValue = sout2_nValue.data();

        MPI_Bcast(mpi_sout2_nX,              iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout2_nY,              iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout2_nZ,              iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);
        MPI_Bcast(mpi_sout2_nValue,          iGridOutput, mpi_delphi_real, 0, MPI_COMM_WORLD);

        fEnergy_Temp             = 0.0;
        fEnergy_SolvToChgIn_Temp = 0.0;

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

                for (n = 0; n < iGridOutput; n++)
                {
                    dx = prgfgCrgPoseA_nX[i] - sout2_nX[n];
                    dy = prgfgCrgPoseA_nY[i] - sout2_nY[n];
                    dz = prgfgCrgPoseA_nZ[i] - sout2_nZ[n];
                    fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                    fEnergy_Temp2 += sout2_nValue[n] / fDistance;
                }

                fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
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
//                /*
//                 * calculation for the solvent contribution.
//                 */
//                for (n = 0; n < iGridOutput; n++)
//                {
//                    dx = prgfgCrgPoseA_nX[i] - sout2_nX[n];
//                    dy = prgfgCrgPoseA_nY[i] - sout2_nY[n];
//                    dz = prgfgCrgPoseA_nZ[i] - sout2_nZ[n];
//                    fDistance = sqrt(dx * dx + dy * dy + dz * dz);
//                    fEnergy_Temp2 += sout2_nValue[n] / fDistance;
//                }
//
//                fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
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

                /*
                 * calculation for the solvent contribution.
                 */
                for (n = 0; n < iGridOutput; n++)
                {
                    dx = prgfgCrgPoseA_nX[i] - sout2_nX[n];
                    dy = prgfgCrgPoseA_nY[i] - sout2_nY[n];
                    dz = prgfgCrgPoseA_nZ[i] - sout2_nZ[n];
                    fDistance = sqrt(dx * dx + dy * dy + dz * dz);
                    fEnergy_Temp2 += sout2_nValue[n] / fDistance;
                }

                fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
            }
        }

        MPI_Reduce(&fEnergy_Temp,             &fEnergy_Temp,             1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&fEnergy_SolvToChgIn_Temp, &fEnergy_SolvToChgIn_Temp, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);

        /*
         * Clean memory occupied by these vectors
         */
        vector<delphi_real>().swap(prgfgCrgPoseA_nX);
        vector<delphi_real>().swap(prgfgCrgPoseA_nY);
        vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
        vector<delphi_real>().swap(prggvAtomicCrg_nValue);
        vector<delphi_real>().swap(prgfAtomEps_Val);
        vector<delphi_real>().swap(sout2_nX);
        vector<delphi_real>().swap(sout2_nY);
        vector<delphi_real>().swap(sout2_nZ);
        vector<delphi_real>().swap(sout2_nValue);
    } /* end of splitting master and slave processes */
}

#endif
