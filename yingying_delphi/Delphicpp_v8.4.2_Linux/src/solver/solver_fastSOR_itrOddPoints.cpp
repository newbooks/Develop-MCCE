/*
 * solver_fastSOR_itrOddPoints.cpp
 *
 *  Created on: Feb 7, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

#ifdef PARALLEL_OMP
#include <omp.h>
#endif

/*
 * -------------------- iterate over odd points (Sequential/OMP version) --------------------
 */
void CDelphiFastSOR::itrOddPoints(const int& forWhom, const int& flag) 
{
    delphi_integer n, ix, iy, iz;
    delphi_integer star, fin;
    delphi_real temp1, temp2, temp3, temp4;
    delphi_integer itemp1, itemp2, itemp3, itemp4;

    //cout << "### oddpoints phimap1: " << flag << endl;
#ifdef PARALLEL_OMP
    int omp_num_threads,omp_thread_id;

    /*
     * set number of threads = number of processors
     */
    //omp_set_num_threads(2);
    omp_set_num_threads(omp_get_num_procs());

#pragma omp parallel default(shared) private(omp_thread_id,n,ix,iy,star,fin,temp1,temp2,temp3)
    {
        delphi_integer omp_index;

        omp_thread_id = omp_get_thread_num();

        if (0 == omp_thread_id) omp_num_threads = omp_get_num_threads();

        //cout << "thread " << omp_thread_id << " of " << omp_num_threads << " is alive\n";
#endif 
        
    /* the following loops are about four times faster than the original loop over all grid points for
     * several reasons, the biggest being that we are only solving laplace's equation (unless salt is present),
     * which numerically much simpler, hence faster. we put all we leave out, back in below, ending up with
     * an equivalent calculation, but much faster.
     */
    if (fZero < abs(fIonStrength))  //----- the main loop is as below:
    {
#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif
        for (n = 1; n < iGrid - 1; n++) 
        {
            star = sta1[n];
            fin = fi1[n];
            for (ix = star; ix <= fin; ix++) 
            {
                temp1 = phimap2[ix - 1] + phimap2[(ix - 1) - 1];
                temp2 = phimap2[(ix - 1) + lat1] + phimap2[(ix - 1) - lat2];
                temp3 = phimap2[(ix - 1) + long1] + phimap2[(ix - 1) - long2];
                //phimap1[ix-1] = phimap1[ix-1]*om1 + (qmap1[ix-1]+temp1+temp2+temp3)*prgfSaltMap1[ix-1];
                phimap1[ix - 1] = phimap1[ix - 1] * om1 + (qmap1[ix - 1] + temp1 + temp2 + temp3) * prgfSaltMap1[ix - 1];
            }
        }
    } 
    else //----- if there is no salt then the main loop is executed without sf saving about 15% in execution time
    {
#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif
        for (n = 1; n < iGrid - 1; n++) 
        {
            star = sta1[n];
            fin = fi1[n];
            for (ix = star; ix <= fin; ix++) 
            {
                temp1 = phimap2[ix - 1] + phimap2[(ix - 1) - 1];
                temp2 = phimap2[(ix - 1) + lat1] + phimap2[(ix - 1) - lat2];
                temp3 = phimap2[(ix - 1) + long1] + phimap2[(ix - 1) - long2];
                phimap1[ix - 1] = phimap1[ix - 1] * om1 + (temp1 + temp2 + temp3) * sixth;
            }
        }
    }

#ifdef PARALLEL_OMP
    //#pragma omp barrier
#endif      
    
    /*
     * first we add back the dielectric boundary points, by recalculating them individually. note this is still
     * vectorised by means of a gathering load by the compiler.
     */
    if (iGaussian != 0)
    {
        if (fZero < abs(fIonStrength)) // If there is ion, Gaussian
        {
#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif
            for (n = 0; n < iDielecBndyEven; n++) 
            {
                ix = prgiBndyDielecIndex[n];

                //We need to recalculate the boudary points
                //Here we only calculate the pure linear part, and then add back the nonliear part

                delphi_real eps1 = gaussianBoundaryDielec[n][0];
                delphi_real eps2 = gaussianBoundaryDielec[n][1];
                delphi_real eps3 = gaussianBoundaryDielec[n][2];
                delphi_real eps4 = gaussianBoundaryDielec[n][3];
                delphi_real eps5 = gaussianBoundaryDielec[n][4];
                delphi_real eps6 = gaussianBoundaryDielec[n][5];

                delphi_real phi1 = phimap2[(ix - 1) - 1];
                delphi_real phi2 = phimap2[ix - 1];
                delphi_real phi3 = phimap2[(ix - 1) - lat2];
                delphi_real phi4 = phimap2[(ix - 1) + lat1];
                delphi_real phi5 = phimap2[(ix - 1) - long2];
                delphi_real phi6 = phimap2[(ix - 1) + long1];

                delphi_real myLastPhi = phimap1[ix - 1] - (qmap1[ix - 1] + temp1 + temp2 + temp3) * prgfSaltMap1[ix - 1];

                delphi_real myDensity = gaussianBoundaryDensity[n];

                delphi_real myExpSolvE = calcExpSolvE(myDensity);

                delphi_real myNonlinearCorrection = gaussianBoundaryNonlinear[n];

                delphi_real numerator = (eps1 * phi1 + eps2 * phi2 + eps3 * phi3 + eps4 * phi4 + eps5 * phi5 + eps6 * phi6) / fEPKT;
                delphi_real demonimator = (eps1 + eps2 + eps3 + eps4 + eps5 + eps6) / fEPKT + fDebFct * myExpSolvE;

                phimap1[ix - 1] = myLastPhi + (numerator / demonimator + myNonlinearCorrection) * (1 - (om1));
            }
        } 
        else  //if there is no ion, Gaussian
        {
#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif
            for (n = 0; n < iDielecBndyEven; n++) 
            {
                ix = prgiBndyDielecIndex[n];

                //We need to recalculate the boudary points
                //Here we only calculate the pure linear part, and then add back the nonliear part

                delphi_real eps1 = gaussianBoundaryDielec[n][0];
                delphi_real eps2 = gaussianBoundaryDielec[n][1];
                delphi_real eps3 = gaussianBoundaryDielec[n][2];
                delphi_real eps4 = gaussianBoundaryDielec[n][3];
                delphi_real eps5 = gaussianBoundaryDielec[n][4];
                delphi_real eps6 = gaussianBoundaryDielec[n][5];

                delphi_real phi1 = phimap2[(ix - 1) - 1];
                delphi_real phi2 = phimap2[ix - 1];
                delphi_real phi3 = phimap2[(ix - 1) - lat2];
                delphi_real phi4 = phimap2[(ix - 1) + lat1];
                delphi_real phi5 = phimap2[(ix - 1) - long2];
                delphi_real phi6 = phimap2[(ix - 1) + long1];

                delphi_real myLastPhi = phimap1[ix - 1] - (phi1 + phi2 + phi3 + phi4 + phi5 + phi6) * sixth;

                delphi_real numerator   = eps1 * phi1 + eps2 * phi2 + eps3 * phi3 + eps4 * phi4 + eps5 * phi5 + eps6 * phi6;
                delphi_real demonimator = eps1 + eps2 + eps3 + eps4 + eps5 + eps6;

                phimap1[ix - 1] = myLastPhi + (numerator / demonimator) * (1 - (om1));
            }
        }
    } 
    else // if not Gaussian
    {
#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif
        for (n = 0; n < iDielecBndyEven; n++) 
        {
            ix = prgiBndyDielecIndex[n];
            temp1 = phimap2[(ix - 1) - 1] * prgfBndyDielec[n][0] + phimap2[ix - 1] * prgfBndyDielec[n][1];
            temp2 = phimap2[(ix - 1) - lat2] * prgfBndyDielec[n][2] + phimap2[(ix - 1) + lat1] * prgfBndyDielec[n][3];
            temp3 = phimap2[(ix - 1) - long2] * prgfBndyDielec[n][4] + phimap2[(ix - 1) + long1] * prgfBndyDielec[n][5];
            phimap1[ix - 1] += temp1 + temp2 + temp3;
        }
    }  
    
    /*
     * Now reset boundary values altered in above loops.
     */
#ifdef PARALLEL_OMP
    star = (iGrid+1)/2; fin = (iGrid*(iGrid-1)-2)/2; omp_index = iGrid*(iGrid+1)/2-iGrid+1; //iy = iGrid*(iGrid+1)/2-iGrid+1;
#pragma omp for schedule(auto)
    for (n = 0; n < fin-star+1; n++)
    {
        iy = omp_index+(n+1)*iGrid;
        phimap1[iy-1] = bndx1[n];
        phimap1[iy+((iGrid+1)/2-1)-1] = bndx2[n];
    }
#else
    star = (iGrid + 1) / 2;
    fin = (iGrid * (iGrid - 1) - 2) / 2;
    iy = iGrid * (iGrid + 1) / 2 - iGrid + 1;
    for (n = 0; n < fin - star + 1; n++) 
    {
        iy = iy + iGrid;
        phimap1[iy - 1] = bndx1[n];
        phimap1[iy + ((iGrid + 1) / 2 - 1) - 1] = bndx2[n];
    }
#endif 
    
    /*
     * next we add back an adjustment to all the charged grid points due to the charge assigned. the compiler
     * directive just reassures the vector compiler that all is well as far as recurrence is concerned, i.e. it
     * would think there is a recurrence below, where as in fact there is none.
     */
    if (0 != forWhom) 
    {
        if (iGaussian != 0)
        {
            if (fZero < abs(fIonStrength)) // If there is ion, Gaussian
            {
#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif
                for (n = 0; n < iCrgedGridEven; n++) 
                {
                    ix = prgiCrgPose[n];

                    delphi_real eps1 = gaussianChargeDielec[n][0];
                    delphi_real eps2 = gaussianChargeDielec[n][1];
                    delphi_real eps3 = gaussianChargeDielec[n][2];
                    delphi_real eps4 = gaussianChargeDielec[n][3];
                    delphi_real eps5 = gaussianChargeDielec[n][4];
                    delphi_real eps6 = gaussianChargeDielec[n][5];

                    delphi_real myDensity = gaussianChargeDensity[n];
                    delphi_real myCharge  = prgfCrgValG[n];

                    delphi_real myExpSolvE = calcExpSolvE(myDensity);

                    delphi_real myNonlinearCorrection = gaussianChargeNonlinear[n];

                    delphi_real numerator   = myCharge * f4Pi * fScale;
                    delphi_real demonimator = (eps1 + eps2 + eps3 + eps4 + eps5 + eps6) / fEPKT + fDebFct * myExpSolvE;

                    phimap1[ix - 1] = phimap1[ix - 1] + (numerator / demonimator + myNonlinearCorrection) * (1 - (om1));
                }
            } 
            else // If there is no ion, Gaussian
            {
#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif
                for (n = 0; n < iCrgedGridEven; n++) 
                {
                    ix = prgiCrgPose[n];

                    delphi_real eps1 = gaussianChargeDielec[n][0];
                    delphi_real eps2 = gaussianChargeDielec[n][1];
                    delphi_real eps3 = gaussianChargeDielec[n][2];
                    delphi_real eps4 = gaussianChargeDielec[n][3];
                    delphi_real eps5 = gaussianChargeDielec[n][4];
                    delphi_real eps6 = gaussianChargeDielec[n][5];

                    delphi_real myCharge = prgfCrgValG[n];

                    delphi_real numerator = myCharge * f4Pi * fScale;
                    delphi_real demonimator = (eps1 + eps2 + eps3 + eps4 + eps5 + eps6) / fEPKT;

                    phimap1[ix - 1] = phimap1[ix - 1] + (numerator / demonimator) * (1 - (om1));
                }
            }
        }
        else 
        {
#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif
            for (n = 0; n < iCrgedGridEven; n++) 
            {
                ix = prgiCrgPose[n];
                phimap1[ix - 1] += prgfCrgValA[n];
            }
        }
    }

#ifdef PARALLEL_OMP
} // end of #pragma omp parallel
#endif
    
    /*
     * if periodic boundary condition option, force periodicity using wrap around update of boundary values:
     *    2nd slice-->last
     *    last-1 slice-->first
     */
    if (rgbPeriodicBndy[2]) //----- z periodicity
    {
        for (iz = 0; iz < (iGrid - 2) * (iGrid - 2); iz += 2) 
        {
            temp1  = ibndz[iz];
            itemp1 = (delphi_integer) temp1;
            temp2  = temp1 + idif1z;
            itemp2 = (delphi_integer) temp2;
            temp3  = temp2 + inc1za;
            itemp3 = (delphi_integer) temp3;
            temp4  = temp1 + inc1zb;
            itemp4 = (delphi_integer) temp4;
            phimap1[itemp1 - 1] = phimap2[itemp2 - 1];
            phimap1[itemp3 - 1] = phimap2[itemp4 - 1];
        }
    }

    if (rgbPeriodicBndy[1]) //----- y periodicity
    {
        for (iy = 0; iy < (iGrid - 2) * (iGrid - 2); iy += 2) 
        {
            temp1  = ibndy[iy];
            itemp1 = (delphi_integer) temp1;
            temp2  = temp1 + idif1y;
            itemp2 = (delphi_integer) temp2;
            temp3  = temp2 + inc1ya;
            itemp3 = (delphi_integer) temp3;
            temp4  = temp1 + inc1yb;
            itemp4 = (delphi_integer) temp4;
            phimap1[itemp1 - 1] = phimap2[itemp2 - 1];
            phimap1[itemp3 - 1] = phimap2[itemp4 - 1];
        }
    }

    if (rgbPeriodicBndy[0]) //----- x periodicity
    {
        for (ix = 0; ix < (iGrid - 2) * (iGrid - 2); ix += 2) 
        {
            temp1  = ibndx[ix];
            itemp1 = (delphi_integer) temp1;
            temp2  = temp1 + idif1x;
            itemp2 = (delphi_integer) temp2;
            temp3  = temp2 + inc1xa;
            itemp3 = (delphi_integer) temp3;
            temp4  = temp1 + inc1xb;
            itemp4 = (delphi_integer) temp4;
            phimap1[itemp1 - 1] = phimap2[itemp2 - 1];
            phimap1[itemp3 - 1] = phimap2[itemp4 - 1];
        }
    }
    
//    if (1 == flag)
//    {
//        string strTestFile = "rank1_solver_itrOddPoints.dat";
//        ofstream ofTestStream(strTestFile.c_str());
//        ofTestStream << boolalpha;
//        ofTestStream << fixed << setprecision(7);
//
//        ofTestStream << "flag = " << flag << endl;
//        
//        ix = 0;
//        for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
//        {
//            ofTestStream << "mpi_phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
//        {
//            ofTestStream << "mpi_phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ofTestStream.close();
//    } 
}

#ifdef PARALLEL_MPI

/*
 * -------------------- iterate over odd points (MPI version) --------------------
 */
void CDelphiFastSOR::mpi_itrOddPoints(const int& forWhom, const int& flag)
{
    delphi_integer n, ix, iy, iz;
    delphi_integer star, fin;
    delphi_real temp1, temp2, temp3, temp4;
    delphi_integer itemp1, itemp2, itemp3, itemp4;
    
    if (0 == mpi_rank) /* master process */
    {
        /*
         * periodic boundary conditions - master scatters phimap1 with updated boundary values to all slave processes
         */
        if (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2])
        {
            MPI_Scatterv(mpi_phimap1, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real, 
                         MPI_IN_PLACE, 0, mpi_delphi_real, 
                         0, MPI_COMM_WORLD);
            MPI_Scatterv(mpi_phimap2, mpi_sendcounts2l, mpi_senddispls2l, mpi_delphi_real, 
                         MPI_IN_PLACE, 0, mpi_delphi_real, 
                         0, MPI_COMM_WORLD);
            MPI_Scatterv(mpi_phimap2, mpi_sendcounts2r, mpi_senddispls2r, mpi_delphi_real, 
                         MPI_IN_PLACE, 0, mpi_delphi_real, 
                         0, MPI_COMM_WORLD);
        }

        /*
         * waiting for slave processes finish their work ......
         */
        MPI_Barrier (MPI_COMM_WORLD);

        /*
         * periodic boundary conditions - master collects computed phimap from all slave processes and resets boundary values
         */
        if (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2])
        {
            /*
             * master process collects computed phimap1
             */
            MPI_Gatherv(MPI_IN_PLACE, 0, mpi_delphi_real, mpi_phimap1, mpi_recvcounts1, mpi_recvdispls1, 
                        mpi_delphi_real, 0, MPI_COMM_WORLD);

            /*
             * master process resets boundary values in phimap2 with updated phimap1
             * if periodic boundary condition option force periodicity using wrap around update of boundary values:
             *    2nd slice-->last    last-1 slice-->first
             */
            if (rgbPeriodicBndy[2]) //----- z periodicity
            {
                for (iz = 1; iz < (iGrid - 2) * (iGrid - 2); iz += 2) 
                {
                    temp1  = ibndz[iz];
                    itemp1 = (delphi_integer) temp1;
                    temp2  = temp1 + idif1z;
                    itemp2 = (delphi_integer) temp2;
                    temp3  = temp2 + inc1za;
                    itemp3 = (delphi_integer) temp3;
                    temp4  = temp1 + inc1zb;
                    itemp4 = (delphi_integer) temp4;
                    phimap1[itemp1 - 1] = phimap2[itemp2 - 1];
                    phimap1[itemp3 - 1] = phimap2[itemp4 - 1];
                }
            }

            if (rgbPeriodicBndy[1]) //----- y periodicity
            {
                for (iy = 1; iy < (iGrid - 2) * (iGrid - 2); iy += 2) 
                {
                    temp1  = ibndy[iy];
                    itemp1 = (delphi_integer) temp1;
                    temp2  = temp1 + idif1y;
                    itemp2 = (delphi_integer) temp2;
                    temp3  = temp2 + inc1ya;
                    itemp3 = (delphi_integer) temp3;
                    temp4  = temp1 + inc1yb;
                    itemp4 = (delphi_integer) temp4;
                    phimap1[itemp1 - 1] = phimap2[itemp2 - 1];
                    phimap1[itemp3 - 1] = phimap2[itemp4 - 1];
                }
            }

            if (rgbPeriodicBndy[0]) //----- x periodicity
            {
                for (ix = 1; ix < (iGrid - 2) * (iGrid - 2); ix += 2) 
                {
                    temp1  = ibndx[ix];
                    itemp1 = (delphi_integer) temp1;
                    temp2  = temp1 + idif1x;
                    itemp2 = (delphi_integer) temp2;
                    temp3  = temp2 + inc1xa;
                    itemp3 = (delphi_integer) temp3;
                    temp4  = temp1 + inc1xb;
                    itemp4 = (delphi_integer) temp4;
                    phimap1[itemp1 - 1] = phimap2[itemp2 - 1];
                    phimap1[itemp3 - 1] = phimap2[itemp4 - 1];
                }
            }

            /*
             * master process scatters phimap2 with corrected boundary values to all slave processes
             */
            MPI_Scatterv(mpi_phimap1, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real, 
                         MPI_IN_PLACE, 0, mpi_delphi_real, 0, MPI_COMM_WORLD);
        }
    }
    else /* slave processes */
    {
        /*
         * periodic boundary conditions
         *
         * slave processes receive phimap1 with updated boundary values from master process
         */
        if (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2])
        {
            MPI_Scatterv(mpi_phimap1 + mpi_wrphimap1_start, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real, 
                         mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_wrlen1[mpi_rank], mpi_delphi_real, 
                         0, MPI_COMM_WORLD);
            MPI_Scatterv(mpi_phimap2 + mpi_wrphimap2_start, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real, 
                         mpi_phimap2 + mpi_wrstar1[mpi_rank] - (iGrid * iGrid + 1) / 2 - 1, mpi_wrstar2[mpi_rank] - (mpi_wrstar1[mpi_rank] - (iGrid * iGrid + 1) / 2), mpi_delphi_real, 
                         0, MPI_COMM_WORLD);
            MPI_Scatterv(mpi_phimap2 + mpi_wrphimap2_start, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real, 
                         mpi_phimap2 + mpi_wrfinl2[mpi_rank], mpi_wrfinl1[mpi_rank] + (iGrid * iGrid - 1) / 2 - mpi_wrfinl2[mpi_rank], mpi_delphi_real, 
                         0, MPI_COMM_WORLD);
        }
        /*
         * non-periodic boundary conditions
         *
         * slaves exchange boundary values
         */
        else
        {
            if (1 < mpi_num_workers) /* more than one slave processes */
            {
                if (1 == mpi_rank) /* the 1st slave process */
                {
                    MPI_Win_post(mpi_postgroup, 0, mpi_towin2);
                    MPI_Win_wait(mpi_towin2);
                }
                else if (1 < mpi_rank && mpi_rank < mpi_num_workers) /* the slave processes in between */
                {
                    MPI_Win_post(mpi_postgroup, 0, mpi_towin2);
                    
                    MPI_Win_start(mpi_startgroup, 0, mpi_towin2);
                    MPI_Put(mpi_phimap2 + mpi_wrstar2[mpi_rank] - 1, mpi_ltosize2, mpi_delphi_real, 
                            mpi_startrank, mpi_todispl2, mpi_ltosize2, mpi_delphi_real, mpi_towin2);
                    MPI_Get(mpi_phimap2 + mpi_lfromstart2, mpi_lfromsize2, mpi_delphi_real, 
                            mpi_startrank, mpi_zerodispl, mpi_lfromsize2, mpi_delphi_real, mpi_towin2);
                    MPI_Win_complete(mpi_towin2);
                    
                    MPI_Win_wait(mpi_towin2);
                }
                else /* the last slave process */
                {
                    MPI_Win_start(mpi_startgroup, 0, mpi_towin2);
                    MPI_Put(mpi_phimap2 + mpi_wrstar2[mpi_rank] - 1, mpi_ltosize2, mpi_delphi_real, 
                            mpi_startrank, mpi_todispl2, mpi_ltosize2, mpi_delphi_real, mpi_towin2);
                    MPI_Get(mpi_phimap2 + mpi_lfromstart2, mpi_lfromsize2, mpi_delphi_real, 
                            mpi_startrank, mpi_zerodispl, mpi_lfromsize2, mpi_delphi_real, mpi_towin2);
                    MPI_Win_complete(mpi_towin2);
                }
            }
        } // ---------- end of periodic and non-periodic boundary conditions
        
        /*
         * waiting for slave processes finish their work ......
         */
        MPI_Barrier (MPI_COMM_WORLD);

        /* the following loops are about four times faster than the original loop over all grid points for
         * several reasons, the biggest being that we are only solving laplace's equation (unless salt is present),
         * which numerically much simpler, hence faster. we put all we leave out, back in below, ending up with
         * an equivalent calculation, but much faster.
         */
        if (fZero < abs(fIonStrength))  //----- the main loop is as below:
        {
            for (n = 1; n <= mpi_wrnstafi1[mpi_rank]; n++)
            {
                star = sta1[n];
                fin  =  fi1[n];
                for (ix = star; ix <= fin; ix++)
                {
                    temp1 = mpi_phimap2[(ix - 1)        ] + mpi_phimap2[(ix - 1) - 1    ];
                    temp2 = mpi_phimap2[(ix - 1) + lat1 ] + mpi_phimap2[(ix - 1) - lat2 ];
                    temp3 = mpi_phimap2[(ix - 1) + long1] + mpi_phimap2[(ix - 1) - long2];
                    mpi_phimap1[ix - 1] = mpi_phimap1[ix - 1] * om1 + (mpi_qmap1[ix - 1] + temp1 + temp2 + temp3) * mpi_sf1[ix - 1];
                }
            }
        }
        else //----- if there is no salt then the main loop is executed without sf saving about 15% in execution time
        {
            for (n = 1; n <= mpi_wrnstafi1[mpi_rank]; n++)
            {
                star = sta1[n];
                fin  =  fi1[n];
                for (ix = star; ix <= fin; ix++)
                {
                    temp1 = mpi_phimap2[(ix - 1)        ] + mpi_phimap2[(ix - 1) - 1    ];
                    temp2 = mpi_phimap2[(ix - 1) + lat1 ] + mpi_phimap2[(ix - 1) - lat2 ];
                    temp3 = mpi_phimap2[(ix - 1) + long1] + mpi_phimap2[(ix - 1) - long2];
                    mpi_phimap1[ix - 1] = mpi_phimap1[ix - 1] * om1 + (temp1 + temp2 + temp3) * sixth;
                }
            }
        }
        
        /*
         * first we add back the dielectric boundary points, by recalculating them individually. note this is still
         * vectorized by means of a gathering load by the compiler.
         */
        
        /*
         * Gaussian based runs are not parallelized yet!
         */
        if (iGaussian != 0)
        {
//          if (fZero < abs(fIonStrength)) // If there is ion, Gaussian
//          {
//              for (n = 0; n < iDielecBndyEven; n++) 
//              {
//                  ix = prgiBndyDielecIndex[n];
//
//                  //We need to recalculate the boudary points
//                  //Here we only calculate the pure linear part, and then add back the nonliear part
//
//                  delphi_real eps1 = gaussianBoundaryDielec[n][0];
//                  delphi_real eps2 = gaussianBoundaryDielec[n][1];
//                  delphi_real eps3 = gaussianBoundaryDielec[n][2];
//                  delphi_real eps4 = gaussianBoundaryDielec[n][3];
//                  delphi_real eps5 = gaussianBoundaryDielec[n][4];
//                  delphi_real eps6 = gaussianBoundaryDielec[n][5];
//
//                  delphi_real phi1 = phimap2[(ix - 1) - 1];
//                  delphi_real phi2 = phimap2[ix - 1];
//                  delphi_real phi3 = phimap2[(ix - 1) - lat2];
//                  delphi_real phi4 = phimap2[(ix - 1) + lat1];
//                  delphi_real phi5 = phimap2[(ix - 1) - long2];
//                  delphi_real phi6 = phimap2[(ix - 1) + long1];
//
//                  delphi_real myLastPhi = phimap1[ix - 1] - (qmap1[ix - 1] + temp1 + temp2 + temp3) * prgfSaltMap1[ix - 1];
//
//                  delphi_real myDensity = gaussianBoundaryDensity[n];
//
//                  delphi_real myExpSolvE = calcExpSolvE(myDensity);
//
//                  delphi_real myNonlinearCorrection = gaussianBoundaryNonlinear[n];
//
//                  delphi_real numerator = (eps1 * phi1 + eps2 * phi2 + eps3 * phi3 + eps4 * phi4 + eps5 * phi5 + eps6 * phi6) / fEPKT;
//                  delphi_real demonimator = (eps1 + eps2 + eps3 + eps4 + eps5 + eps6) / fEPKT + fDebFct * myExpSolvE;
//
//                  phimap1[ix - 1] = myLastPhi + (numerator / demonimator + myNonlinearCorrection) * (1 - (om1));
//              }
//          } 
//          else  //if there is no ion, Gaussian
//          {
//              for (n = 0; n < iDielecBndyEven; n++) 
//              {
//                  ix = prgiBndyDielecIndex[n];
//
//                  //We need to recalculate the boudary points
//                  //Here we only calculate the pure linear part, and then add back the nonliear part
//
//                  delphi_real eps1 = gaussianBoundaryDielec[n][0];
//                  delphi_real eps2 = gaussianBoundaryDielec[n][1];
//                  delphi_real eps3 = gaussianBoundaryDielec[n][2];
//                  delphi_real eps4 = gaussianBoundaryDielec[n][3];
//                  delphi_real eps5 = gaussianBoundaryDielec[n][4];
//                  delphi_real eps6 = gaussianBoundaryDielec[n][5];
//
//                  delphi_real phi1 = phimap2[(ix - 1) - 1];
//                  delphi_real phi2 = phimap2[ix - 1];
//                  delphi_real phi3 = phimap2[(ix - 1) - lat2];
//                  delphi_real phi4 = phimap2[(ix - 1) + lat1];
//                  delphi_real phi5 = phimap2[(ix - 1) - long2];
//                  delphi_real phi6 = phimap2[(ix - 1) + long1];
//
//                  delphi_real myLastPhi = phimap1[ix - 1] - (phi1 + phi2 + phi3 + phi4 + phi5 + phi6) * sixth;
//
//                  delphi_real numerator   = eps1 * phi1 + eps2 * phi2 + eps3 * phi3 + eps4 * phi4 + eps5 * phi5 + eps6 * phi6;
//                  delphi_real demonimator = eps1 + eps2 + eps3 + eps4 + eps5 + eps6;
//
//                  phimap1[ix - 1] = myLastPhi + (numerator / demonimator) * (1 - (om1));
//              }
//          }
        }
        else // if not Gaussian
        {
            for (n = 0; n < mpi_wricount2a; n++)
            {
                ix = prgiBndyDielecIndex[n];
                temp1 = mpi_phimap2[(ix - 1) -     1] * prgfBndyDielec[n][0] + mpi_phimap2[ix       -     1] * prgfBndyDielec[n][1];
                temp2 = mpi_phimap2[(ix - 1) -  lat2] * prgfBndyDielec[n][2] + mpi_phimap2[(ix - 1) +  lat1] * prgfBndyDielec[n][3];
                temp3 = mpi_phimap2[(ix - 1) - long2] * prgfBndyDielec[n][4] + mpi_phimap2[(ix - 1) + long1] * prgfBndyDielec[n][5];
                mpi_phimap1[ix - 1] += temp1 + temp2 + temp3;
            }           
        }

        /*
         * Now reset boundary values altered in above loops.
         */
        star = (iGrid + 1) / 2;
        fin = (iGrid * (iGrid - 1) - 2) / 2;
        iy = iGrid * (iGrid + 1) / 2 - iGrid + 1;
        
        for (n = 0; n < fin - star + 1; n++) 
        {
            iy += iGrid;
            
            if (mpi_wrphimap1_start <= (iy - 1) && (iy - 1) <= mpi_wrphimap1_end)
                mpi_phimap1[iy - 1] = bndx1[n];
            
            if (mpi_wrphimap1_start <= (iy + ((iGrid + 1) / 2 - 1) - 1) && (iy + ((iGrid + 1) / 2 - 1) - 1) <= mpi_wrphimap1_end)
                mpi_phimap1[iy + ((iGrid + 1) / 2 - 1) - 1] = bndx2[n];
        }   
                
        /*
         * next we add back an adjustment to all the charged grid points due to the charge assigned. the compiler
         * directive just reassures the vector compiler that all is well as far as recurrence is concerned, i.e. it
         * would think there is a recurrence below, where as in fact there is none.
         */
        if (0 != forWhom)
        {
            /*
             * Gaussian based runs are not parallelized yet!
             */
            if (iGaussian != 0)
            {
//              if (fZero < abs(fIonStrength)) // If there is ion, Gaussian
//              {
//                  for (n = 0; n < iCrgedGridEven; n++) 
//                  {
//                      ix = prgiCrgPose[n];
//
//                      delphi_real eps1 = gaussianChargeDielec[n][0];
//                      delphi_real eps2 = gaussianChargeDielec[n][1];
//                      delphi_real eps3 = gaussianChargeDielec[n][2];
//                      delphi_real eps4 = gaussianChargeDielec[n][3];
//                      delphi_real eps5 = gaussianChargeDielec[n][4];
//                      delphi_real eps6 = gaussianChargeDielec[n][5];
//
//                      delphi_real myDensity = gaussianChargeDensity[n];
//                      delphi_real myCharge  = prgfCrgValG[n];
//
//                      delphi_real myExpSolvE = calcExpSolvE(myDensity);
//
//                      delphi_real myNonlinearCorrection = gaussianChargeNonlinear[n];
//
//                      delphi_real numerator   = myCharge * f4Pi * fScale;
//                      delphi_real demonimator = (eps1 + eps2 + eps3 + eps4 + eps5 + eps6) / fEPKT + fDebFct * myExpSolvE;
//
//                      phimap1[ix - 1] = phimap1[ix - 1] + (numerator / demonimator + myNonlinearCorrection) * (1 - (om1));
//                  }
//              } 
//              else // If there is no ion, Gaussian
//              {
//                  for (n = 0; n < iCrgedGridEven; n++) 
//                  {
//                      ix = prgiCrgPose[n];
//
//                      delphi_real eps1 = gaussianChargeDielec[n][0];
//                      delphi_real eps2 = gaussianChargeDielec[n][1];
//                      delphi_real eps3 = gaussianChargeDielec[n][2];
//                      delphi_real eps4 = gaussianChargeDielec[n][3];
//                      delphi_real eps5 = gaussianChargeDielec[n][4];
//                      delphi_real eps6 = gaussianChargeDielec[n][5];
//
//                      delphi_real myCharge = prgfCrgValG[n];
//
//                      delphi_real numerator = myCharge * f4Pi * fScale;
//                      delphi_real demonimator = (eps1 + eps2 + eps3 + eps4 + eps5 + eps6) / fEPKT;
//
//                      phimap1[ix - 1] = phimap1[ix - 1] + (numerator / demonimator) * (1 - (om1));
//                  }
//              }               
            }
            else
            {
                for (n = 0; n < mpi_wricount1a; n++)
                {
                    ix = prgiCrgPose[n];
                    mpi_phimap1[ix - 1] += prgfCrgValA[n];
                }               
            }
        }
        
        /*
         * send back computed phimap1 to master and then recv updated phimap1 w/ corrected bdy values
         */
        if (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2])
        {
                MPI_Gatherv( mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_wrlen1[mpi_rank], mpi_delphi_real, 
                             mpi_phimap1 + mpi_wrphimap1_start, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real, 
                             0, MPI_COMM_WORLD);
                MPI_Scatterv(mpi_phimap1 + mpi_wrphimap1_start, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real, 
                             mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_wrlen1[mpi_rank], mpi_delphi_real, 
                             0, MPI_COMM_WORLD);
        }

    } // ---------- end of splitting work on master and slave processes
    
//    if (1 == flag)
//    {
//        if ( ! (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2]))
//        {
//            if (0 == mpi_rank) /* master process */
//            {            
//                MPI_Gatherv( MPI_IN_PLACE, 0, mpi_delphi_real, 
//                             phimap1.data(), mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real,
//                             0, MPI_COMM_WORLD);
//                MPI_Gatherv( MPI_IN_PLACE, 0, mpi_delphi_real, 
//                             phimap2.data(), mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real,
//                             0, MPI_COMM_WORLD);
//            }
//            else /* slave processes */
//            {
//                MPI_Gatherv( mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_wrlen1[mpi_rank], mpi_delphi_real, 
//                             mpi_phimap1 + mpi_wrphimap1_start, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real,
//                             0, MPI_COMM_WORLD);
//                MPI_Gatherv( mpi_phimap2 + mpi_wrstar2[mpi_rank] - 1, mpi_wrlen2[mpi_rank], mpi_delphi_real, 
//                             mpi_phimap2 + mpi_wrphimap2_start, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real,
//                             0, MPI_COMM_WORLD);
//            }
//        }
//        
//        string strTestFile;
//        
//        if (0 == mpi_rank) strTestFile = "rank0_solver_itrOddPoints.dat";
//        if (1 == mpi_rank) strTestFile = "rank1_solver_itrOddPoints.dat";
//        if (2 == mpi_rank) strTestFile = "rank2_solver_itrOddPoints.dat";
//        if (3 == mpi_rank) strTestFile = "rank3_solver_itrOddPoints.dat";
//        
//        ofstream ofTestStream(strTestFile.c_str());
//        ofTestStream << boolalpha;
//        ofTestStream << fixed << setprecision(7);
//
//        ofTestStream << "flag = " << flag << endl;
//        
//        ix = mpi_wrphimap1_start;
//        for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
//        {
//            ofTestStream << "mpi_phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = mpi_wrphimap2_start;
//        for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
//        {
//            ofTestStream << "mpi_phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ofTestStream.close();
//    }   
}

#endif
