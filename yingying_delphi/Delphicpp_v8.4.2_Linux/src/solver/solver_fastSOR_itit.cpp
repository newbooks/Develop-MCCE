/*
 * solver_fastSOR_itit.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

void CDelphiFastSOR::itit() 
{
    delphi_real rmsch, rmsch2, rmxch, rmxch2;
    int itr, ires;
    delphi_real grden, grdn[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    delphi_integer ix, iy, iz;
    delphi_real maxres = (fRmsc > fMaxc) ? fRmsc : fMaxc;
    maxres = (maxres > fGridConverge) ? maxres : fGridConverge;
    vector<delphi_real> rmsl(nxran, 0.0), rmaxl(nxran, 0.0);
    string strLine60 = " ----------------------------------------------------------------";

    cout << strLine60 << endl;
    if (0.0 < fGridConverge)
        cout << "      " << " rms-change   max change    grid energy    #iterations" << endl;
    else
        cout << "      " << " rms-change   max change       #iterations" << endl;
    cout << strLine60 << endl;


    if (0 == iConvergeFract) 
    {
        iIterateInterval = 10;
        iConvergeFract   = 1;
    }

    if (iIterateInterval > iLinIterateNum) iIterateInterval = iLinIterateNum;

    initOddEvenItr(1); // forWhom = 1

    if (iGaussian != 0)
    {
        //for linear iteration, the pure non-linear part equals to zero 
        gaussianBoundaryNonlinear.assign(iDielecBndyOdd, 0.0);
        gaussianChargeNonlinear.assign(iCrgedGridSum, 0.0);

        //for linear iteration, there is no ion. Thus we won't need a density here
        //if(fabs(fIonStrength)<fZero) gaussianBoundaryDensity.assign(iDielecBndyOdd, 1.0);
    }

    if (debug_solver) 
    {
        cout << "gaussianBoundaryDielec.size= "  << gaussianBoundaryDielec.size()  << endl;
        cout << "gaussianBoundaryDensity.size= " << gaussianBoundaryDensity.size() << endl;
        cout << "gaussianChargeDielec.size= "    << gaussianChargeDielec.size()    << endl;
        cout << "gaussianChargeDensity.size= "   << gaussianChargeDensity.size()   << endl;
        cout << "iDielecBndyEven= " << iDielecBndyEven << endl;
        cout << "iDielecBndyOdd= " << iDielecBndyOdd << endl;
        cout << "iCrgedGridSum= " << iCrgedGridSum << endl;
    }

    itr = 1;
    ires = 0;
    do 
    {
        rmsch = 0.0;
        rmxch = 0.0;

        /*
         * iterate over odd points
         */
        itrOddPoints(1, itr); // forWhom = 1
        
        if (bFixedRelaxParam) 
        {
            int itr2 = 2 * itr - 1;
            om3 = 1.0 / (1.0 - om2 * fSpec * 0.25);
            if (fZero > om1) om3 = 1.0 / (1.0 - om2 * fSpec * 0.5);
            om4 = om3 / om2;
            om2 = om3;
            om1 = 1.0 - om2;

            if (0.0 < fIonStrength) 
            {
                if (1 == itr2 % 2) 
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
                        *it = (*it) * om4;
                } 
                else 
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it)
                        *it = (*it) * om4;
                }
            }

            for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)
                *it = (*it) * om4;

            for (delphi_integer iy = 0; iy < iDielecBndyOdd; iy++)
                for (delphi_integer ix = 0; ix < 6; ix++)
                    prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix] * om4;

            sixth = sixth * om4;
        }
        
        /*
         * Next update phimap2 using the new phimap1
         */
        itrEvenPoints(1, itr); // forWhom = 1     
        
        if (bFixedRelaxParam) 
        {
            int itr2 = 2 * itr;
            om3 = 1.0 / (1.0 - om2 * fSpec * 0.25);
            if (fZero > om1) om3 = 1.0 / (1.0 - om2 * fSpec * 0.5);
            om4 = om3 / om2;
            om2 = om3;
            om1 = 1.0 - om2;

            if (0.0 < fIonStrength) 
            {
                if (1 == itr2 % 2) 
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
                        *it = (*it) * om4;
                } 
                else 
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it)
                        *it = (*it) * om4;
                }
            }

            for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)
                *it = (*it) * om4;

            for (delphi_integer iy = 0; iy < iDielecBndyOdd; iy++)
                for (delphi_integer ix = 0; ix < 6; ix++)
                    prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix] * om4;

            sixth = sixth * om4;
        }

        #ifdef DEBUG_DELPHI_SOLVER_ITIT
        if (1 == itr)
        {
            string strTestFile = "test_itit.dat";
            ofstream ofTestStream(strTestFile.c_str());
            ofTestStream << boolalpha;
            ofTestStream << fixed << setprecision(7);

            ix = 1;
            for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
            {
                ofTestStream << "phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
                ix++;
            }

            ix = 1;
            for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
            {
                ofTestStream << "phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
                ix++;
            }

            ofTestStream.close();
        }
        #endif // DEBUG_DELPHI_SOLVER_ITIT
        
        /*
         * we also save time by only checking convergence every 10 iterations, rather than every single iteration.
         * store phi2 in phi3 to compare against next iteration
         */
        if (iIterateInterval - 1 == itr % iIterateInterval) // itr = 9,19,29,...
        {
            for (ix = 1; ix < iHalfGridNum; ix += iConvergeFract) 
                prgfPhiMap[ix] = phimap2[ix];
        }

        if (0.0 < fGridConverge) 
        {
            grden = 0.0;

            for (ix = 0; ix < iCrgedGridEven; ix++) 
            {
                iy     = prgiCrgPose[ix];
                grden += phimap1[iy - 1] * prgfCrgValG[ix];
            }

            for (ix = iCrgedGridEven; ix < iCrgedGridSum; ix++) 
            {
                iy     = prgiCrgPose[ix];
                grden += phimap2[iy - 1] * prgfCrgValG[ix];
            }

            grdn[itr % 5] = grden / 2.0;
            if (10 < itr) 
            {
                bool igt = true;
                
                for (int i = 0; i < 5; i++)
                    for (int j = 0; j < 5; j++)
                        if (abs(grdn[j] - grdn[i]) > fGridConverge) igt = false;
                
                if (igt) 
                {
                    cout << grdn[0] << " " << grdn[1] << " " << grdn[2] << " " << grdn[3] << " " << grdn[4] << endl;
                    ires = 1;
                }
            }
        }

        if (0 == itr % iIterateInterval || 1 == ires) //----- check to see if accuracy is sufficient
        {            
            delphi_real rnorm2 = 0.0, temp2;

            for (ix = 1; ix < iHalfGridNum; ix += iConvergeFract) 
            {
                temp2 = prgfPhiMap[ix] - phimap2[ix];
                rnorm2 += temp2 * temp2;

                //if(rmxch<abs(temp2)) flag_itr=ix; //Lin Li: to see which grid is the key grid
                rmxch = max(rmxch, abs(temp2));
            }

            rmsch = sqrt((delphi_real) iConvergeFract * rnorm2 / ((iGrid - 2) * (iGrid - 2) * (iGrid - 2)));
            //rnormch = sqrt(rnorm2);
            rmsch2 = rmsch;
            rmxch2 = rmxch;

            if (0.0 < fGridConverge)
                cout << "        " << scientific << rmsch2 << "  " << rmxch2 << "  " << grden << "  at  " 
                     << setw(5) << left << itr << " iterations\n";
            else
                cout << "        " << scientific << rmsch2 << "  " << rmxch2 << "  at  " << setw(5) 
                     << left << itr << " iterations\n";

            if (fRmsc > rmsch || fMaxc > rmxch) ires = 1;

            if (bLogGraph) 
            {
                int ibin;
                for (int j = itr - 9; j <= itr; j++) 
                {
                    ibin = (j - 1) * (60 - 1) / (iLinIterateNum - 1) + 1;
                    rmsl[ibin - 1] = rmsch;
                    rmaxl[ibin - 1] = rmxch;
                }
            }

//            //if (9 == itr)
//            {
//                string strTestFile = "rank1_solver_itit.dat";
//                ofstream ofTestStream(strTestFile.c_str());
//                ofTestStream << boolalpha;
//                ofTestStream << fixed << setprecision(7);
//
//                ix = 0;
//                for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
//                {
//                    ofTestStream << "mpi_phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//                    ix++;
//                }
//
//                ix = 0;
//                for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
//                {
//                    ofTestStream << "mpi_phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//                    ix++;
//                }
//
//                ofTestStream.close();
//                
//                return;
//            }
        } 
        
        itr++;
        
        /*
         * check to see if accuracy is sufficient
         */
        if (1.0e-7 > maxres) 
        {
            if (iLinIterateNum >= itr && 0 == ires)
                continue;
            else
                break;
        } 
        else 
        {
            if ((iLinIterateNum >= itr || bAutoConverge) && (0 == ires))
                continue;
            else
                break;
        }
        
    } while (true);  
    
    //Argo: Printing the last rms-change of the iteration
    //Needed for GAUSSIAN runs since many were found to diverge
    cout << strLine60 << endl;
    cout << infoString << "Iteration ended with final rms-change of " << scientific << rmsch2 << endl;
    
    if (rmsch2 > fMaxc) 
    {
        //cout << " WARNING !!! Run probably DIVERGED" << endl;
        //cout << " WARNING !!! Final rms-change is " << rmsch2 << " greater than " << fMaxc << endl;
        //cout << " WARNING !!! Try increasing the scale or lower the maxc value or both" << endl;
        stringstream str_Diverged;
        str_Diverged << "Run probably DIVERGED! Final rms-change is " << rmsch2 << " greater than " << fMaxc << endl;
        str_Diverged << "             Try increasing the scale or lower the maxc value or both" << endl;
        CDiverged warn_div(str_Diverged);
    } 
    else 
    {
        cout << infoString << "Run converged in the provided limits" << endl;
    }

    postItr(rmaxl, rmsl);

    /*
     * code phimap corner, for use in transference from irises to convex and via versa
     */
    {
        delphi_real ap1, ap2, ap3, ap4;
        ap1 = prgfPhiMap[0];
        ap2 = ap1 * 10000.0;
        ap3 = (int) ap2;
        if (0 < ap3)
            ap4 = (ap3 + 0.8) / 10000.0;
        else
            ap4 = (ap3 - 0.8) / 10000.0;
        prgfPhiMap[0] = ap4;
    }

    #ifdef DEBUG_DELPHI_SOLVER_ITIT
    {
        string strTestFile = "test_itit.dat";
        ofstream ofTestStream(strTestFile.c_str());
        ofTestStream << boolalpha;
        ofTestStream << fixed << setprecision(7);

        const delphi_real *** phimap = pdc->getKey_constPtr<delphi_real>("phimap",iGrid,iGrid,iGrid); // const pointer to 3D phimap
        for (iz = 0; iz < iGrid; iz++)
            for (iy = 0; iy < iGrid; iy++)
                for (ix = 0; ix < iGrid; ix++)
                    ofTestStream << "phimap[" << setw(6) << right << iz+1 << "," << setw(6) 
                                 << right << iy+1 << "," << setw(6) << right << ix+1 << "] = " 
                                 << setw(12) << right << phimap[iz][iy][ix] << endl;

        ofTestStream.close();
    }
    #endif // DEBUG_DELPHI_SOLVER_ITIT

}


#ifdef PARALLEL_MPI

/*
 * subroutine mpi_itit() run on all processes
 */
void CDelphiFastSOR::mpi_itit()
{
    delphi_real rmsch, rmsch2, rmxch, rmxch2;
    int itr, ires;
    delphi_real grden, grdn[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    delphi_integer ix, iy, iz;
    delphi_real maxres = (fRmsc > fMaxc) ? fRmsc : fMaxc;
    maxres = (maxres > fGridConverge) ? maxres : fGridConverge;
    vector<delphi_real> rmsl(nxran, 0.0), rmaxl(nxran, 0.0);
    string strLine60 = " ----------------------------------------------------------------";

    vector<delphi_real> realbuff;
    delphi_real *mpi_realbuff;
    
    if (0 == mpi_rank)
    {
        cout << strLine60 << endl;
        if (0.0 < fGridConverge)
            cout << "      " << " rms-change   max change    grid energy    #iterations" << endl;
        else
            cout << "      " << " rms-change   max change       #iterations" << endl;
        cout << strLine60 << endl;
    }

    if (0 == iConvergeFract)
    {
        iIterateInterval = 10;
        iConvergeFract   = 1;
    }

    if (iIterateInterval > iLinIterateNum) iIterateInterval = iLinIterateNum;

    mpi_initOddEvenItr(1); // forWhom = 1
    
    /*
     *  Gaussian based runs are not parallelized yet!
     */
//    if (iGaussian != 0)
//    {
//        gaussianBoundaryNonlinear.assign(iDielecBndyOdd, 0.0);
//        gaussianChargeNonlinear.assign(iCrgedGridSum, 0.0);  
//    }    
    
    /*
     * iteration starts
     */
    itr  = 1;
    ires = 0;
    do
    {
        rmsch = 0.0;
        rmxch = 0.0;

        /*
         * iterate over odd points
         */
        mpi_itrOddPoints(1, itr); // forWhom = 1    
        
        if (bFixedRelaxParam)
        {
            int itr2 = 2 * itr - 1;
            om3 = 1.0 / (1.0 - om2 * fSpec * 0.25);
            if (fZero > om1) om3 = 1.0 / (1.0 - om2 * fSpec * 0.5);
            om4 = om3 / om2;
            om2 = om3;
            om1 = 1.0 - om2;

            if (0.0 < fIonStrength)
            {
                if (1 == itr2 % 2)
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
                        *it = (*it) * om4;
                }
                else
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it)
                        *it = (*it) * om4;
                }
            }

            for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)
                *it = (*it) * om4;

            for (delphi_integer iy = 0; iy < prgfBndyDielec.size(); iy++)
                for (delphi_integer ix = 0; ix < 6; ix++)
                    prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix] * om4;

            sixth = sixth * om4;
        }
        
        /*
         * Synchronization before getting into the loop
         */
        MPI_Barrier (MPI_COMM_WORLD);

        /*
         * Next update phimap2 using the new phimap1
         */
        mpi_itrEvenPoints(1, itr); // forWhom = 1
        
        if (bFixedRelaxParam)
        {
            int itr2 = 2 * itr;
            om3 = 1.0 / (1.0 - om2 * fSpec * 0.25);
            if (fZero > om1) om3 = 1.0 / (1.0 - om2 * fSpec * 0.5);
            om4 = om3 / om2;
            om2 = om3;
            om1 = 1.0 - om2;

            if (0.0 < fIonStrength)
            {
                if (1 == itr2 % 2)
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
                        *it = (*it) * om4;
                }
                else
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it)
                        *it = (*it) * om4;
                }
            }

            for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)
                *it = (*it) * om4;

            for (delphi_integer iy = 0; iy < prgfBndyDielec.size(); iy++)
                for (delphi_integer ix = 0; ix < 6; ix++)
                    prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix] * om4;

            sixth = sixth * om4;
        }  
        
        /*
         * we also save time by only checking convergence every 10 iterations, rather than every single iteration.
         * store phi2 in phi3 to compare against next iteration
         */
        if ( (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2]) && 0 == mpi_rank) // if periodic bdy, only master process works
        {
            if (iIterateInterval - 1 == itr % iIterateInterval) // itr = 9,19,29,...
            {
                for (ix = 1; ix < iHalfGridNum; ix += iConvergeFract)
                    prgfPhiMap[ix] = phimap2[ix];
            }

            if (0.0 < fGridConverge)
            {
                grden = 0.0;

                for (ix = 0; ix < iCrgedGridEven; ix++)
                {
                    iy     = prgiCrgPose[ix];
                    grden += phimap1[iy - 1] * prgfCrgValG[ix];
                }

                for (ix = iCrgedGridEven; ix < iCrgedGridSum; ix++)
                {
                    iy     = prgiCrgPose[ix];
                    grden += phimap2[iy - 1] * prgfCrgValG[ix];
                }

                grdn[itr % 5] = grden / 2.0; /* modified to save on grdn dimension */
                if (10 < itr)
                {
                    bool igt = true;
                    for (int i = 0; i < 5; i++)
                        for (int j = 0; j < 5; j++)
                            if (abs(grdn[j] - grdn[i]) > fGridConverge)
                                igt = false;
                    if (igt)
                    {
                        cout << grdn[0] << " " << grdn[1] << " " << grdn[2] << " " << grdn[3] << " " << grdn[4] << endl;
                        ires = 1;
                    }
                }
            }

            if (0 == itr % iIterateInterval || 1 == ires) //----- check to see if accuracy is sufficient
            {
                delphi_real rnorm2 = 0.0, temp2;

                for (ix = 1; ix < iHalfGridNum; ix += iConvergeFract)
                {
                    temp2   = prgfPhiMap[ix] - phimap2[ix];
                    rnorm2 += temp2 * temp2;
                    rmxch   = max(rmxch, abs(temp2));
                }

                rmsch  = sqrt( (delphi_real) iConvergeFract * rnorm2 / ((iGrid - 2) * (iGrid - 2) * (iGrid - 2)) );
                rmsch2 = rmsch;
                rmxch2 = rmxch;

                if (0.0 < fGridConverge)
                    cout << "        " << scientific << rmsch2 << "  " << rmxch2 << "  " << grden << "  at  " 
                         << setw(5) << left << itr << " iterations\n";
                else
                    cout << "        " << scientific << rmsch2 << "  " << rmxch2 << "  at  " 
                         << setw(5) << left << itr << " iterations\n";

                if (fRmsc > rmsch || fMaxc > rmxch)
                    ires = 1;

                if (bLogGraph)
                {
                    int ibin;
                    for (int j = itr - 9; j <= itr; j++)
                    {
                        ibin            = (j - 1) * (60 - 1) / (iLinIterateNum - 1) + 1;
                        rmsl[ibin - 1]  = rmsch;
                        rmaxl[ibin - 1] = rmxch;
                    }
                }
            }
        }
        else // if not periodic bdy, all processes work together...
        {
            delphi_real mpi_sumall;

            if (0 == mpi_rank) /* master process */
            {
                if (0.0 < fGridConverge)
                {
                    grden = 0.0;
                    MPI_Allreduce(&grden, &mpi_sumall, 1, mpi_delphi_real, MPI_SUM, MPI_COMM_WORLD);
                    grden = mpi_sumall;

                    grdn[itr % 5] = grden / 2.0;
                    if (10 < itr)
                    {
                        bool igt = true;
                        for (int i = 0; i < 5; i++)
                            for (int j = 0; j < 5; j++)
                                if (abs(grdn[j] - grdn[i]) > fGridConverge)
                                    igt = false;
                        if (igt)
                        {
                            cout << grdn[0] << " " << grdn[1] << " " << grdn[2] << " " << grdn[3] << " " << grdn[4] << endl;
                            ires = 1;
                        }

                        MPI_Bcast(&ires, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
                    }
                }

                /*
                 * check to see if accuracy is sufficient
                 */
                if (0 == itr % iIterateInterval || 1 == ires)
                {                    
                    delphi_real rnorm2 = 0.0, temp2;

                    MPI_Allreduce(&rnorm2, &mpi_sumall, 1, mpi_delphi_real, MPI_SUM, MPI_COMM_WORLD);
                    rnorm2 = mpi_sumall;

                    MPI_Allreduce(&rmxch, &mpi_sumall, 1, mpi_delphi_real, MPI_MAX, MPI_COMM_WORLD);
                    rmxch = mpi_sumall;

                    rmsch  = sqrt( (delphi_real) iConvergeFract * rnorm2 / ((iGrid - 2) * (iGrid - 2) * (iGrid - 2)) );
                    rmsch2 = rmsch;
                    rmxch2 = rmxch;

                    if (0.0 < fGridConverge)
                        cout << "        " << scientific << rmsch2 << "  " << rmxch2 << "  " << grden << "  at  " 
                             << setw(5) << left << itr << " iterations\n";
                    else
                        cout << "        " << scientific << rmsch2 << "  " << rmxch2 << "  at  " 
                             << setw(5) << left << itr << " iterations\n";

                    if (fRmsc > rmsch || fMaxc > rmxch) ires = 1;

                    if (bLogGraph)
                    {
                        int ibin;
                        for (int j = itr - 9; j <= itr; j++)
                        {
                            ibin            = (j - 1) * (60 - 1) / (iLinIterateNum - 1) + 1;
                            rmsl[ibin - 1]  = rmsch;
                            rmaxl[ibin - 1] = rmxch;
                        }
                    }
                    
//                    //if (9 == itr)
//                    {
//                        if ( ! (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2]))
//                        {
//                            if (0 == mpi_rank) /* master process */
//                            {            
//                                MPI_Gatherv( MPI_IN_PLACE, 0, mpi_delphi_real, 
//                                             phimap1.data(), mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real,
//                                             0, MPI_COMM_WORLD);
//                                MPI_Gatherv( MPI_IN_PLACE, 0, mpi_delphi_real, 
//                                             phimap2.data(), mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real,
//                                             0, MPI_COMM_WORLD);
//                            }
//                            else /* slave processes */
//                            {
//                                MPI_Gatherv( mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_wrlen1[mpi_rank], mpi_delphi_real, 
//                                             mpi_phimap1 + mpi_wrphimap1_start, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real,
//                                             0, MPI_COMM_WORLD);
//                                MPI_Gatherv( mpi_phimap2 + mpi_wrstar2[mpi_rank] - 1, mpi_wrlen2[mpi_rank], mpi_delphi_real, 
//                                             mpi_phimap2 + mpi_wrphimap2_start, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real,
//                                             0, MPI_COMM_WORLD);
//                            }
//                        }
//                        
//                        string strTestFile;
//                        
//                        if (0 == mpi_rank) strTestFile = "rank0_solver_itit.dat";
//                        if (1 == mpi_rank) strTestFile = "rank1_solver_itit.dat";
//                        if (2 == mpi_rank) strTestFile = "rank2_solver_itit.dat";
//                        if (3 == mpi_rank) strTestFile = "rank3_solver_itit.dat";
//                        
//                        ofstream ofTestStream(strTestFile.c_str());
//                        ofTestStream << boolalpha;
//                        ofTestStream << fixed << setprecision(7);
//
//                        ix = mpi_wrphimap1_start;
//                        for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
//                        {
//                            ofTestStream << "mpi_phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//                            ix++;
//                        }
//
//                        ix = mpi_wrphimap2_start;
//                        for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
//                        {
//                            ofTestStream << "mpi_phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//                            ix++;
//                        }
//
//                        ofTestStream.close();
//                        
//                        return;
//                    }
                }
            }
            else /* slave processes */
            {
                /*
                 * we also save time by only checking convergence every 10 iterations, rather than every single iteration.
                 * store phi2 in phi3 to compare against next iteration
                 */
                if (iIterateInterval - 1 == itr % iIterateInterval) // itr = 9,19,29,...
                {
                    realbuff.assign(phimap2.begin(), phimap2.end());
                    mpi_realbuff = realbuff.data() - mpi_wrphimap2_start;
                }

                if (0.0 < fGridConverge)
                {
                    grden = 0.0;

                    for (ix = 0; ix < mpi_wricount1a; ix++)
                    {
                        iy     = prgiCrgPose[ix];
                        grden += phimap1[iy - 1] * prgfCrgValG[ix];
                    }

                    for (ix = mpi_wricount1a; ix < mpi_wricount1b; ix++)
                    {
                        iy     = prgiCrgPose[ix];
                        grden += phimap2[iy - 1] * prgfCrgValG[ix];
                    }

                    MPI_Allreduce(&grden, &mpi_sumall, 1, mpi_delphi_real, MPI_SUM, MPI_COMM_WORLD);
                    grden = mpi_sumall;

                    if (10 < itr)
                        MPI_Bcast(&ires, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
                }

                /*
                 * check to see if accuracy is sufficient
                 */
                if (0 == itr % iIterateInterval || 1 == ires)
                {                    
                    delphi_real rnorm2 = 0.0, temp2;

                    for (ix = mpi_wrstar2[mpi_rank] - 1; ix < mpi_wrfinl2[mpi_rank]; ix += iConvergeFract)
                    {
                        temp2   = mpi_realbuff[ix] - mpi_phimap2[ix];
                        rnorm2 += temp2 * temp2;
                        rmxch   = max(rmxch, abs(temp2));
                    }

                    MPI_Allreduce(&rnorm2, &mpi_sumall, 1, mpi_delphi_real, MPI_SUM, MPI_COMM_WORLD);
                    rnorm2 = mpi_sumall;

                    MPI_Allreduce(&rmxch, &mpi_sumall, 1, mpi_delphi_real, MPI_MAX, MPI_COMM_WORLD);
                    rmxch = mpi_sumall;

                    rmsch  = sqrt( (delphi_real) iConvergeFract * rnorm2 / ((iGrid - 2) * (iGrid - 2) * (iGrid - 2)) );
                    rmsch2 = rmsch;
                    rmxch2 = rmxch;

                    if (fRmsc > rmsch || fMaxc > rmxch) ires = 1;
                    
//                    //if (9 == itr)
//                    {
//                        if ( ! (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2]))
//                        {
//                            if (0 == mpi_rank) /* master process */
//                            {            
//                                MPI_Gatherv( MPI_IN_PLACE, 0, mpi_delphi_real, 
//                                             phimap1.data(), mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real,
//                                             0, MPI_COMM_WORLD);
//                                MPI_Gatherv( MPI_IN_PLACE, 0, mpi_delphi_real, 
//                                             phimap2.data(), mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real,
//                                             0, MPI_COMM_WORLD);
//                            }
//                            else /* slave processes */
//                            {
//                                MPI_Gatherv( mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_wrlen1[mpi_rank], mpi_delphi_real, 
//                                             mpi_phimap1 + mpi_wrphimap1_start, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real,
//                                             0, MPI_COMM_WORLD);
//                                MPI_Gatherv( mpi_phimap2 + mpi_wrstar2[mpi_rank] - 1, mpi_wrlen2[mpi_rank], mpi_delphi_real, 
//                                             mpi_phimap2 + mpi_wrphimap2_start, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real,
//                                             0, MPI_COMM_WORLD);
//                            }
//                        }
//                        
//                        string strTestFile;
//                        
//                        if (0 == mpi_rank) strTestFile = "rank0_solver_itit.dat";
//                        if (1 == mpi_rank) strTestFile = "rank1_solver_itit.dat";
//                        if (2 == mpi_rank) strTestFile = "rank2_solver_itit.dat";
//                        if (3 == mpi_rank) strTestFile = "rank3_solver_itit.dat";
//                        
//                        ofstream ofTestStream(strTestFile.c_str());
//                        ofTestStream << boolalpha;
//                        ofTestStream << fixed << setprecision(7);
//
//                        ix = mpi_wrphimap1_start;
//                        for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
//                        {
//                            ofTestStream << "mpi_phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//                            ix++;
//                        }
//
//                        ix = mpi_wrphimap2_start;
//                        for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
//                        {
//                            ofTestStream << "mpi_phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//                            ix++;
//                        }
//
//                        ofTestStream.close();
//                        
//                        return;
//                    }
                }
            } // end of splitting between master and slave processes
        }
        
        itr++;
        
        /*
         * check to see if accuracy is sufficient
         */
        {
            int mpi_command;

            if (0 == mpi_rank) /* master process */
            {
                if (1.0e-7 > maxres)
                {
                    if (iLinIterateNum >= itr && 0 == ires)
                    {
                        mpi_command  = 1;
                        MPI_Bcast(&mpi_command, 1, MPI_INT, 0, MPI_COMM_WORLD);
                        continue;
                    }
                    else
                    {
                        mpi_command = 0;
                        MPI_Bcast(&mpi_command, 1, MPI_INT, 0, MPI_COMM_WORLD);
                        break;
                    }
                }
                else
                {
                    if ((iLinIterateNum >= itr || bAutoConverge) && (0 == ires))
                    {
                        mpi_command = 1;
                        MPI_Bcast(&mpi_command, 1, MPI_INT, 0, MPI_COMM_WORLD);
                        continue;
                    }
                    else
                    {
                        mpi_command = 0;
                        MPI_Bcast(&mpi_command, 1, MPI_INT, 0, MPI_COMM_WORLD);
                        break;
                    }
                }
            }
            else /* slave processes */
            {
                MPI_Bcast(&mpi_command, 1, MPI_INT, 0, MPI_COMM_WORLD);
                if (0 == mpi_command) break;
            }
        }
        
    } while (true);
    
    if ( ! (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2]))
    {
        if (0 == mpi_rank) /* master process */
        {            
            MPI_Gatherv( MPI_IN_PLACE, 0, mpi_delphi_real, 
                         phimap1.data(), mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real,
                         0, MPI_COMM_WORLD);
            MPI_Gatherv( MPI_IN_PLACE, 0, mpi_delphi_real, 
                         phimap2.data(), mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real,
                         0, MPI_COMM_WORLD);
        }
        else /* slave processes */
        {
            MPI_Gatherv( mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_wrlen1[mpi_rank], mpi_delphi_real, 
                         mpi_phimap1 + mpi_wrphimap1_start, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real,
                         0, MPI_COMM_WORLD);
            MPI_Gatherv( mpi_phimap2 + mpi_wrstar2[mpi_rank] - 1, mpi_wrlen2[mpi_rank], mpi_delphi_real, 
                         mpi_phimap2 + mpi_wrphimap2_start, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real,
                         0, MPI_COMM_WORLD);
        }
    }
    
    if (1 < mpi_num_workers)
    {
        MPI_Win_free(&mpi_towin1);
        MPI_Win_free(&mpi_towin2);
        MPI_Group_free(&mpi_wholegroup);
        MPI_Group_free(&mpi_postgroup);
        MPI_Group_free(&mpi_startgroup);
    }
    
    if (0 == mpi_rank)
    {
        /*
         * Argo: printing the last rms-change of the iteration. Needed for GAUSSIAN runs since many were found to diverge
         */
        cout << strLine60 << endl;
        cout << infoString << "Iteration ended with final rms-change of " << scientific << rmsch2 << endl;

        if (rmsch2 > fMaxc)
        {
            stringstream str_Diverged;
            str_Diverged << "Run probably DIVERGED! Final rms-change is " << rmsch2 << " greater than " << fMaxc << endl;
            str_Diverged << "             Try increasing the scale or lower the maxc value or both" << endl;
            CDiverged warn_div(str_Diverged);
        }
        else
        {
            cout << infoString << "Run converged in the provided limits" << endl;
        }

        postItr(rmaxl, rmsl);

        /*
         * code phimap corner, for use in transference from irises to convex and via versa
         */
        {
            delphi_real ap1, ap2, ap3, ap4;
            ap1 = prgfPhiMap[0];
            ap2 = ap1 * 10000.0;
            ap3 = (int) ap2;
            if (0 < ap3)
                ap4 = (ap3 + 0.8) / 10000.0;
            else
                ap4 = (ap3 - 0.8) / 10000.0;
            prgfPhiMap[0] = ap4;
        }
    }
}

#endif
