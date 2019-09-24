/*
 * solver_fastSOR_initOddEvenItr.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

void CDelphiFastSOR::initOddEvenItr(const int& forWhom) 
{
    int ix, iy, iz;
    delphi_integer ip, iq;
    delphi_integer iadd1, iadd2, itemp1, itemp2, star, fin;

    /*
     * setup periodic boundaries
     */
    if (rgbPeriodicBndy[0]) 
    {
        for (iz = 1; iz < iGrid - 1; iz++) 
        {
            iadd1 = iz * iGrid * iGrid;
            for (iy = 1; iy < iGrid - 1; iy++) 
            {
                iadd2 = (iadd1 + iy * iGrid + 2) / 2;
                ibndx.push_back(iadd2);
            }
        }

        idif1x = (iGrid - 2) / 2;
        idif2x = idif1x + 1;
        inc1xa = 1;
        inc1xb = 0;
        inc2xa = 0;
        inc2xb = 1;
    }

    if (rgbPeriodicBndy[1]) 
    {
        for (iz = 1; iz < iGrid - 1; iz++) 
        {
            iadd1 = iz * iGrid * iGrid;
            for (ix = 1; ix < iGrid - 1; ix++) 
            {
                iadd2 = (iadd1 + ix + 2) / 2;
                ibndy.push_back(iadd2);
            }
        }

        idif1y = iGrid * (iGrid - 2) / 2;
        idif2y = idif1y + 1;
        inc1ya = iGrid / 2 + 1;
        inc1yb = inc1ya - 1;
        inc2ya = inc1yb;
        inc2yb = inc1ya;
    }

    if (rgbPeriodicBndy[2]) 
    {
        for (ix = 1; ix < iGrid - 1; ix++) 
        {
            iadd1 = ix + 2;
            for (iy = 1; iy < iGrid - 1; iy++) 
            {
                iadd2 = (iadd1 + iy * iGrid) / 2;
                ibndz.push_back(iadd2);
            }
        }

        idif1z = iGrid * iGrid * (iGrid - 2) / 2;
        idif2z = idif1z + 1;
        inc1za = iGrid * iGrid / 2 + 1;
        inc1zb = inc1za;
        inc2za = inc1zb;
        inc2zb = inc1za;
    }

    /*
     * setup start and stop vectors
     */
    sta1.assign(iGrid, 0); sta1[1] = (iGrid * iGrid + iGrid + 4) / 2;
    fi1.assign(iGrid,  0);  fi1[1] = iGrid * iGrid - (iGrid + 1) / 2;
    sta2.assign(iGrid, 0); sta2[1] = sta1[1] - 1;
    fi2.assign(iGrid,  0);  fi2[1] = fi1[1];

    itemp1  = iGrid + 2;
    itemp2  = iGrid * iGrid - iGrid - 2;
    for (delphi_integer i = 2; i < iGrid - 1; i++) 
    {
        sta1[i] =  fi1[i - 1] + itemp1;
        sta2[i] =  fi2[i - 1] + itemp1;
        fi1[i]  = sta1[i - 1] + itemp2;
        fi2[i]  = sta2[i - 1] + itemp2;
    }

    lat1 = (iGrid - 1) / 2; long1 = (iGrid * iGrid - 1) / 2;
    lat2 = (iGrid + 1) / 2; long2 = (iGrid * iGrid + 1) / 2;

    /*
     * split prgfPhiMap to odd and even vectors
     */
    for (ip = 0; ip < iHalfGridNum; ip++) 
    {
        iq = ip * 2;
        if (prgfPhiMap.size() > iq)     phimap1[ip] = prgfPhiMap[iq];
        if (prgfPhiMap.size() > iq + 1) phimap2[ip] = prgfPhiMap[iq + 1];
    }

    /*
     * setup vectors for restoring x boundary values
     */
    star = (iGrid + 1) / 2;
    iq   = iGrid * (iGrid + 1) / 2 - iGrid + 1;
    fin  = (iGrid * (iGrid - 1) - 2) / 2;
    for (ip = 0; ip < fin - star + 1; ip++) 
    {
        iq += iGrid;
        bndx1.push_back(phimap1[iq - 1]);
        bndx2.push_back(phimap1[iq + ((iGrid + 1) / 2 - 1) - 1]);
    }

    star = (iGrid + 2) / 2;
    iq   = iGrid * (iGrid + 2) / 2 - iGrid + 1;
    fin  = (iGrid * (iGrid - 1) - 1) / 2;
    for (ip = 0; ip < fin - star + 1; ip++) 
    {
        iq += iGrid;
        bndx3.push_back(phimap2[iq - 1]);
        bndx4.push_back(phimap2[iq + ((iGrid + 1) / 2 - 1) - 1]);
    }

    if (0 == forWhom || (1 == forWhom && bFixedRelaxParam)) 
    {
        om2 = 1.0;
        sixth = fSixth;
    } 
    else 
    {
        om2 = 2.0 / (1.0 + sqrt(1.0 - fSpec));

        for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
            *it = (*it) * om2;

        for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it)
            *it = (*it) * om2;

        for (vector<delphi_real>::iterator it = prgfCrgValA.begin();  it != prgfCrgValA.end(); ++it)
            *it = (*it) * om2;

        for (iy = 0; iy < iDielecBndyOdd; iy++)
            for (ix = 0; ix < 6; ix++) 
                prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix] * om2;

        sixth = fSixth * om2;
    }

    om1 = 1.0 - om2;
     
    qmap1.assign(iHalfGridNum, 0.0);
    qmap2.assign(iHalfGridNum, 0.0);
    
    /*
     * test code
     */
//    {
//        string strTestFile = "rank1_solver_initOddEvenItr.dat";
//        
//        ofstream ofTestStream(strTestFile.c_str());
//        ofTestStream << boolalpha;
//        ofTestStream << fixed << setprecision(7);
//
//        ix = 0;
//        for (vector<delphi_integer>::iterator it = ibndx.begin(); it != ibndx.end(); ++it)
//        {
//            ofTestStream << "ibndx[" << setw(6) << right << ix << "] = " << setw(8) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_integer>::iterator it = ibndy.begin(); it != ibndy.end(); ++it)
//        {
//            ofTestStream << "ibndy[" << setw(6) << right << ix << "] = " << setw(8) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_integer>::iterator it = ibndz.begin(); it != ibndz.end(); ++it)
//        {
//            ofTestStream << "ibndz[" << setw(6) << right << ix << "] = " << setw(8) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = bndx1.begin(); it != bndx1.end(); ++it)
//        {
//            ofTestStream << "bndx1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = bndx2.begin(); it != bndx2.end(); ++it)
//        {
//            ofTestStream << "bndx2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = bndx3.begin(); it != bndx3.end(); ++it)
//        {
//            ofTestStream << "bndx3[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = bndx4.begin(); it != bndx4.end(); ++it)
//        {
//            ofTestStream << "bndx4[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); it++)
//        {
//            ofTestStream << "mpi_phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); it++)
//        {
//            ofTestStream << "mpi_phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); it++)
//        {
//            ofTestStream << "mpi_sf1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }        
//        
//        ix = 0;
//        for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); it++)
//        {
//            ofTestStream << "mpi_sf2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }  
//        
//        ofTestStream.close();
//    } // end of test code
}

#ifdef PARALLEL_MPI

/*
 * subroutine mpi_initOddEvenItr(forWhom)
 */
void CDelphiFastSOR::mpi_initOddEvenItr(const int& forWhom)
{
    int ix, iy, iz, iw;
    delphi_integer ip, iq;
    delphi_integer itemp1, itemp2, iadd1, iadd2, ihgd2, star, fin;

    int mpi_size;
    MPI_Info mpi_info;
    
    if (0 == mpi_rank) /* master process performs some initialization */
    {
        /*
         * setup periodic boundaries
         */
        if (rgbPeriodicBndy[0])
        {
            for (iz = 1; iz < iGrid - 1; iz++)
            {
                iadd1 = iz * iGrid * iGrid;
                for (iy = 1; iy < iGrid - 1; iy++)
                {
                    iadd2 = (iadd1 + iy * iGrid + 2) / 2;
                    ibndx.push_back(iadd2);
                }
            }

            idif1x = (iGrid - 2) / 2;
            idif2x = idif1x + 1;
            inc1xa = 1;
            inc1xb = 0;
            inc2xa = 0;
            inc2xb = 1;
        }

        if (rgbPeriodicBndy[1])
        {
            for (iz = 1; iz < iGrid - 1; iz++)
            {
                iadd1 = iz * iGrid * iGrid;
                for (ix = 1; ix < iGrid - 1; ix++)
                {
                    iadd2 = (iadd1 + ix + 2) / 2;
                    ibndy.push_back(iadd2);
                }
            }

            idif1y = iGrid * (iGrid - 2) / 2;
            idif2y = idif1y + 1;
            inc1ya = iGrid / 2 + 1;
            inc1yb = inc1ya - 1;
            inc2ya = inc1yb;
            inc2yb = inc1ya;
        }

        if (rgbPeriodicBndy[2])
        {
            for (ix = 1; ix < iGrid - 1; ix++)
            {
                iadd1 = ix + 2;
                for (iy = 1; iy < iGrid - 1; iy++)
                {
                    iadd2 = (iadd1 + iy * iGrid) / 2;
                    ibndz.push_back(iadd2);
                }
            }

            idif1z = iGrid * iGrid * (iGrid - 2) / 2;
            idif2z = idif1z + 1;
            inc1za = iGrid * iGrid / 2 + 1;
            inc1zb = inc1za;
            inc2za = inc1zb;
            inc2zb = inc1za;
        }

        /*
         * setup start and stop vectors are moved in below to be performed on all processes
         */

        /*
         * split prgfPhiMap to odd and even vectors
         */
        for (ip = 0; ip < iHalfGridNum; ip++)
        {
            iq = ip * 2;
            if (prgfPhiMap.size() > iq)     phimap1[ip] = prgfPhiMap[iq];
            if (prgfPhiMap.size() > iq + 1) phimap2[ip] = prgfPhiMap[iq + 1];
        }

        /*
         * setup vectors for restoring x boundary values
         */
        star = (iGrid + 1) / 2;
        iq   = iGrid * (iGrid + 1) / 2 - iGrid + 1;
        fin  = (iGrid * (iGrid - 1) - 2) / 2;
        for (ip = 0; ip < fin - star + 1; ip++)
        {
            iq += iGrid;
            bndx1.push_back(phimap1[iq - 1]);
            bndx2.push_back(phimap1[iq + ((iGrid + 1) / 2 - 1) - 1]);
        }

        star = (iGrid + 2) / 2;
        iq   = iGrid * (iGrid + 2) / 2 - iGrid + 1;
        fin  = (iGrid * (iGrid - 1) - 1) / 2;
        for (ip = 0; ip < fin - star + 1; ip++)
        {
            iq += iGrid;
            bndx3.push_back(phimap2[iq - 1]);
            bndx4.push_back(phimap2[iq + ((iGrid + 1) / 2 - 1) - 1]);
        }

        if (0 == forWhom || (1 == forWhom && bFixedRelaxParam))
        {
            om2 = 1.0;
            sixth = fSixth;
        }
        else
        {
            om2 = 2.0 / (1.0 + sqrt(1.0 - fSpec));

            for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it) 
                    *it = (*it) * om2;
            
            for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it) 
                    *it = (*it) * om2;
            for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)   
                    *it = (*it) * om2;

            for (iy = 0; iy < iDielecBndyOdd; iy++)
                for (ix = 0; ix < 6; ix++)
                    prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix] * om2;

            sixth = fSixth * om2;
        }

    }
    else  /* slave processes perform some initialization */
    {
        if (0 == forWhom || (1 == forWhom && bFixedRelaxParam))
        {
            om2 = 1.0;
            sixth = fSixth;
        } 
        else
        {
            om2 = 2.0 / (1.0 + sqrt(1.0 - fSpec));
            sixth = fSixth * om2;
        }
    }

    om1 = 1.0 - om2;

    /*
     * STEP 1: setup start and stop vectors
     *         indices of sta1-2, fi1-2 actually start with index 1 in C++!
     *         sta1-2 and fi1-2 are duplicated on all processes in order to avoid sending/receiving costs
     */
    sta1.assign(iGrid, 0); sta1[1] = (iGrid * iGrid + iGrid + 4) / 2;
     fi1.assign(iGrid, 0);  fi1[1] = iGrid * iGrid - (iGrid + 1) / 2;
    sta2.assign(iGrid, 0); sta2[1] = sta1[1] - 1;
     fi2.assign(iGrid, 0);  fi2[1] = fi1[1];

    itemp1 = iGrid + 2;
    itemp2 = iGrid * iGrid - iGrid - 2;
    for (delphi_integer i = 2; i < iGrid - 1; i++)
    {
        sta1[i] =  fi1[i - 1] + itemp1;
         fi1[i] = sta1[i - 1] + itemp2;
        sta2[i] =  fi2[i - 1] + itemp1;
         fi2[i] = sta2[i - 1] + itemp2;
    }

    lat1 = (iGrid - 1) / 2; long1 = (iGrid * iGrid - 1) / 2;
    lat2 = (iGrid + 1) / 2; long2 = (iGrid * iGrid + 1) / 2;

    /*
     * STEP 2: calculate
     *         1) SAME on all processes:
     *            mpi_wrstar1-2[mpi_num_procs], mpi_wrfinl1-2[mpi_num_procs], mpi_wrlen1-2[mpi_num_procs], mpi_wrnstafi1-2[mpi_num_procs]
     *            for MPI_GATHERV  use: mpi_recvcounts1-2[mpi_num_procs], mpi_recvdispls1-2[mpi_num_procs]
     *            for MPI_SCATTERV use: mpi_senddispls2l[mpi_num_procs], mpi_sendcounts2l[mpi_num_procs],
     *                                  mpi_senddispls2r[mpi_num_procs], mpi_sendcounts2r[mpi_num_procs],
     *                                  mpi_senddispls1l[mpi_num_procs], mpi_sendcounts1l[mpi_num_procs],
     *                                  mpi_senddispls1r[mpi_num_procs], mpi_sendcounts1r[mpi_num_procs]
     *         2) VARIOUS on individual process:
     *            sta1-2[iGrid], fi1-2[iGrid]
     */
    mpi_mrlen1 = fi1[iGrid - 2] - sta1[1] + 1; mpi_wrlen1[0] = mpi_mrlen1 / mpi_num_workers + 1; mpi_wrsplit1 = mpi_mrlen1 % mpi_num_workers;
    mpi_mrlen2 = fi2[iGrid - 2] - sta2[1] + 1; mpi_wrlen2[0] = mpi_mrlen2 / mpi_num_workers + 1; mpi_wrsplit2 = mpi_mrlen2 % mpi_num_workers;

    if (0 == mpi_rank) /* master process checks whether there are too many slave processes  */
    {
        if (  mpi_mrlen1 < mpi_num_workers || mpi_mrlen2 < mpi_num_workers || mpi_wrlen1[0] < (iGrid * iGrid + 1) / 2 
           || mpi_wrlen2[0] < (iGrid * iGrid + 1) / 2)
            throw CTooManyProcs(mpi_num_procs);
    }

    mpi_wrstar1[0] = 0; // mpi_wrfinl1[0] = 0; mpi_wrnstafi1[0] = 0;
    mpi_wrstar2[0] = 0; // mpi_wrfinl2[0] = 0; mpi_wrnstafi2[0] = 0;
    for (int mpi_worker = 1; mpi_worker != mpi_num_procs; mpi_worker++)
    {
        /*
         * SAME mpi_wrstar1[mpi_num_procs], mpi_wrfinl1[mpi_num_procs], mpi_wrlen1[mpi_num_procs], mpi_wrnstafi1[mpi_num_procs]
         */
        if (mpi_worker <= mpi_wrsplit1)
        {
            mpi_wrstar1[mpi_worker] = (sta1[1] + (mpi_worker - 1) * mpi_wrlen1[0]);
            mpi_wrfinl1[mpi_worker] = mpi_wrstar1[mpi_worker] + mpi_wrlen1[0] - 1;
        }
        else
        {
            mpi_wrstar1[mpi_worker] = sta1[1] + mpi_wrsplit1 * mpi_wrlen1[0] + (mpi_worker - mpi_wrsplit1 - 1) * (mpi_wrlen1[0] - 1);
            mpi_wrfinl1[mpi_worker] = mpi_wrstar1[mpi_worker] + (mpi_wrlen1[0] - 1) - 1;
        }

        mpi_wrlen1[mpi_worker] = mpi_wrfinl1[mpi_worker] - mpi_wrstar1[mpi_worker] + 1;

        for (iw = 1; iw < iGrid - 2; iw++)
        {
            if (sta1[iw] <= mpi_wrstar1[mpi_worker] && mpi_wrstar1[mpi_worker] <= sta1[iw + 1])
            {
                ix = iw;
                break;
            }
        }
        if (ix < 1) ix = 1;
        
        for (iw = 1; iw < iGrid - 2; iw++)
        {
            if (fi1[iw] <= mpi_wrfinl1[mpi_worker] && mpi_wrfinl1[mpi_worker] <= fi1[iw + 1])
            {
                iy = iw + 1;
                break;
            }
        }
        if (iy > iGrid - 2) iy = iGrid - 3;
        
        /*
         * INDIVIDUAL sta1 and fi1
         */
        if (mpi_worker == mpi_rank)
        {
            vector<delphi_integer> tmp_sta1, tmp_fi1;

            tmp_sta1.assign(iGrid, 0);
            tmp_fi1.assign(iGrid, 0);
            iz = 1;
            for (iw = ix; iw <= iy; iw++)
            {
                tmp_sta1[iz] = sta1[iw];
                tmp_fi1[iz]  =  fi1[iw];
                iz++;
            }
            tmp_sta1[1] = mpi_wrstar1[mpi_rank]; tmp_fi1[iz-1] = mpi_wrfinl1[mpi_rank];
            
            sta1 = tmp_sta1;
            fi1  = tmp_fi1;
            
            iz -= 1;
            mpi_wrnstafi1[mpi_rank] = iz;   
        }

        /*
         * SAME mpi_wrstar2[mpi_num_procs], mpi_wrfinl2[mpi_num_procs], mpi_wrlen2[mpi_num_procs], mpi_wrnstafi2[mpi_num_procs]
         */
        if (mpi_worker <= mpi_wrsplit2)
        {
            mpi_wrstar2[mpi_worker] = sta2[1] + (mpi_worker - 1) * mpi_wrlen2[0];
            mpi_wrfinl2[mpi_worker] = mpi_wrstar2[mpi_worker] + mpi_wrlen2[0] - 1;
        }
        else
        {
            mpi_wrstar2[mpi_worker] = sta2[1] + mpi_wrsplit2 * mpi_wrlen2[0] + (mpi_worker - mpi_wrsplit2 - 1) * (mpi_wrlen2[0] - 1);
            mpi_wrfinl2[mpi_worker] = mpi_wrstar2[mpi_worker] + (mpi_wrlen2[0] - 1) - 1;
        }

        mpi_wrlen2[mpi_worker] = mpi_wrfinl2[mpi_worker] - mpi_wrstar2[mpi_worker] + 1;

        for (iw = 1; iw < iGrid - 2; iw++)
        {
            if (sta2[iw] <= mpi_wrstar2[mpi_worker] && mpi_wrstar2[mpi_worker] <= sta2[iw + 1])
            {
                ix = iw;
                break;
            }
        }
        if (ix < 1) ix = 1;
        
        for (iw = 1; iw < iGrid - 2; iw++)
        {
            if (fi2[iw] <= mpi_wrfinl2[mpi_worker] && mpi_wrfinl2[mpi_worker] <= fi2[iw + 1])
            {
                iy = iw + 1;
                break;
            }
        }
        if (iy > iGrid - 2) iy = iGrid - 3;
        
        /*
         * INDIVIDUAL sta2 and fi2
         */
        if (mpi_worker == mpi_rank)
        {
            vector<delphi_integer> tmp_sta2, tmp_fi2;

            tmp_sta2.assign(iGrid, 0);
            tmp_fi2.assign(iGrid, 0);
            iz = 1;
            for (iw = ix; iw <= iy; iw++)
            {
                tmp_sta2[iz] = sta2[iw];
                tmp_fi2[iz]  =  fi2[iw];
                iz++;
            }
            tmp_sta2[1] = mpi_wrstar2[mpi_rank]; tmp_fi2[iz-1] = mpi_wrfinl2[mpi_rank];
            
            sta2 = tmp_sta2;
            fi2  = tmp_fi2;
            
            iz -= 1;
            mpi_wrnstafi2[mpi_rank] = iz;  
        }
 
        /*
         * STEP 3: send sf1 and sf2 piece to each slave processes
         *
         *         prgfSaltMap1(sf1) and prgfSaltMap2(sf2) are set in setDielecBndySaltMap and need to be updated
         *         here on each slave process.
         *
         * Notice: indices of sf1-2 start at 0 in C++, while they are sf1(WRstar1:WRfinl1),sf2(WRstar2:WRfinl2) in Fortran!
         */
        if (fZero < abs(fIonStrength)) // when fIonStrength == 0.0, prgfSaltMap11 and prgfSaltMap12 are of size 0
        {
            if (0 == mpi_rank) /* master process */
            {
                mpi_sf1 = prgfSaltMap1.data();
                MPI_Send(mpi_sf1 + mpi_wrstar1[mpi_worker] - 1, mpi_wrlen1[mpi_worker], mpi_delphi_real, mpi_worker, mpi_worker, MPI_COMM_WORLD);

                mpi_sf2 = prgfSaltMap2.data();
                MPI_Send(mpi_sf2 + mpi_wrstar2[mpi_worker] - 1, mpi_wrlen2[mpi_worker], mpi_delphi_real, mpi_worker, mpi_worker, MPI_COMM_WORLD);
            }
            else if (mpi_worker == mpi_rank)
            {
                prgfSaltMap1.assign(mpi_wrlen1[mpi_rank], 0.0);
                mpi_sf1 = prgfSaltMap1.data() - mpi_wrstar1[mpi_rank]; // mpi_sf1 is shifted mpi_wrstar1[mpi_rank] units to the negative
                                                                       // direction of prgfSaltMap1.begin() so that mpi_sf1 is
                                                                       // accommodate to sf1(WRstar1:WRfinl1) in Fortran
                MPI_Recv(mpi_sf1 + mpi_wrstar1[mpi_rank], mpi_wrlen1[mpi_rank], mpi_delphi_real, 0, mpi_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                prgfSaltMap2.assign(mpi_wrlen2[mpi_rank], 0.0);
                mpi_sf2 = prgfSaltMap2.data() - mpi_wrstar2[mpi_rank]; // mpi_sf2 is shifted mpi_wrstar2[mpi_rank] units to the negative
                                                                       // direction of prgfSaltMap2.begin() so that mpi_sf2 is
                                                                       // accommodate to sf2(WRstar2:WRfinl2) in Fortran
                MPI_Recv(mpi_sf2 + mpi_wrstar2[mpi_rank], mpi_wrlen2[mpi_rank], mpi_delphi_real, 0, mpi_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
 
        /*
         * SAME mpi_recvcounts1-2 and mpi_recvdispls1-2 for MPI_GATHERV used on master
         */
        mpi_recvdispls1[mpi_worker] = mpi_wrstar1[mpi_worker] - 1; mpi_recvcounts1[mpi_worker] = mpi_wrlen1[mpi_worker];
        mpi_recvdispls2[mpi_worker] = mpi_wrstar2[mpi_worker] - 1; mpi_recvcounts2[mpi_worker] = mpi_wrlen2[mpi_worker];

        /*
         * SAME mpi_senddispls2l, mpi_sendcounts2l, mpi_senddispls2r, mpi_sendcounts2r, mpi_senddispls1l, mpi_sendcounts1l,
         * mpi_senddispls1r, mpi_sendcounts1r for MPI_SCATTERV used on master
         */
        if (mpi_wrstar1[mpi_worker] - (iGrid * iGrid + 1) / 2 <= mpi_wrstar2[mpi_worker])
        {
            mpi_senddispls2l[mpi_worker] = mpi_wrstar1[mpi_worker] - (iGrid * iGrid + 1) / 2 - 1;
            mpi_sendcounts2l[mpi_worker] = mpi_wrstar2[mpi_worker] - (mpi_wrstar1[mpi_worker] - (iGrid * iGrid + 1) / 2);
        }
        else
        {
            mpi_sendcounts2l[mpi_worker] = 0;
            mpi_senddispls2l[mpi_worker] = 0;
        }

        if (mpi_wrfinl2[mpi_worker] <= mpi_wrfinl1[mpi_worker] + (iGrid * iGrid - 1) / 2)
        {
            mpi_senddispls2r[mpi_worker] = mpi_wrfinl2[mpi_worker];
            mpi_sendcounts2r[mpi_worker] = mpi_wrfinl1[mpi_worker] + (iGrid * iGrid - 1) / 2 - mpi_wrfinl2[mpi_worker];
        }
        else
        {
            mpi_senddispls2r[mpi_worker] = 0;
            mpi_sendcounts2r[mpi_worker] = 0;
        }

        if (mpi_wrstar2[mpi_worker] - (iGrid * iGrid - 1) / 2 <= mpi_wrstar1[mpi_worker])
        {
            mpi_senddispls1l[mpi_worker] = mpi_wrstar2[mpi_worker] - (iGrid * iGrid - 1) / 2 - 1;
            mpi_sendcounts1l[mpi_worker] = mpi_wrstar1[mpi_worker] - (mpi_wrstar2[mpi_worker] - (iGrid * iGrid - 1) / 2);
        }
        else
        {
            mpi_senddispls1l[mpi_worker] = 0;
            mpi_sendcounts1l[mpi_worker] = 0;
        }

        if (mpi_wrfinl1[mpi_worker] <= mpi_wrfinl2[mpi_worker] + (iGrid * iGrid + 1) / 2)
        {
            mpi_senddispls1r[mpi_worker] = mpi_wrfinl1[mpi_worker];
            mpi_sendcounts1r[mpi_worker] = mpi_wrfinl2[mpi_worker] + (iGrid * iGrid + 1) / 2 - mpi_wrfinl1[mpi_worker];
        }
        else
        {
            mpi_senddispls1r[mpi_worker] = 0;
            mpi_sendcounts1r[mpi_worker] = 0;
        }
    } // end of for (int mpi_worker = 1; mpi_worker != mpi_num_procs; mpi_worker++)    

    /*
     * STEP 4: construct local db, idpos, qval, iqpos and gval
     */
    if (0 == mpi_rank)
    {
        mpi_size = prgfBndyDielec.size();                                                                
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);                                                    
        for (int k = 0; k < prgfBndyDielec.size(); k++)
            MPI_Bcast(prgfBndyDielec[k].data(), 6, mpi_delphi_real, 0, MPI_COMM_WORLD);                // db    - broadcast(world, prgfBndyDielec, 0);         
        
        mpi_size = prgiBndyDielecIndex.size();
        MPI_Bcast(&mpi_size,                                           1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(prgiBndyDielecIndex.data(), prgiBndyDielecIndex.size(), MPI_INT, 0, MPI_COMM_WORLD); // idpos - broadcast(world, prgiBndyDielecIndex, 0);
        
        mpi_size = prgfCrgValA.size();
        MPI_Bcast(&mpi_size,                           1,         MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(prgfCrgValA.data(), prgfCrgValA.size(), mpi_delphi_real, 0, MPI_COMM_WORLD);         // qval  - broadcast(world, prgfCrgValA, 0);  
        
        mpi_size = prgiCrgPose.size();
        MPI_Bcast(&mpi_size,                           1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(prgiCrgPose.data(), prgiCrgPose.size(), MPI_INT, 0, MPI_COMM_WORLD);                 // iqpos - broadcast(world, prgiCrgPose, 0);
        
        mpi_size = prgfCrgValG.size();
        MPI_Bcast(&mpi_size,                           1,         MPI_INT, 0, MPI_COMM_WORLD);          
        MPI_Bcast(prgfCrgValG.data(), prgfCrgValG.size(), mpi_delphi_real, 0, MPI_COMM_WORLD);         // gval  - broadcast(world, prgfCrgValG, 0);                        
    }
    else
    {
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);                      
        vector< vector<delphi_real> > tdb(mpi_size, vector<delphi_real>(6));  
        for (int k = 0; k < mpi_size; k++)                                          // temporary db - broadcast(world, tdb,    0);
            MPI_Bcast(tdb[k].data(), 6, mpi_delphi_real, 0, MPI_COMM_WORLD);  
        
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        vector<delphi_integer> tidpos(mpi_size, 0);  
        MPI_Bcast(tidpos.data(), tidpos.size(), MPI_INT, 0, MPI_COMM_WORLD);       // temporary idpos - broadcast(world, tidpos, 0);
        
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        vector<delphi_real> tqval(mpi_size, 0.0);  
        MPI_Bcast(tqval.data(), tqval.size(), mpi_delphi_real, 0, MPI_COMM_WORLD); // temporary qval - broadcast(world, tqval,  0);
        
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        vector<delphi_integer> tiqpos(mpi_size, 0);
        MPI_Bcast(tiqpos.data(), tiqpos.size(), MPI_INT, 0, MPI_COMM_WORLD);       // temporary iqpos - broadcast(world, tiqpos, 0);
        
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        vector<delphi_real> tgval(mpi_size, 0.0);    
        MPI_Bcast(tgval.data(), tgval.size(), mpi_delphi_real, 0, MPI_COMM_WORLD); // temporary gval - broadcast(world, tgval,  0);
        
        vector< vector<delphi_real> >& db = prgfBndyDielec;                         // local db
        vector<delphi_integer>& idpos     = prgiBndyDielecIndex;                    // local idpos
        vector<delphi_real>& qval         = prgfCrgValA;                            // local qval
        vector<delphi_integer>& iqpos     = prgiCrgPose;                            // local iqpos
        vector<delphi_real>& gval         = prgfCrgValG;                            // local gval

        /*
         * if pts on dielectric bdy fall into this WR
         */
        for (int k = 0; k < iDielecBndyEven; k++) // icount1b
        {
            ix = tidpos[k];
            if (mpi_wrstar1[mpi_rank] <= ix && ix <= mpi_wrfinl1[mpi_rank])
            {
                idpos.push_back(ix);
                db.push_back(tdb[k]);
            }
        }    
        
        mpi_wricount2a = idpos.size();
        
        for (int k = iDielecBndyEven; k < iDielecBndyOdd; k++) //icount2b
        {
            ix = tidpos[k];
            if (mpi_wrstar1[mpi_rank] <= ix && ix <= mpi_wrfinl1[mpi_rank])
            {
                idpos.push_back(ix);
                db.push_back(tdb[k]);
            }
        }
        
        mpi_wricount2b = idpos.size(); 
        
        for (int k = 0; k < iCrgedGridEven; k++) // icount1a
        {
            ix = tiqpos[k];
            if (mpi_wrstar1[mpi_rank] <= ix && ix <= mpi_wrfinl1[mpi_rank])
            {
                iqpos.push_back(ix);
                qval.push_back(tqval[k]);
                gval.push_back(tgval[k]);
            }
        }
        
        mpi_wricount1a = iqpos.size();

        for (int k = iCrgedGridEven; k < iCrgedGridSum; k++) // icount1b
        {
            ix = tiqpos[k];
            if (mpi_wrstar2[mpi_rank] <= ix && ix <= mpi_wrfinl2[mpi_rank])
            {
                iqpos.push_back(ix);
                qval.push_back(tqval[k]);
                gval.push_back(tgval[k]);
            }
        }

        mpi_wricount1b = iqpos.size();

        vector<vector<delphi_real> >().swap(tdb);   // clear tdb reallocating
        vector<delphi_integer>().swap(tidpos);      // clear tidpos reallocating
        vector<delphi_real>().swap(tqval);          // clear tqval reallocating
        vector<delphi_integer>().swap(tiqpos);      // clear tiqpos reallocating
        vector<delphi_real>().swap(tgval);          // clear tgval reallocating
    }

    /*
     * STEP 5: set up variables for talking to neighbors
     */
    if (0 == mpi_rank)
    {
        /*
         * for phimap1 updating on the 1st worker
         */
        mpi_ltostart2 = mpi_wrstar1[1] - (iGrid * iGrid + 1) / 2 - 1;
        mpi_ltosize2  = mpi_wrstar2[1] - mpi_ltostart2;
        if (0 > mpi_ltosize2) mpi_ltosize2 = 0;

        /*
         * for phimap2 updating on the 1st worker
         */
        mpi_ltostart1 = mpi_wrstar2[1] - (iGrid * iGrid + 1) / 2 - 1;
        mpi_ltosize1  = mpi_wrstar1[1] - mpi_ltostart1;
        if (0 > mpi_ltosize1) mpi_ltosize1 = 0;

        /*
         * for phimap1 updating on the last worker
         */
        mpi_rtostart2 = mpi_wrfinl1[mpi_num_workers];
        mpi_rtosize2  = mpi_wrfinl1[mpi_num_workers] + (iGrid * iGrid + 1) / 2 - mpi_wrfinl2[mpi_num_workers];
        if (0 > mpi_rtosize2) mpi_rtosize2 = 0;

        /*
         * for phimap2 updating on the last worker
         */
        mpi_rtostart1 = mpi_wrfinl2[mpi_num_workers];
        mpi_rtosize1 = mpi_wrfinl2[mpi_num_workers] + (iGrid * iGrid + 1) / 2 - mpi_wrfinl1[mpi_num_workers];
        if (0 > mpi_rtosize1) mpi_rtosize1 = 0;
    
        mpi_phimap1 = phimap1.data(); mpi_wrphimap1_start = 0;
        mpi_phimap2 = phimap2.data(); mpi_wrphimap2_start = 0;
    }
    else
    {
        /*
         * for updating phimap1 on the left of this worker
         */
        mpi_lfromstart2 = mpi_wrstar1[mpi_rank] - (iGrid * iGrid + 1) / 2 - 1;
        mpi_lfromsize2  = mpi_wrstar2[mpi_rank] - mpi_lfromstart2;
        if (0 > mpi_lfromsize2) mpi_lfromsize2 = 0;

        /*
         * for updating phimap2 on the left of this worker
         */
        mpi_lfromstart1 = mpi_wrstar2[mpi_rank] - (iGrid * iGrid + 1) / 2 - 1;
        mpi_lfromsize1  = mpi_wrstar1[mpi_rank] - mpi_lfromstart1;
        if (0 > mpi_lfromsize1) mpi_lfromsize1 = 0;

        /*
         * for updating phimap1 on the right of this worker
         */
        mpi_rfromstart2 = mpi_wrfinl2[mpi_rank];
        mpi_rfromsize2  = mpi_wrfinl1[mpi_rank] + (iGrid * iGrid + 1) / 2 - mpi_wrfinl2[mpi_rank];
        if (0 > mpi_rfromsize2) mpi_rfromsize2 = 0;

        /*
         * for updating phimap2 on the right of this worker
         */
        mpi_rfromstart1 = mpi_wrfinl1[mpi_rank];
        mpi_rfromsize1  = mpi_wrfinl2[mpi_rank] + (iGrid * iGrid + 1) / 2 - mpi_wrfinl1[mpi_rank];
        if (0 > mpi_rtosize1) mpi_rtosize1 = 0;

        /*
         * allocate memory for phimap1 and phimap2 on slaves
         */
        if (rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2] || mpi_num_workers == 1) 
        {
            mpi_wrphimap1_start = 0;
            mpi_wrphimap1_end   = iHalfGridNum - 1;
            phimap1.assign(iHalfGridNum, 0.0);
            mpi_phimap1 = phimap1.data();
            
            mpi_wrphimap2_start = 0;
            mpi_wrphimap2_end   = iHalfGridNum - 1;
            phimap2.assign(iHalfGridNum, 0.0);
            mpi_phimap2 = phimap2.data();
        }
        else
        {
            mpi_wrphimap1_start = min(mpi_wrstar1[mpi_rank], mpi_lfromstart1); 
            mpi_wrphimap1_end   = max(mpi_wrfinl1[mpi_rank], mpi_rfromstart1 + mpi_rfromsize1);
            phimap1.assign(mpi_wrphimap1_end - mpi_wrphimap1_start + 1, 0.0);
            mpi_phimap1 = phimap1.data() - mpi_wrphimap1_start; // mpi_phimap1 is shifted mpi_wrphimap1_start units to 
                                                                // the negative direction of phimap1.begin() so that mpi_phimap1
                                                                // is accommodate to 
                                                                // phimap1(min(WRstar1, lfromstart1):max(WRfinl1, rfromstart1+rfromsize1)) 
                                                                // in Fortran 
            
            mpi_wrphimap2_start = min(mpi_wrstar2[mpi_rank], mpi_lfromstart2); 
            mpi_wrphimap2_end   = max(mpi_wrfinl2[mpi_rank], mpi_rfromstart2 + mpi_rfromsize2);
            phimap2.assign(mpi_wrphimap2_end - mpi_wrphimap2_start + 1, 0.0);
            mpi_phimap2 = phimap2.data() - mpi_wrphimap2_start; // mpi_phimap2 is shifted mpi_wrphimap2_start units to 
                                                                // the negative direction of phimap2.begin() so that mpi_phimap2
                                                                // is accommodate to 
                                                                // phimap2(min(WRstar2, lfromstart2):max(WRfinl2, rfromstart2+rfromsize2)) 
                                                                // in Fortran            
        }  
        
        if (1 < mpi_num_workers) /* more than one worker */
        {
            if (1 == mpi_rank)
            {
                /*
                 * communicate with the next worker on the right
                 */
                mpi_recver = mpi_rank + 1; mpi_sendtag = mpi_msgtag(mpi_rank, mpi_recver);
                mpi_sender = mpi_rank + 1; mpi_recvtag = mpi_msgtag(mpi_sender, mpi_rank);

                MPI_Sendrecv(&mpi_rfromstart2, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtostart2,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_rfromsize2 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtosize2 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_rfromstart1, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtostart1,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_rfromsize1 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtosize1 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
            }
            else if (1 < mpi_rank && mpi_rank <= mpi_num_workers / 2)
            {
                /*
                 * communicate with previous worker on the left
                 */
                mpi_recver  = mpi_rank - 1; mpi_sendtag = mpi_msgtag(mpi_rank, mpi_recver);
                mpi_sender  = mpi_rank - 1; mpi_recvtag = mpi_msgtag(mpi_sender, mpi_rank);

                MPI_Sendrecv(&mpi_lfromstart2, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_ltostart2,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_lfromsize2 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_ltosize2 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_lfromstart1, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_ltostart1,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_lfromsize1 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_ltosize1 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);

                /*
                 * communicate with the next worker on the right
                 */
                mpi_recver  = mpi_rank + 1; mpi_sendtag = mpi_msgtag(mpi_rank, mpi_recver);
                mpi_sender  = mpi_rank + 1; mpi_recvtag = mpi_msgtag(mpi_sender, mpi_rank);

                MPI_Sendrecv(&mpi_rfromstart2, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtostart2,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_rfromsize2 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtosize2 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_rfromstart1, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtostart1,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_rfromsize1 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtosize1 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
            }
            else if (mpi_num_workers / 2 < mpi_rank && mpi_rank < mpi_num_workers)
            {
                /*
                 * communicate with the next worker on the right
                 */
                mpi_recver  = mpi_rank + 1; mpi_sendtag = mpi_msgtag(mpi_rank, mpi_recver);
                mpi_sender  = mpi_rank + 1; mpi_recvtag = mpi_msgtag(mpi_sender, mpi_rank);

                MPI_Sendrecv(&mpi_rfromstart2, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtostart2,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_rfromsize2 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtosize2 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_rfromstart1, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtostart1,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_rfromsize1 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_rtosize1 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);

                /*
                 * communicate with previous worker on the left
                 */
                mpi_recver  = mpi_rank - 1; mpi_sendtag = mpi_msgtag(mpi_rank, mpi_recver);
                mpi_sender  = mpi_rank - 1; mpi_recvtag = mpi_msgtag(mpi_sender, mpi_rank);

                MPI_Sendrecv(&mpi_lfromstart2, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_ltostart2,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_lfromsize2 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_ltosize2 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_lfromstart1, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_ltostart1,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_lfromsize1 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_ltosize1 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag, 
                             MPI_COMM_WORLD, &mpi_stat);
            }
            else
            {
                /*
                 * communicate with previous worker on the left
                 */
                mpi_recver  = mpi_rank - 1; mpi_sendtag = mpi_msgtag(mpi_rank, mpi_recver);
                mpi_sender  = mpi_rank - 1; mpi_recvtag = mpi_msgtag(mpi_sender, mpi_rank);

                MPI_Sendrecv(&mpi_lfromstart2, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag, 
                             &mpi_ltostart2,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag,
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_lfromsize2 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag,
                             &mpi_ltosize2 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag,
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_lfromstart1, 1, mpi_delphi_integer, mpi_recver, mpi_sendtag,
                             &mpi_ltostart1,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag,
                             MPI_COMM_WORLD, &mpi_stat);
                MPI_Sendrecv(&mpi_lfromsize1 , 1, mpi_delphi_integer, mpi_recver, mpi_sendtag,
                             &mpi_ltosize1 ,   1, mpi_delphi_integer, mpi_sender, mpi_recvtag,
                             MPI_COMM_WORLD, &mpi_stat);
            }
        }
    }

    /*
     * STEP 6: create MPI groups for shared memory access
     *
     * We use active target synchronization MPI calls MPI_Win_start(), MPI_Win_complete(), MPI_ Win_post(),
     * and MPI_Win_wait() calls. The scope of synchronization can be restricted to only a pair of
     * communicating processes (which may be more efficient when communication occurs with a small number
     * of neighboring processes).
     *
     * Refer http://www.linux-mag.com/id/1793/ for more details and examples.
     */
    if (1 < mpi_num_workers) // more than one slave process
    {
        /*
         * Compute the rank of the process, which will access our memory,and the rank exposing the memory
         * that we are going to access:
         * If we have rank x, then we would like to access rank x-1 and to let rank x+1 access our exposed memory.
         */
        MPI_Type_size(mpi_delphi_real, &mpi_sizeofdouble);
        
        if (0 == mpi_rank)
        {
            /*
             * create window encompassing all of x array
             */
            MPI_Info_create(&mpi_info);
            MPI_Info_set(mpi_info, "no_locks", "true");

            /*
             * No memory to be shared on master
             */
            MPI_Win_create(NULL, 0, 1, mpi_info, MPI_COMM_WORLD, &mpi_towin1);
            MPI_Win_create(NULL, 0, 1, mpi_info, MPI_COMM_WORLD, &mpi_towin2);
        }
        else /* slave processes */
        {
            /*
             * create window encompassing all of x array
             */
            MPI_Info_create(&mpi_info);
            MPI_Info_set(mpi_info, "no_locks", "false");

            mpi_todispl1 = (mpi_lfromsize1-1) * mpi_sizeofdouble;
            mpi_todispl2 = (mpi_lfromsize2-1) * mpi_sizeofdouble;

            if (mpi_rank < mpi_num_workers)
            {
                /*
                 * Create the window (same size for all processes.
                 * Note that we use a displacement of 1. This means that displacements in Puts and Gets
                 * have to be calculated as byte offsets.
                 */
                mpi_towin1start = min(mpi_rfromstart1, mpi_rtostart1);
                mpi_towin1size  = (mpi_rfromsize1 + mpi_rtosize1) * mpi_sizeofdouble; // in bytes
                MPI_Win_create(mpi_phimap1 + mpi_towin1start, mpi_towin1size, 1, mpi_info, MPI_COMM_WORLD, &mpi_towin1);

                mpi_towin2start = min(mpi_rfromstart2, mpi_rtostart2);
                mpi_towin2size  = (mpi_rfromsize2 + mpi_rtosize2) * mpi_sizeofdouble; // in bytes
                MPI_Win_create(mpi_phimap2 + mpi_towin2start, mpi_towin2size, 1, mpi_info, MPI_COMM_WORLD, &mpi_towin2);

//                mpi_recver  = mpi_rank + 1; mpi_sendtag = mpi_msgtag(mpi_rank, mpi_recver);
//                MPI_Send(&mpi_towin1, 1, MPI_INT, mpi_recver, mpi_sendtag, MPI_COMM_WORLD);
//                MPI_Send(&mpi_towin2, 1, MPI_INT, mpi_recver, mpi_sendtag, MPI_COMM_WORLD);
//
//                if (1 < mpi_rank)
//                {
//                    mpi_sender  = mpi_rank - 1; mpi_recvtag = mpi_msgtag(mpi_sender, mpi_rank);
//                    MPI_Recv(&mpi_fromwin1, 1, MPI_INT, mpi_sender, mpi_recvtag, MPI_COMM_WORLD, &mpi_stat);
//                    MPI_Recv(&mpi_fromwin2, 1, MPI_INT, mpi_sender, mpi_recvtag, MPI_COMM_WORLD, &mpi_stat);
//                }        
            }
            else /* No memory to be shared on the last worker */
            {
                MPI_Win_create(NULL, 0, 1, mpi_info, MPI_COMM_WORLD, &mpi_towin1);
                MPI_Win_create(NULL, 0, 1, mpi_info, MPI_COMM_WORLD, &mpi_towin2);

//                mpi_sender  = mpi_rank - 1; mpi_recvtag = mpi_msgtag(mpi_sender, mpi_rank);
//                MPI_Recv(&mpi_fromwin1, 1, MPI_INT, mpi_sender, mpi_recvtag, MPI_COMM_WORLD, &mpi_stat);
//                MPI_Recv(&mpi_fromwin2, 1, MPI_INT, mpi_sender, mpi_recvtag, MPI_COMM_WORLD, &mpi_stat);            
            }     
        } // end of running on slave processes          
        
        /*
         * create communication groups
         */        
        mpi_startrank = (mpi_num_procs + mpi_rank - 1)%mpi_num_procs;     // start-complete are one pair   
        mpi_postrank  = (mpi_rank + 1)%mpi_num_procs;                     // post-wait are another pair      
                
        /*
         * Get the group of MPI_COMM_WORLD and create groups for the Post and Start ranks
         */
        MPI_Comm_group(MPI_COMM_WORLD, &mpi_wholegroup);
        MPI_Group_incl(mpi_wholegroup, 1, &mpi_startrank, &mpi_startgroup);
        MPI_Group_incl(mpi_wholegroup, 1, &mpi_postrank,  &mpi_postgroup );
      
    } // end of if ( 1 < mpi_num_workers)

    /*
     * STEP 7: broadcast boundary values, bndx1-4 to all slaves 
     */
    mpi_size = 0;
        
    if (0 == mpi_rank)
    {
        mpi_size = bndx1.size();
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(bndx1.data(), mpi_size, mpi_delphi_real, 0, MPI_COMM_WORLD);
        
        mpi_size = bndx2.size();
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(bndx2.data(), mpi_size, mpi_delphi_real, 0, MPI_COMM_WORLD);
        
        mpi_size = bndx3.size();
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(bndx3.data(), mpi_size, mpi_delphi_real, 0, MPI_COMM_WORLD);
        
        mpi_size = bndx4.size();
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(bndx4.data(), mpi_size, mpi_delphi_real, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        bndx1.assign(mpi_size, 0.0);
        MPI_Bcast(bndx1.data(), mpi_size, mpi_delphi_real, 0, MPI_COMM_WORLD); // bndx1 - broadcast(world, bndx1, 0); 
            
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        bndx2.assign(mpi_size, 0.0);
        MPI_Bcast(bndx2.data(), mpi_size, mpi_delphi_real, 0, MPI_COMM_WORLD); // bndx2 - broadcast(world, bndx2, 0); 
            
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        bndx3.assign(mpi_size, 0.0);
        MPI_Bcast(bndx3.data(), mpi_size, mpi_delphi_real, 0, MPI_COMM_WORLD); // bndx3 - broadcast(world, bndx3, 0); 
            
        MPI_Bcast(&mpi_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        bndx4.assign(mpi_size, 0.0);
        MPI_Bcast(bndx4.data(), mpi_size, mpi_delphi_real, 0, MPI_COMM_WORLD); // bndx4 - broadcast(world, bndx4, 0); 
    }

    /*
     * STEP 8: scatter pieces of phimap1-2 to all workers
     */
    if (0 == mpi_rank)
    {
        MPI_Scatterv(mpi_phimap1, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real, 
                     MPI_IN_PLACE, 0, mpi_delphi_real, 
                     0, MPI_COMM_WORLD);
        MPI_Scatterv(mpi_phimap2, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real, 
                     MPI_IN_PLACE, 0, mpi_delphi_real, 
                     0, MPI_COMM_WORLD);

        if (!(rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2]))
        {
            /*
             * communicate with the 1st worker
             */
            mpi_recver = 1; mpi_sendtag = mpi_msgtag(mpi_rank, mpi_recver);
            MPI_Send(mpi_phimap2 + mpi_ltostart2, mpi_ltosize2, mpi_delphi_real, mpi_recver, mpi_sendtag, MPI_COMM_WORLD);
            MPI_Send(mpi_phimap1 + mpi_ltostart1, mpi_ltosize1, mpi_delphi_real, mpi_recver, mpi_sendtag, MPI_COMM_WORLD);

            /*
             * communicate with the last worker
             */
            mpi_recver = mpi_num_workers; mpi_sendtag = mpi_msgtag(mpi_rank, mpi_recver);
            MPI_Send(mpi_phimap2 + mpi_rtostart2, mpi_rtosize2, mpi_delphi_real, mpi_recver, mpi_sendtag, MPI_COMM_WORLD);
            MPI_Send(mpi_phimap1 + mpi_rtostart1, mpi_rtosize1, mpi_delphi_real, mpi_recver, mpi_sendtag, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Scatterv(mpi_phimap1 + mpi_wrphimap1_start, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real, 
                     mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_wrlen1[mpi_rank], mpi_delphi_real, 
                     0, MPI_COMM_WORLD);
        MPI_Scatterv(mpi_phimap2 + mpi_wrphimap2_start, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real, 
                     mpi_phimap2 + mpi_wrstar2[mpi_rank] - 1, mpi_wrlen2[mpi_rank], mpi_delphi_real, 
                     0, MPI_COMM_WORLD);

        if (!(rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2]) )
        {
            if (1 == mpi_rank) // 1st worker communicates with master
            {
                mpi_sender  = 0; mpi_recvtag = mpi_msgtag(mpi_sender, mpi_rank);
                MPI_Recv(mpi_phimap2 + mpi_lfromstart2, mpi_lfromsize2, mpi_delphi_real, mpi_sender, mpi_recvtag, MPI_COMM_WORLD, &mpi_stat);
                MPI_Recv(mpi_phimap1 + mpi_lfromstart1, mpi_lfromsize1, mpi_delphi_real, mpi_sender, mpi_recvtag, MPI_COMM_WORLD, &mpi_stat);
            }
            
            if (mpi_num_workers == mpi_rank) // last worker communicates with master
            {
                mpi_sender  = 0; mpi_recvtag = mpi_msgtag(mpi_sender, mpi_rank);
                MPI_Recv(mpi_phimap2 + mpi_rfromstart2, mpi_rfromsize2, mpi_delphi_real, mpi_sender, mpi_recvtag, MPI_COMM_WORLD, &mpi_stat);
                MPI_Recv(mpi_phimap1 + mpi_rfromstart1, mpi_rfromsize1, mpi_delphi_real, mpi_sender, mpi_recvtag, MPI_COMM_WORLD, &mpi_stat);                   
            }
        }
    }
    
    /*
     * qmap1 and qmap2 are new
     * 
     * on master: qmap1.assign(iHalfGridNum, 0.0); qmap2.assign(iHalfGridNum, 0.0);
     */
    qmap1.assign(phimap1.size(), 0.0); 
    mpi_wrqmap1_start = mpi_wrphimap1_start; 
    mpi_qmap1         = qmap1.data() - mpi_wrqmap1_start;
    
    qmap2.assign(phimap2.size(), 0.0);
    mpi_wrqmap2_start = mpi_wrphimap2_start; 
    mpi_qmap2         = qmap2.data() - mpi_wrqmap2_start;  
    
    /*
     * construction for nonlinear solver only
     */
    if (2 == forWhom) 
    {
        if (0 != mpi_rank)
        {
            debmap1.assign(phimap1.size(), 0.0);
            debmap2.assign(phimap2.size(), 0.0); 
        }
 
        mpi_wrdebmap1_start = mpi_wrphimap1_start; 
        mpi_debmap1         = debmap1.data() - mpi_wrdebmap1_start;
        
        mpi_wrdebmap2_start = mpi_wrphimap2_start; 
        mpi_debmap2         = debmap2.data() - mpi_wrdebmap2_start; 
        
        
        if (0 == mpi_rank)
        {       
            MPI_Scatterv(mpi_debmap1, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real, 
                         MPI_IN_PLACE, 0, mpi_delphi_real, 
                         0, MPI_COMM_WORLD);
            MPI_Scatterv(mpi_debmap2, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real, 
                         MPI_IN_PLACE, 0, mpi_delphi_real, 
                         0, MPI_COMM_WORLD);
        }
        else
        {            
            MPI_Scatterv(mpi_debmap1 + mpi_wrdebmap1_start, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real, 
                         mpi_debmap1 + mpi_wrstar1[mpi_rank] - 1, mpi_wrlen1[mpi_rank], mpi_delphi_real, 
                         0, MPI_COMM_WORLD);
            MPI_Scatterv(mpi_debmap2 + mpi_wrdebmap2_start, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real, 
                         mpi_debmap2 + mpi_wrstar2[mpi_rank] - 1, mpi_wrlen2[mpi_rank], mpi_delphi_real, 
                         0, MPI_COMM_WORLD); 
        }
    }
    
    /*
     * test code
     */
//    {
//        string strTestFile;
//        
//        if (0 == mpi_rank) strTestFile = "rank0_solver_initOddEvenItr.dat";
//        if (1 == mpi_rank) strTestFile = "rank1_solver_initOddEvenItr.dat";
//        if (2 == mpi_rank) strTestFile = "rank2_solver_initOddEvenItr.dat";
//        if (3 == mpi_rank) strTestFile = "rank3_solver_initOddEvenItr.dat";
//        
//        ofstream ofTestStream(strTestFile.c_str());
//        ofTestStream << boolalpha;
//        ofTestStream << fixed << setprecision(7);
//        
//        ix = mpi_wrphimap1_start;
//        for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); it++)
//        {
//            ofTestStream << "mpi_phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }     
//        
//        ix = mpi_wrphimap2_start;
//        for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); it++)
//        {
//            ofTestStream << "mpi_phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }     
//        
//        ofTestStream.close();
//        
//        //cout << "mpi_rank = " << mpi_rank << ", mpi_startrank = " << mpi_startrank << ", mpi_postrank = " << mpi_postrank << endl;       
//        
//        if (0 != mpi_rank && 1 < mpi_num_workers)
//        {
//            if(1 == mpi_rank) /* the 1st slave process */
//            {
//                /*
//                 * phimap1
//                 */
//                //mpi_phimap1[7447] = 7447.1; mpi_phimap1[7448] = 7448.1; mpi_phimap1[7449] = 7449.1; 
//                                       
//                MPI_Win_post(mpi_postgroup, 0, mpi_towin1);
//                MPI_Win_wait(mpi_towin1);
//                
//                /*
//                 * phimap2
//                 */
//                mpi_phimap2[7447] = 7447.1; mpi_phimap2[7448] = 7448.1; mpi_phimap2[7449] = 7449.1;
//                
//                MPI_Win_post(mpi_postgroup, 0, mpi_towin2);                    
//                MPI_Win_wait(mpi_towin2); 
//            }
//            else if (1 < mpi_rank && mpi_rank < mpi_num_workers) /* the slave processes in between */
//            {   
//                /*
//                 * phimap1
//                 */
//                MPI_Win_post(mpi_postgroup, 0, mpi_towin1);
//                
//                MPI_Win_start(mpi_startgroup, 0, mpi_towin1);
//                MPI_Put(mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_ltosize1, mpi_delphi_real, 
//                        mpi_startrank, mpi_todispl1, mpi_ltosize1, mpi_delphi_real, mpi_towin1);
//                MPI_Get(mpi_phimap1 + mpi_lfromstart1, mpi_lfromsize1, mpi_delphi_real, 
//                        mpi_startrank, mpi_zerodispl, mpi_lfromsize1, mpi_delphi_real, mpi_towin1);
//                MPI_Win_complete(mpi_towin1);
//                
//                MPI_Win_wait(mpi_towin1);
//                
//                /*
//                 * phimap2
//                 */
//                MPI_Win_post(mpi_postgroup, 0, mpi_towin2);
//                
//                MPI_Win_start(mpi_startgroup, 0, mpi_towin2);
//                MPI_Put(mpi_phimap2 + mpi_wrstar2[mpi_rank] - 1, mpi_ltosize2, mpi_delphi_real, 
//                        mpi_startrank, mpi_todispl2, mpi_ltosize2, mpi_delphi_real, mpi_towin2);
//                MPI_Get(mpi_phimap2 + mpi_lfromstart2, mpi_lfromsize2, mpi_delphi_real, 
//                        mpi_startrank, mpi_zerodispl, mpi_lfromsize2, mpi_delphi_real, mpi_towin2);
//                MPI_Win_complete(mpi_towin2);
//                
//                MPI_Win_wait(mpi_towin2);
//            }
//            else // the last slave
//            {
//                /*
//                 * phimap1
//                 */
//                //mpi_phimap1[7447] = 7447.2; mpi_phimap1[7448] = 7448.2; mpi_phimap1[7449] = 7449.2; 
//                
//                MPI_Win_start(mpi_startgroup, 0, mpi_towin1);
//                MPI_Put(mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_ltosize1, mpi_delphi_real, 
//                        mpi_startrank, mpi_todispl1,  mpi_ltosize1, mpi_delphi_real, mpi_towin1);
//                MPI_Get(mpi_phimap1 + mpi_lfromstart1, mpi_lfromsize1, mpi_delphi_real, 
//                        mpi_startrank, mpi_zerodispl, mpi_lfromsize1, mpi_delphi_real, mpi_towin1);
//                MPI_Win_complete(mpi_towin1);
//                
//                /*
//                 * phimap2
//                 */
//                mpi_phimap2[7447] = 7447.2; mpi_phimap2[7448] = 7448.2; mpi_phimap2[7449] = 7449.2; 
//                
//                MPI_Win_start(mpi_startgroup, 0, mpi_towin2);
//                MPI_Put(mpi_phimap2 + mpi_wrstar2[mpi_rank] - 1, mpi_ltosize2, mpi_delphi_real, 
//                        mpi_startrank, mpi_todispl2, mpi_ltosize2, mpi_delphi_real, mpi_towin2);
//                MPI_Get(mpi_phimap2 + mpi_lfromstart2, mpi_lfromsize2, mpi_delphi_real, 
//                        mpi_startrank, mpi_zerodispl, mpi_lfromsize2, mpi_delphi_real, mpi_towin2);
//                MPI_Win_complete(mpi_towin2);
//            }
//        } // end of running on all slave processes
//        
//        //cout << "mpi_rank = " << mpi_rank << " finishes..." << endl; 
//        
//        if (0 == mpi_rank) /* master process */
//        {            
//            MPI_Gatherv( MPI_IN_PLACE, 0, mpi_delphi_real, 
//                         phimap1.data(), mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real,
//                         0, MPI_COMM_WORLD);
//            MPI_Gatherv( MPI_IN_PLACE, 0, mpi_delphi_real, 
//                         phimap2.data(), mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real,
//                         0, MPI_COMM_WORLD);
//        }
//        else /* slave processes */
//        {
//            MPI_Gatherv( mpi_phimap1 + mpi_wrstar1[mpi_rank] - 1, mpi_wrlen1[mpi_rank], mpi_delphi_real, 
//                         mpi_phimap1 + mpi_wrphimap1_start, mpi_recvcounts1, mpi_recvdispls1, mpi_delphi_real,
//                         0, MPI_COMM_WORLD);
//            MPI_Gatherv( mpi_phimap2 + mpi_wrstar2[mpi_rank] - 1, mpi_wrlen2[mpi_rank], mpi_delphi_real, 
//                         mpi_phimap2 + mpi_wrphimap2_start, mpi_recvcounts2, mpi_recvdispls2, mpi_delphi_real,
//                         0, MPI_COMM_WORLD);
//        }
//        
//        string strTestFile2;
//
//        if (0 == mpi_rank) strTestFile2 = "rank0_solver_initOddEvenItr2.dat";
//        if (1 == mpi_rank) strTestFile2 = "rank1_solver_initOddEvenItr2.dat";
//        if (2 == mpi_rank) strTestFile2 = "rank2_solver_initOddEvenItr2.dat";
//        if (3 == mpi_rank) strTestFile2 = "rank3_solver_initOddEvenItr2.dat";
//        
//        ofstream ofTestStream2(strTestFile2.c_str());
//        ofTestStream2 << boolalpha;
//        ofTestStream2 << fixed << setprecision(7);   
// 
//        ix = mpi_wrphimap1_start;
//        for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); it++)
//        {
//            ofTestStream2 << "mpi_phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        } 
//        
//        ix = mpi_wrphimap2_start;
//        for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); it++)
//        {
//            ofTestStream2 << "mpi_phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//        
//        ofTestStream2.close();
//               
//
//        ix = 0;
//        for (vector<delphi_integer>::iterator it = ibndx.begin(); it != ibndx.end(); ++it)
//        {
//            ofTestStream << "ibndx[" << setw(6) << right << ix << "] = " << setw(8) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_integer>::iterator it = ibndy.begin(); it != ibndy.end(); ++it)
//        {
//            ofTestStream << "ibndy[" << setw(6) << right << ix << "] = " << setw(8) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_integer>::iterator it = ibndz.begin(); it != ibndz.end(); ++it)
//        {
//            ofTestStream << "ibndz[" << setw(6) << right << ix << "] = " << setw(8) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = bndx1.begin(); it != bndx1.end(); ++it)
//        {
//            ofTestStream << "bndx1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = bndx2.begin(); it != bndx2.end(); ++it)
//        {
//            ofTestStream << "bndx2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = bndx3.begin(); it != bndx3.end(); ++it)
//        {
//            ofTestStream << "bndx3[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = 0;
//        for (vector<delphi_real>::iterator it = bndx4.begin(); it != bndx4.end(); ++it)
//        {
//            ofTestStream << "bndx4[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = mpi_wrphimap1_start;
//        for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); it++)
//        {
//            ofTestStream << "mpi_phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//
//        ix = mpi_wrphimap2_start;
//        for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); it++)
//        {
//            ofTestStream << "mpi_phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }
//        
//        if (0 == mpi_rank)
//            ix = 0;
//        else
//            ix = mpi_wrstar1[mpi_rank] - 1;
//        for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); it++)
//        {
//            ofTestStream << "mpi_sf1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }        
//        
//        if (0 == mpi_rank)
//            ix = 0;
//        else
//            ix = mpi_wrstar2[mpi_rank] - 1;
//        for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); it++)
//        {
//            ofTestStream << "mpi_sf2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//            ix++;
//        }         
//
//    } // end of test code      
    
    //mpi_phimap1 += 1; mpi_phimap2 += 1;
    mpi_sf1 += 1; mpi_sf2 += 1;
    
    /*
     * Synchronization before getting into the loop
     */
    MPI_Barrier (MPI_COMM_WORLD);
    
    /*
     * Now mpi_phimap1[mpi_wrphimap1_start] is at phimap1.begin() and mpi_phimap2[mpi_wrphimap2_start] is at phimap2.begin()
     *     mpi_sf1[mpi_wrstar1[mpi_rank]] is at prgfSaltMap1.begin() and mpi_sf2[mpi_wrstar2[mpi_rank]] is at prgfSaltMap2.begin()
     */
}

#endif
