/*
 * energy_react.cpp
 *
 * Author: Lin Wang, lwang3@g.clemson.edu
 *
 * This file is function of reaction field energy (solvation) calculation.
 * It's called by energy_run function if energy function solvation(s or sol) flag set to be TRUE.
 *
 * ARGO: INCLUDES CHANGES IN OUTPUT FORMATTING ONLY
 * ARGO: NO TECHNICAL CHANGE
 */

#include "energy.h"

void CDelphiEnergy::energy_react(delphi_real& fEnergy_Solvation, delphi_real& fEnergy_AnalySurf, int& iisitpot)
{
    bool bRadiusWarn = false;
    int ix, iy, iz;
    int i, j, ii, jj, qq;
    delphi_real fRadius, fCost, fVal;
    delphi_real fEnergy_Temp = 0, fEnergy_TotalCharge = 0;
    delphi_real dist, fEnergy_Temp1, fFact, ptemp; // fEnergy_Temp2, fEnergy_Temp3,
    delphi_real temp, temp1, temp2, temp3, spt1, fConstSixth;
    delphi_real dx, dy, dz;

    vector<delphi_real> spdiv(iTotalBdyGridNum), spot(iTotalBdyGridNum), schrg_omp(iTotalBdyGridNum);
    vector<SGridValue<delphi_real> > cgrid(iTotalBdyGridNum);

    SGrid<delphi_integer> ixyz;

    fFact                = 0.9549296586 / (2.0 * fScale * fEPKT);
    fConstSixth          = 1.0 / 6.0;
    fEnergy_Temp1        = 0.0;

    /*
     * Define these variables in this scope for OpenMP
     */
    vector<delphi_real> prgfgSurfCrgA_nX(iTotalBdyGridNum);
    vector<delphi_real> prgfgSurfCrgA_nY(iTotalBdyGridNum);
    vector<delphi_real> prgfgSurfCrgA_nZ(iTotalBdyGridNum);
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> prggvAtomicCrg_nValue(iCrgGridNum);

    for (i = 0; i < iTotalBdyGridNum; i++)
    {
        ix = prgigBndyGrid[i].nX;
        iy = prgigBndyGrid[i].nY;
        iz = prgigBndyGrid[i].nZ;

        temp1 = prgfPhimap[(iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + ix]   + prgfPhimap[(iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 2)];
        temp2 = prgfPhimap[(iz - 1) * iGrid * iGrid + (iy) * iGrid + (ix - 1)] + prgfPhimap[(iz - 1) * iGrid * iGrid + (iy - 2) * iGrid + (ix - 1)];
        temp3 = prgfPhimap[(iz) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1)] + prgfPhimap[(iz - 2) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1)];
        temp  = prgfPhimap[(iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1)] - (temp1 + temp2 + temp3) * fConstSixth;

        spdiv[i]          = temp;
        cgrid[i].nGrid.nX = ix;
        cgrid[i].nGrid.nY = iy;
        cgrid[i].nGrid.nZ = iz;
        cgrid[i].nValue   = temp;
    }

    if (iCrgBdyGrid != 0)
    {
        /* The folloing Loop for variable cgrid moved to above and vector push_back has discarded for higher performance.
         SGridValue<delphi_real> fCGridValue_Temp;
         for(i=0;i<iTotalBdyGridNum;i++)
         {
             fCGridValue_Temp.nGrid.nX = prgigBndyGrid[i].nX;
             fCGridValue_Temp.nGrid.nY = prgigBndyGrid[i].nY;
             fCGridValue_Temp.nGrid.nZ = prgigBndyGrid[i].nZ;
             fCGridValue_Temp.nValue = spdiv[i];
             cgrid.push_back(fCGridValue_Temp);
         }
         */

        for (i = 0; i < iCrgBdyGrid; i++)
        {
            ix   = prgdgvCrgBndyGrid[i].fgCoord.nX;
            iy   = prgdgvCrgBndyGrid[i].fgCoord.nY;
            iz   = prgdgvCrgBndyGrid[i].fgCoord.nZ;
            fVal = prgdgvCrgBndyGrid[i].fVal1;

            for (j = 0; j < iTotalBdyGridNum; j++)
                if (cgrid[j].nGrid.nX == ix && cgrid[j].nGrid.nY == iy && cgrid[j].nGrid.nZ == iz)
                    cgrid[j].nValue = cgrid[j].nValue - fVal;
        }

        for (i = 0; i < iTotalBdyGridNum; i++)
            spdiv[i] = cgrid[i].nValue;

        vector<SGridValue<delphi_real> >().swap(cgrid);
    }

    fEnergy_Temp1 = 0.0;
    schrg.reserve(iTotalBdyGridNum);

    for (i = 0; i < iTotalBdyGridNum; i++)
    {
        temp = spdiv[i] * fFact;
        schrg[i] = temp;
        schrg_omp[i] = temp;
        fEnergy_Temp1 = fEnergy_Temp1 + temp;
    } // omp

    if (bReactFieldlnFRC || bSolvEng || bNonlinearEng || bSurfEngOut || bSurfCrgOut)
    {
        if (bSolvEng || bNonlinearEng)
        {
            fEnergy_Temp        = 0.0;
            fEnergy_TotalCharge = 0.0;
            fCost               = 0.0;


            /*
             * Next TWO LOOPS for OpenMP
             */
            for (i = 0; i < iTotalBdyGridNum; i++)
            {
                prgfgSurfCrgA_nX[i] = prgfgSurfCrgA[i].nX;
                prgfgSurfCrgA_nY[i] = prgfgSurfCrgA[i].nY;
                prgfgSurfCrgA_nZ[i] = prgfgSurfCrgA[i].nZ;
            }

            for (j = 0; j < iCrgGridNum; j++)
            {
                prgfgCrgPoseA_nX[j] = prgfgCrgPoseA[j].nX;
                prgfgCrgPoseA_nY[j] = prgfgCrgPoseA[j].nY;
                prgfgCrgPoseA_nZ[j] = prgfgCrgPoseA[j].nZ;
                prggvAtomicCrg_nValue[j] = prggvAtomicCrg[j].nValue;
            }

            // ++++++++ Starting parallel computing and evoke OpenMP +++++++++++ //
            // Refer OpenMP website for more details about tutorials and manuals //
            // http://openmp.org/wp/openmp-specifications                        //
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

#ifdef PARALLEL_OMP

#pragma omp parallel shared(prgfgSurfCrgA_nX,prgfgSurfCrgA_nY,prgfgSurfCrgA_nZ,prgfgCrgPoseA_nX,prgfgCrgPoseA_nY,prgfgCrgPoseA_nZ,prggvAtomicCrg_nValue,spot,schrg_omp) private(j,dx,dy,dz,dist,ptemp)
            {

#pragma omp for reduction( + : fEnergy_TotalCharge, fEnergy_Temp)

#endif
            for (i = 0; i < iTotalBdyGridNum; i++)
            {
                ptemp = 0.0;

                for (j = 0; j < iCrgGridNum; j++)
                {
                    dx = prgfgSurfCrgA_nX[i] - prgfgCrgPoseA_nX[j];
                    dy = prgfgSurfCrgA_nY[i] - prgfgCrgPoseA_nY[j];
                    dz = prgfgSurfCrgA_nZ[i] - prgfgCrgPoseA_nZ[j];

                    dist = sqrt(dx * dx + dy * dy + dz * dz);
                    ptemp = ptemp + prggvAtomicCrg_nValue[j] / dist;
                } // vectorization

                spot[i]              = ptemp;
                fEnergy_Temp        += spot[i] * schrg_omp[i];
                fEnergy_TotalCharge += schrg_omp[i];
            }
#ifdef PARALLEL_OMP
        }
#endif

            fEnergy_Solvation = fEnergy_Temp * fEPKT / 2.0;

            #ifdef VERBOSE
            if(iTotalBdyGridNum==0 || iGrid<5)
            {
                cout << " WARNING !!! Midpoints are out side the cube and delphi cannot determine the molecular surface." << endl;
                cout << " WARNING !!! Please enlarge the gsize or decrease the perfil value." << endl;
            }
            #endif

            cout << enerString << left << setw(MAXWIDTH) << "Corrected reaction field energy" << " : "
                    << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_Solvation << " kT" << endl;

            ergr = fEnergy_Solvation;

            if (bEngOut)
            {
                ofstream ofEnergyFile;
                ofEnergyFile.open(strEnergyFile, std::fstream::app);
                ofEnergyFile << "Corrected reaction field energy : " << fEnergy_Solvation << " kT\n";
                ofEnergyFile.close();
            }
        }

        if (bSurfCrgOut)
        {
            cout << " Writing surface charge file " << strScrgFile << endl;
            ofstream ofScrgFile;
            ofScrgFile.open(strScrgFile);

            if (iSurfCrgFormatOut == 1 || iSurfCrgFormatOut == 2)
            {
                ofScrgFile << "DELPHI FORMAT PDB" << endl;
                ofScrgFile << "FORMAT NUMBER = " << iSurfCrgFormatOut << endl;
                ofScrgFile << "       bgp#  atom SC   res#      pos                               scrg           surf energy" << endl;
            }

            /*
             eBuffz function has removed from C++ DelPhi
             if(bBuffz){
             lim_min = 2+ieBuffz.nMin;
             lim_max = iGrid-1-ieBuffz.nMax;
             }
             */

            for (i = 0; i < iTotalBdyGridNum; i++) {
                ixyz = prgigBndyGrid[i];

                /*
                 eBuffz function has removed from C++ DelPhi
                 if(bBuffz){
                 ido = 1;
                 if(optORLT<int>(ixyz,lim_min) || optORGT<int>(ixyz,lim_max)){
                 ido = 0;
                 }
                 if(ido==0) continue;
                 }
                 */
                fEnergy_Temp1 = schrg[i];
                spt1 = spot[i] * fEnergy_Temp1 * fEPKT / 2.0;

                if (iSurfCrgFormatOut == 0)
                    ofScrgFile << setw(NUMWIDTH) << fixed << right << i + 1 << "  " << prgfgSurfCrgA[i].nX << "  " << prgfgSurfCrgA[i].nY << "  "
                            << prgfgSurfCrgA[i].nZ << "  " << fEnergy_Temp1 << endl;

                if (iSurfCrgFormatOut == 1 || iSurfCrgFormatOut == 2)
                {
                    jj = atsurf[i];    //from surface construction class
                    qq = atoi( prgapAtomPdb[jj - 1].getAtInf().substr(11,4).c_str());

                    ofScrgFile << "ATOM  " << setw(5) << fixed << right << i + 1 << " " << setw(5) << fixed << right << jj << " SC "
                            << setw(6) << fixed << right << qq << "      " << setw(5) << fixed << right << prgfgSurfCrgA[i].nX
                            << "  " << setw(5) << fixed << right << prgfgSurfCrgA[i].nY << "  " << setw(5) << fixed
                            << right << prgfgSurfCrgA[i].nZ << "  " << setw(5) << fixed << right << scientific << fEnergy_Temp1
                            << "  " << setw(5) << fixed << right << scientific << spt1 << endl;
                }
            }
        }

        if (bSurfEngOut)
        {
            cout << " Writing surface energy file = surfen.dat" << endl;
            ofstream surfen;
            surfen.open("surfen.dat");
            //fEnergy_Temp2 = 0.0;

            for (i = 0; i < iTotalBdyGridNum; i++)
            {
                fEnergy_Temp1 = 0.000;

                surfen << setw(NUMWIDTH) << fixed << right << i + 1 << "  " << prgfgSurfCrgA[i].nX << "  " << prgfgSurfCrgA[i].nY
                        << "  " << prgfgSurfCrgA[i].nZ << "  " << fEnergy_Temp1 << endl;
            }

            surfen.close();
        }
    }

    /*
     * Clean memory occupied by these vectors
     */
    vector<delphi_real>().swap(prgfgSurfCrgA_nX);
    vector<delphi_real>().swap(prgfgSurfCrgA_nY);
    vector<delphi_real>().swap(prgfgSurfCrgA_nZ);
    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);

    vector<delphi_real>().swap(spdiv);
    vector<delphi_real>().swap(spot);
    vector<delphi_real>().swap(schrg_omp);
}


#ifdef PARALLEL_MPI

/*
 * MPI version of energy_react function
 */
void CDelphiEnergy::mpi_energy_react(delphi_real& fEnergy_Solvation, delphi_real& fEnergy_AnalySurf, int& iisitpot)
{
    bool bRadiusWarn = false;
    int ix, iy, iz;
    int i, j, ii, jj, qq;
    delphi_real fRadius, fCost, fVal;
    delphi_real fEnergy_Temp = 0, fEnergy_TotalCharge = 0;
    delphi_real dist, fEnergy_Temp1, fFact, ptemp; // fEnergy_Temp2, fEnergy_Temp3,
    delphi_real temp, temp1, temp2, temp3, spt1, fConstSixth;
    delphi_real dx, dy, dz;

    vector<delphi_real> spdiv(iTotalBdyGridNum), spot(iTotalBdyGridNum), schrg_omp(iTotalBdyGridNum);
    vector<SGridValue<delphi_real> > cgrid(iTotalBdyGridNum);

    SGrid<delphi_integer> ixyz;

    fFact                = 0.9549296586 / (2.0 * fScale * fEPKT);
    fConstSixth          = 1.0 / 6.0;
    fEnergy_Temp1        = 0.0;

    /*
     * Define these variables in this scope for OpenMP
     */
    vector<delphi_real> prgfgSurfCrgA_nX(iTotalBdyGridNum);
    vector<delphi_real> prgfgSurfCrgA_nY(iTotalBdyGridNum);
    vector<delphi_real> prgfgSurfCrgA_nZ(iTotalBdyGridNum);
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> prggvAtomicCrg_nValue(iCrgGridNum);

    if (0 == mpi_rank) /* master process */
    {
        for (i = 0; i < iTotalBdyGridNum; i++)
        {
            ix = prgigBndyGrid[i].nX;
            iy = prgigBndyGrid[i].nY;
            iz = prgigBndyGrid[i].nZ;

            temp1 = prgfPhimap[(iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + ix]   + prgfPhimap[(iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 2)];
            temp2 = prgfPhimap[(iz - 1) * iGrid * iGrid + (iy) * iGrid + (ix - 1)] + prgfPhimap[(iz - 1) * iGrid * iGrid + (iy - 2) * iGrid + (ix - 1)];
            temp3 = prgfPhimap[(iz) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1)] + prgfPhimap[(iz - 2) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1)];
            temp  = prgfPhimap[(iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1)] - (temp1 + temp2 + temp3) * fConstSixth;

            spdiv[i]          = temp;
            cgrid[i].nGrid.nX = ix;
            cgrid[i].nGrid.nY = iy;
            cgrid[i].nGrid.nZ = iz;
            cgrid[i].nValue   = temp;
        }

        if (iCrgBdyGrid != 0)
        {
            for (i = 0; i < iCrgBdyGrid; i++)
            {
                ix   = prgdgvCrgBndyGrid[i].fgCoord.nX;
                iy   = prgdgvCrgBndyGrid[i].fgCoord.nY;
                iz   = prgdgvCrgBndyGrid[i].fgCoord.nZ;
                fVal = prgdgvCrgBndyGrid[i].fVal1;

                for (j = 0; j < iTotalBdyGridNum; j++)
                    if (cgrid[j].nGrid.nX == ix && cgrid[j].nGrid.nY == iy && cgrid[j].nGrid.nZ == iz)
                        cgrid[j].nValue = cgrid[j].nValue - fVal;
            }

            for (i = 0; i < iTotalBdyGridNum; i++)
                spdiv[i] = cgrid[i].nValue;

            vector<SGridValue<delphi_real> >().swap(cgrid);
        }

        fEnergy_Temp1 = 0.0;
        schrg.reserve(iTotalBdyGridNum);

        for (i = 0; i < iTotalBdyGridNum; i++)
        {
            temp          = spdiv[i] * fFact;
            schrg[i]      = temp;
            schrg_omp[i]  = temp;
            fEnergy_Temp1 = fEnergy_Temp1 + temp;
        }

        if (bReactFieldlnFRC || bSolvEng || bNonlinearEng || bSurfEngOut || bSurfCrgOut)
        {
            if (bSolvEng || bNonlinearEng)
            {
                fEnergy_Temp        = 0.0;
                fEnergy_TotalCharge = 0.0;
                fCost               = 0.0;


                for (i = 0; i < iTotalBdyGridNum; i++)
                {
                    prgfgSurfCrgA_nX[i] = prgfgSurfCrgA[i].nX;
                    prgfgSurfCrgA_nY[i] = prgfgSurfCrgA[i].nY;
                    prgfgSurfCrgA_nZ[i] = prgfgSurfCrgA[i].nZ;
                }

                for (j = 0; j < iCrgGridNum; j++)
                {
                    prgfgCrgPoseA_nX[j] = prgfgCrgPoseA[j].nX;
                    prgfgCrgPoseA_nY[j] = prgfgCrgPoseA[j].nY;
                    prgfgCrgPoseA_nZ[j] = prgfgCrgPoseA[j].nZ;
                    prggvAtomicCrg_nValue[j] = prggvAtomicCrg[j].nValue;
                }

                /*
                 * ---------- parallel computing the following loop ----------
                 */
                {
                    delphi_real * mpi_prgfgSurfCrgA_nX = prgfgSurfCrgA_nX.data();
                    delphi_real * mpi_prgfgSurfCrgA_nY = prgfgSurfCrgA_nY.data();
                    delphi_real * mpi_prgfgSurfCrgA_nZ = prgfgSurfCrgA_nZ.data();

                    MPI_Bcast(mpi_prgfgSurfCrgA_nX, iTotalBdyGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
                    MPI_Bcast(mpi_prgfgSurfCrgA_nY, iTotalBdyGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
                    MPI_Bcast(mpi_prgfgSurfCrgA_nZ, iTotalBdyGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

                    delphi_real * mpi_prgfgCrgPoseA_nX = prgfgCrgPoseA_nX.data();
                    delphi_real * mpi_prgfgCrgPoseA_nY = prgfgCrgPoseA_nY.data();
                    delphi_real * mpi_prgfgCrgPoseA_nZ = prgfgCrgPoseA_nZ.data();

                    MPI_Bcast(mpi_prgfgCrgPoseA_nX, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
                    MPI_Bcast(mpi_prgfgCrgPoseA_nY, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
                    MPI_Bcast(mpi_prgfgCrgPoseA_nZ, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

                    delphi_real * mpi_schrg_omp = schrg_omp.data();

                    MPI_Bcast(mpi_schrg_omp, iTotalBdyGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

                    delphi_real * mpi_prggvAtomicCrg_nValue = prggvAtomicCrg_nValue.data();

                    MPI_Bcast(mpi_prggvAtomicCrg_nValue, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

                    /*
                     * the following loop is given to slaves for parallel computing
                     *
                    for (i = 0; i < iTotalBdyGridNum; i++)
                    {
                        ptemp = 0.0;

                        for (j = 0; j < iCrgGridNum; j++)
                        {
                            dx = prgfgSurfCrgA_nX[i] - prgfgCrgPoseA_nX[j];
                            dy = prgfgSurfCrgA_nY[i] - prgfgCrgPoseA_nY[j];
                            dz = prgfgSurfCrgA_nZ[i] - prgfgCrgPoseA_nZ[j];

                            dist = sqrt(dx * dx + dy * dy + dz * dz);
                            ptemp = ptemp + prggvAtomicCrg_nValue[j] / dist;
                        } // vectorization

                        spot[i]              = ptemp;
                        fEnergy_Temp        += spot[i] * schrg_omp[i];
                        fEnergy_TotalCharge += schrg_omp[i];
                    }*/

                    MPI_Reduce(MPI_IN_PLACE, &fEnergy_Temp,        1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);
                    MPI_Reduce(MPI_IN_PLACE, &fEnergy_TotalCharge, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);
                } /* end of parallel computing */

                fEnergy_Solvation = fEnergy_Temp * fEPKT / 2.0;

                #ifdef VERBOSE
                if(iTotalBdyGridNum==0 || iGrid<5)
                {
                    cout << " WARNING !!! Midpoints are out side the cube and delphi cannot determine the molecular surface." << endl;
                    cout << " WARNING !!! Please enlarge the gsize or decrease the perfil value." << endl;
                }
                #endif

                cout << enerString << left << setw(MAXWIDTH) << "Corrected reaction field energy" << " : "
                     << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_Solvation << " kT" << endl;

                ergr = fEnergy_Solvation;

                if (bEngOut)
                {
                    ofstream ofEnergyFile;
                    ofEnergyFile.open(strEnergyFile, std::fstream::app);
                    ofEnergyFile << "Corrected reaction field energy : " << fEnergy_Solvation << " kT\n";
                    ofEnergyFile.close();
                }
            } /* end of if (bSolvEng || bNonlinearEng) */

            if (bSurfCrgOut)
            {
                cout << " Writing surface charge file " << strScrgFile << endl;
                ofstream ofScrgFile;
                ofScrgFile.open(strScrgFile);

                if (iSurfCrgFormatOut == 1 || iSurfCrgFormatOut == 2)
                {
                    ofScrgFile << "DELPHI FORMAT PDB" << endl;
                    ofScrgFile << "FORMAT NUMBER = " << iSurfCrgFormatOut << endl;
                    ofScrgFile << "       bgp#  atom SC   res#      pos                               scrg           surf energy" << endl;
                }

                for (i = 0; i < iTotalBdyGridNum; i++)
                {
                    ixyz = prgigBndyGrid[i];

                    fEnergy_Temp1 = schrg[i];
                    spt1 = spot[i] * fEnergy_Temp1 * fEPKT / 2.0;

                    if (iSurfCrgFormatOut == 0)
                        ofScrgFile << setw(NUMWIDTH) << fixed << right << i + 1 << "  " << prgfgSurfCrgA[i].nX << "  " << prgfgSurfCrgA[i].nY << "  "
                                << prgfgSurfCrgA[i].nZ << "  " << fEnergy_Temp1 << endl;

                    if (iSurfCrgFormatOut == 1 || iSurfCrgFormatOut == 2)
                    {
                        jj = atsurf[i];    //from surface construction class
                        qq = atoi( prgapAtomPdb[jj - 1].getAtInf().substr(11,4).c_str());

                        ofScrgFile << "ATOM  " << setw(5) << fixed << right << i + 1 << " " << setw(5) << fixed << right << jj << " SC "
                                << setw(6) << fixed << right << qq << "      " << setw(5) << fixed << right << prgfgSurfCrgA[i].nX
                                << "  " << setw(5) << fixed << right << prgfgSurfCrgA[i].nY << "  " << setw(5) << fixed
                                << right << prgfgSurfCrgA[i].nZ << "  " << setw(5) << fixed << right << scientific << fEnergy_Temp1
                                << "  " << setw(5) << fixed << right << scientific << spt1 << endl;
                    }
                }
            } /* end of if (bSurfCrgOut) */

            if (bSurfEngOut)
            {
                cout << " Writing surface energy file = surfen.dat" << endl;
                ofstream surfen;
                surfen.open("surfen.dat");
                //fEnergy_Temp2 = 0.0;

                for (i = 0; i < iTotalBdyGridNum; i++)
                {
                    fEnergy_Temp1 = 0.000;

                    surfen << setw(NUMWIDTH) << fixed << right << i + 1 << "  " << prgfgSurfCrgA[i].nX << "  " << prgfgSurfCrgA[i].nY
                            << "  " << prgfgSurfCrgA[i].nZ << "  " << fEnergy_Temp1 << endl;
                }

                surfen.close();
            } /* end of if (bSurfEngOut) */

        } /* end of if (bReactFieldlnFRC || bSolvEng || bNonlinearEng || bSurfEngOut || bSurfCrgOut) */
    }
    else /* slave processes */
    {
        if (bReactFieldlnFRC || bSolvEng || bNonlinearEng || bSurfEngOut || bSurfCrgOut)
        {
            if (bSolvEng || bNonlinearEng)
            {
                fEnergy_Temp        = 0.0;
                fEnergy_TotalCharge = 0.0;
                
                /*
                 * ---------- parallel computing the following loop ----------
                 */
                delphi_real * mpi_prgfgSurfCrgA_nX = prgfgSurfCrgA_nX.data();
                delphi_real * mpi_prgfgSurfCrgA_nY = prgfgSurfCrgA_nY.data();
                delphi_real * mpi_prgfgSurfCrgA_nZ = prgfgSurfCrgA_nZ.data();

                MPI_Bcast(mpi_prgfgSurfCrgA_nX, iTotalBdyGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
                MPI_Bcast(mpi_prgfgSurfCrgA_nY, iTotalBdyGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
                MPI_Bcast(mpi_prgfgSurfCrgA_nZ, iTotalBdyGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

                delphi_real * mpi_prgfgCrgPoseA_nX = prgfgCrgPoseA_nX.data();
                delphi_real * mpi_prgfgCrgPoseA_nY = prgfgCrgPoseA_nY.data();
                delphi_real * mpi_prgfgCrgPoseA_nZ = prgfgCrgPoseA_nZ.data();

                MPI_Bcast(mpi_prgfgCrgPoseA_nX, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
                MPI_Bcast(mpi_prgfgCrgPoseA_nY, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);
                MPI_Bcast(mpi_prgfgCrgPoseA_nZ, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

                delphi_real * mpi_schrg_omp = schrg_omp.data();
                delphi_real * mpi_spot = spot.data();

                MPI_Bcast(mpi_schrg_omp, iTotalBdyGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

                delphi_real * mpi_prggvAtomicCrg_nValue = prggvAtomicCrg_nValue.data();

                MPI_Bcast(mpi_prggvAtomicCrg_nValue, iCrgGridNum, mpi_delphi_real, 0, MPI_COMM_WORLD);

                /*
                 * the following loop is performed on all slaves for parallel computing
                 */
                if ( iTotalBdyGridNum <= mpi_num_workers )
                {
                    if (mpi_rank <= iTotalBdyGridNum)
                    {
                        i = mpi_rank - 1;

                        ptemp = 0.0;

                        for (j = 0; j < iCrgGridNum; j++)
                        {
                            dx = prgfgSurfCrgA_nX[i] - prgfgCrgPoseA_nX[j];
                            dy = prgfgSurfCrgA_nY[i] - prgfgCrgPoseA_nY[j];
                            dz = prgfgSurfCrgA_nZ[i] - prgfgCrgPoseA_nZ[j];

                            dist = sqrt(dx * dx + dy * dy + dz * dz);
                            ptemp = ptemp + prggvAtomicCrg_nValue[j] / dist;
                        } // vectorization

                        spot[i]            = ptemp;
                        fEnergy_Temp        += spot[i] * schrg_omp[i];
                        fEnergy_TotalCharge += schrg_omp[i];
                    }
                }
                /*
                 * a simple parallelization but does NOT produce high speedups
                 */
//                else
//                {
//                    int increment = mpi_num_workers;
//
//                    for (i = mpi_rank - 1; i < iTotalBdyGridNum; i += increment)
//                    {
//                        ptemp = 0.0;
//
//                        for (j = 0; j < iCrgGridNum; j++)
//                        {
//                            dx = prgfgSurfCrgA_nX[i] - prgfgCrgPoseA_nX[j];
//                            dy = prgfgSurfCrgA_nY[i] - prgfgCrgPoseA_nY[j];
//                            dz = prgfgSurfCrgA_nZ[i] - prgfgCrgPoseA_nZ[j];
//
//                            dist = sqrt(dx * dx + dy * dy + dz * dz);
//                            ptemp = ptemp + prggvAtomicCrg_nValue[j] / dist;
//                        } // vectorization
//
//                        spot[i]              = ptemp;
//                        fEnergy_Temp        += spot[i] * schrg_omp[i];
//                        fEnergy_TotalCharge += schrg_omp[i];
//                    }
//                }
                /*
                 * another parallelization to improve performance by using cache more effectively
                 */
                else
                {
                    int mpi_mystart, mpi_myend, mpi_split = iTotalBdyGridNum % mpi_num_workers;

                    if (mpi_rank <= mpi_split)
                    {
                        mpi_mystart = (mpi_rank - 1) * (iTotalBdyGridNum / mpi_num_workers + 1);
                        mpi_myend   = mpi_mystart + (iTotalBdyGridNum / mpi_num_workers + 1) - 1;
                    }
                    else
                    {
                        mpi_mystart = mpi_split * (iTotalBdyGridNum / mpi_num_workers + 1) + (mpi_rank - mpi_split - 1) * (iTotalBdyGridNum / mpi_num_workers);
                        mpi_myend   = mpi_mystart + (iTotalBdyGridNum / mpi_num_workers) - 1;
                    }

                    if (mpi_myend >= iTotalBdyGridNum - 1) mpi_myend = iTotalBdyGridNum - 1;

                    for (i = mpi_mystart; i <= mpi_myend; i++)
                    {
                        ptemp = 0.0;

                        for (j = 0; j < iCrgGridNum; j++)
                        {
                            dx = prgfgSurfCrgA_nX[i] - prgfgCrgPoseA_nX[j];
                            dy = prgfgSurfCrgA_nY[i] - prgfgCrgPoseA_nY[j];
                            dz = prgfgSurfCrgA_nZ[i] - prgfgCrgPoseA_nZ[j];

                            dist = sqrt(dx * dx + dy * dy + dz * dz);
                            ptemp = ptemp + prggvAtomicCrg_nValue[j] / dist;
                        } // vectorization

                        spot[i]              = ptemp;
                        fEnergy_Temp        += spot[i] * schrg_omp[i];
                        fEnergy_TotalCharge += schrg_omp[i];
                    }
                }

                MPI_Reduce(&fEnergy_Temp,        &fEnergy_Temp,        1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&fEnergy_TotalCharge, &fEnergy_TotalCharge, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);
            } /* end of parallel computing */
        } /* end of if (bReactFieldlnFRC || bSolvEng || bNonlinearEng || bSurfEngOut || bSurfCrgOut) */
    } /* end of splitting master and slave processes */

    /*
     * Clean memory occupied by these vectors
     */
    vector<delphi_real>().swap(prgfgSurfCrgA_nX);
    vector<delphi_real>().swap(prgfgSurfCrgA_nY);
    vector<delphi_real>().swap(prgfgSurfCrgA_nZ);
    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);

    vector<delphi_real>().swap(spdiv);
    vector<delphi_real>().swap(spot);
    vector<delphi_real>().swap(schrg_omp);
}

#endif
