/*
 * energy_run.cpp
 *
 *  Created on: Feb 23, 2014
 *      Author: Lin Wang
 */

#include "energy.h"

#ifdef PARALLEL_OMP
#include <omp.h>
#endif

void CDelphiEnergy::run()
{
    bool ido;
    int i, iGridOutput, iisitpot;
    delphi_integer iw;

    delphi_real fEnergy_Grid         = 0.0;
    delphi_real fEnergy_Solvation    = 0.0;
    delphi_real fEnergy_AnalySurf    = 0.0;
    delphi_real fEnergy_Nonlinear    = 0.0;
    delphi_real fEnergy_Coulombic    = 0.0;
    delphi_real fEnergy_AnalyGrid    = 0.0;
    delphi_real fEnergy_SolvToChgIn  = 0.0;
    delphi_real fEnergy_SolvToChgOut = 0.0; // Solvent contribution to fixed charges  outside the cube is zero.
    delphi_real fEnergy_Total        = 0.0;

    if (inhomo == 1) fIonStrength = 0;
    
    SGrid<delphi_integer> ixyz;

    infoString = " Info> ";
    timeString = " Time> ";
    enerString = " Energy> ";
    MAXWIDTH   = 45;
    NUMWIDTH   = 12;

    /*
     * -------------------- Analytic grid energy --------------------
     */
    if (bEngOut) // output energy calculation results "energy.dat"
    {
        ofstream ofEnergyFile;
        ofEnergyFile.open(strEnergyFile);
        ofEnergyFile.close();
    }

    if (bAnalyEng)
    {
        fEnergy_AnalyGrid = 0.0;
        cout << " WARNING !!! Analytic grid energy is no longer available." << endl;

        if (bEngOut)
        {
            ofstream ofEnergyFile;
            ofEnergyFile.open(strEnergyFile, std::fstream::app);
            ofEnergyFile << " Analytic grid energy is " << fEnergy_AnalyGrid << " kT.\n";
            ofEnergyFile.close();
        }
        exit (EXIT_FAILURE);
    }

    /*
     * -------------------- Total grid energy --------------------
     */
    if (bGridEng)
    {
        fEnergy_Grid = 0.0;
        lim_min = (delphi_integer)2 + ieBuffz.nMin;         //ieBuffz(pdc->getKey_constRef< SExtrema<int> >("bufz"))
        lim_max = iGrid - (delphi_integer)1 - ieBuffz.nMax; //lim%min=2+bufz%min ; lim%max=igrid-1-bufz%max

        for (i = 0; i < iCrgedGridB; i++)
        {
            ixyz = prgigGridCrgPose[i];

            ido = 1;
            if (optORLT<delphi_integer>(ixyz, lim_min) || optORGT<delphi_integer>(ixyz, lim_max)) ido = 0;

            if (ido)
            {
                iw = ((ixyz.nZ - 1) * iGrid * iGrid) + (ixyz.nY - 1) * iGrid + (ixyz.nX - 1);  // i, j, k -> k, j, i (z, y, x)
                fEnergy_Grid = fEnergy_Grid + prgfPhimap[iw] * prgfGridCrg[i];
            }
        }

        fEnergy_Grid = fEnergy_Grid / 2.0;

        if ( iGaussian == 1 && inhomo == 1 && bSolvEng )
        {
            ergsgaussian = fEnergy_Grid;
        }
        else
        {
            cout << enerString << left << setw(MAXWIDTH) << "Total grid energy" << " : " << setw(NUMWIDTH) << right
                    << setprecision(NUMPRECISION) << fixed << fEnergy_Grid << " kT" << endl;
        }

        if (bEngOut)
        {
            ofstream ofEnergyFile;
            ofEnergyFile.open(strEnergyFile, std::fstream::app);
            ofEnergyFile << "Total grid energy      :   " << setw(NUMWIDTH) << setprecision(NUMPRECISION) << fixed << fEnergy_Grid << " kT \n";
            ofEnergyFile.close();
        }

    } /* end of if (bGridEng) for calculating Total grid energy */

    if ( iGaussian == 1 && inhomo == 0 && bSolvEng)
    {
        cout << enerString << left << setw(MAXWIDTH) << "Corrected reaction field energy" << " : " << setw(NUMWIDTH)
                << right << setprecision(NUMPRECISION) << fixed << fEnergy_Grid - ergsgaussian << " kT" << endl;
        ergs = fEnergy_Grid - ergsgaussian;
    }

    /*
     * LinLi,Argo :if iGaussian==1, skip other energy terms
     */
    if (iGaussian == 0)
    {
        if (bGridEng && bAnalyEng)
        {
            #ifdef VERBOSE
            cout << " Difference energy, in kT, is  " << fEnergy_Grid-fEnergy_AnalyGrid << endl;
            cout << " Difference energy, in kcals, is  " << (fEnergy_Grid-fEnergy_AnalyGrid)*0.6 << endl;
            #endif
        }

        /*
         * For polarized membrane (not supported and not tested)
         */
        if (iBndyType == 5)
        {
            cout << " WARNING!!! Not completely tested routine for polarized membrane!!" << endl;
            exit (EXIT_FAILURE);

            if (bNonlinearEng || bAnalySurfEng)
            {
                #ifdef VERBOSE
                cout << " WARNING !!! This option is not yet working with fixed potential difference!" << endl;
                #endif

                fEnergy_AnalySurf = 0.0;
                fEnergy_Nonlinear = 0.0;
            }
        }

        /*
         * -------------------- Reaction field energy calculation --------------------
         */
        if (bReactFieldlnFRC || bSolvEng || bNonlinearEng || bAnalySurfEng || bSurfEngOut || bSurfCrgOut)
        {
            // fEnergy_SolvToChgIn = interaction energy of the solvent and the fixed charges. //
            fEnergy_Solvation   = 0.0;
            fEnergy_AnalySurf   = 0.0;
            fEnergy_Nonlinear   = 0.0;
            fEnergy_SolvToChgIn = 0.0;

            if (bPotentiallnSite) iisitpot = 1;
            else                  iisitpot = 0;

            energy_react(fEnergy_Solvation, fEnergy_AnalySurf, iisitpot); // call reaction field function
        }
    }

    //if( bCoulombEng && ( !bIonsEng || !bNonlinearEng ) )
    if (bCoulombEng && (!bIonsEng || !bNonlinearEng) && !(iGaussian == 1 && inhomo == 1))
    {
        fEnergy_Coulombic = 0.0;

        if (bIonsEng)
        {
            fEnergy_SolvToChgIn = 0.0;
            energy_clbtotal(fEnergy_SolvToChgIn, fEnergy_Coulombic); // call clbtotal function.

            {
                CIonicCalcUnderTest waring;
            }

            #ifdef VERBOSE
            cout << enerString << "Solvent contribution to fixed charges" << endl;

            cout << enerString << left << setw(MAXWIDTH) << "Respectively inside and outside the cube" << " : "
            << setw(NUMWIDTH) << right << fEnergy_SolvToChgIn << "  kT   " << setprecision(NUMPRECISION) << fixed << fEnergy_SolvToChgOut << "  kT" << endl; // where is ergestout??

            cout << enerString << left << setw(MAXWIDTH) << "Total ionic direct contribution" << " : " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << (fEnergy_SolvToChgIn+fEnergy_SolvToChgOut) << "  kT" << endl;
            #endif
        }
        else
        {
            if (iMediaNum == 1)
            {
                energy_clb(fEnergy_Coulombic);  // call clb function.
                fEnergy_Coulombic = fEnergy_Coulombic / fEpsin;
            }
            else
            {
                energy_clbmedia(fEnergy_Coulombic); // call clbmedia function.
            }
        }

        cout << enerString << left << setw(MAXWIDTH) << "Coulombic energy" << " : " << setw(NUMWIDTH) << right
            << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT" << endl;

        if (bEngOut)
        {
            ofstream ofEnergyFile;
            ofEnergyFile.open(strEnergyFile, std::fstream::app);
            ofEnergyFile << "Total coulombic energy  :  " << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT \n";
            ofEnergyFile.close();
        }
    }

    if (bNonlinearEng)
    {
        energy_nonl(fEnergy_Nonlinear, iGridOutput);// call nonlinear function.

        if (bIonsEng)
        {
            fEnergy_Coulombic = 0.0;
            fEnergy_SolvToChgIn = 0.0;
            energy_clbnonl(fEnergy_Coulombic, fEnergy_SolvToChgIn, iGridOutput);
            cout << enerString << left << setw(MAXWIDTH) << "Direct ionic contribution inside the box" << " : "
                << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_SolvToChgIn << " kT" << endl;

            cout << enerString << left << setw(MAXWIDTH) << "Coulombic energy" << " : " << setw(NUMWIDTH) << right
                << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT" << endl;

            if (bEngOut)
            {
                ofstream ofEnergyFile;
                ofEnergyFile.open(strEnergyFile, std::fstream::app);
                ofEnergyFile << "Total coulombic energy :  " << setw(NUMWIDTH) << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT \n";
                ofEnergyFile.close();
            }
        }

        vector<SGridValue<delphi_real> >().swap(sout);
    }

    if (iGaussian == 0) //if iGaussian==1 skip these terms
    {
        if (bSolvEng && bIonsEng)
            cout << " Energy arising from solvent and boundary pol.  " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION)
                << fixed << (fEnergy_Nonlinear + fEnergy_Solvation + fEnergy_SolvToChgIn + fEnergy_SolvToChgOut) << " kT" << "\n";

        if (bNonlinearEng && bGridEng)
            cout << enerString << left << setw(MAXWIDTH) << "Total non linear grid energy" << " : " << setw(NUMWIDTH)
                << right << setprecision(NUMPRECISION) << fixed << (fEnergy_Grid + fEnergy_Nonlinear) << " kT \n";

        fEnergy_Total = fEnergy_Nonlinear + fEnergy_Coulombic + fEnergy_Solvation + fEnergy_SolvToChgIn + fEnergy_SolvToChgOut;

        if (bSolvEng || bCoulombEng)
        {
            cout << enerString << left << setw(MAXWIDTH) << "All required energy terms but grid energy    " << " : "
                << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_Total << " kT" << endl;

            if (bEngOut)
            {
                ofstream ofEnergyFile;
                ofEnergyFile.open(strEnergyFile, std::fstream::app);
                ofEnergyFile << "Total required energy (everything calculated but grid and self_reaction energies: "
                    << setw(NUMWIDTH) << setprecision(NUMPRECISION) << fixed << fEnergy_Total << " kT \n";
                ofEnergyFile.close();
            }
        }

        if (bAnalySurfEng && bAnalyEng && bGridEng)
            cout << enerString << left << setw(MAXWIDTH) << "Excess grid energy" << " : " << setw(NUMWIDTH) << right
                << setprecision(NUMPRECISION) << fixed << (fEnergy_Grid - fEnergy_AnalySurf - fEnergy_AnalyGrid) << endl;

        cout << endl;

        /*
         * Write into data container
         */
        ergg    = fEnergy_Grid;
        ergc    = fEnergy_Coulombic;
        ergs    = fEnergy_Solvation;
        ergions = fEnergy_SolvToChgIn + fEnergy_SolvToChgOut;
    }
    else if (iGaussian == 1)
    {
        ergg = fEnergy_Grid;
        ergc = fEnergy_Coulombic;
        ergs = fEnergy_Solvation;
    }
}

#ifdef PARALLEL_MPI

/*
 * implementation of pure abstract function mpi_run() in class IAbstractModule
 */
void CDelphiEnergy::mpi_run()
{
    bool ido;
    int i, iGridOutput, iisitpot;
    delphi_integer iw;

    delphi_real fEnergy_Grid         = 0.0;
    delphi_real fEnergy_Solvation    = 0.0;
    delphi_real fEnergy_AnalySurf    = 0.0;
    delphi_real fEnergy_Nonlinear    = 0.0;
    delphi_real fEnergy_Coulombic    = 0.0;
    delphi_real fEnergy_AnalyGrid    = 0.0;
    delphi_real fEnergy_SolvToChgIn  = 0.0;
    delphi_real fEnergy_SolvToChgOut = 0.0; // Solvent contribution to fixed charges outside the cube is zero.
    delphi_real fEnergy_Total        = 0.0;

    if (inhomo == 1) fIonStrength = 0;
    
    SGrid<delphi_integer> ixyz;

    infoString = " Info> ";
    timeString = " Time> ";
    enerString = " Energy> ";
    MAXWIDTH   = 45;
    NUMWIDTH   = 12;
    
    //cout << "rank " << mpi_rank << ": in enery mpi_run" << endl;   
    
    if (0 == mpi_rank) /* master process */
    {
        /*
         * -------------------- Analytic grid energy --------------------
         */
        if (bEngOut) // output energy calculation results "energy.dat"
        {
            ofstream ofEnergyFile;
            ofEnergyFile.open(strEnergyFile);
            ofEnergyFile.close();
        }

        if (bAnalyEng)
        {
            fEnergy_AnalyGrid = 0.0;
            cout << " WARNING !!! Analytic grid energy is no longer available." << endl;

            if (bEngOut)
            {
                ofstream ofEnergyFile;
                ofEnergyFile.open(strEnergyFile, std::fstream::app);
                ofEnergyFile << " Analytic grid energy is " << fEnergy_AnalyGrid << " kT.\n";
                ofEnergyFile.close();
            }
            exit (EXIT_FAILURE);
        }

        /*
         * -------------------- parallelize Total grid energy calculations --------------------
         */
        if (bGridEng)
        {
            fEnergy_Grid = 0.0;
            lim_min = (delphi_integer)2 + ieBuffz.nMin;         //ieBuffz(pdc->getKey_constRef< SExtrema<int> >("bufz"))
            lim_max = iGrid - (delphi_integer)1 - ieBuffz.nMax; //lim%min=2+bufz%min ; lim%max=igrid-1-bufz%max

            MPI_Reduce(MPI_IN_PLACE, &fEnergy_Grid, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);

            fEnergy_Grid = fEnergy_Grid / 2.0;

            if ( iGaussian == 1 && inhomo == 1 && bSolvEng )
            {
                ergsgaussian = fEnergy_Grid;
            }
            else
            {
                cout << enerString << left << setw(MAXWIDTH) << "Total grid energy" << " : " << setw(NUMWIDTH) << right
                        << setprecision(NUMPRECISION) << fixed << fEnergy_Grid << " kT" << endl;
            }

            if (bEngOut)
            {
                ofstream ofEnergyFile;
                ofEnergyFile.open(strEnergyFile, std::fstream::app);
                ofEnergyFile << "Total grid energy      :   " << setw(NUMWIDTH) << setprecision(NUMPRECISION) << fixed << fEnergy_Grid << " kT \n";
                ofEnergyFile.close();
            }

        } /* end of if (bGridEng) for calculating otal grid energy */

        if ( iGaussian == 1 && inhomo == 0 && bSolvEng)
        {
            cout << enerString << left << setw(MAXWIDTH) << "Corrected reaction field energy" << " : " << setw(NUMWIDTH)
                    << right << setprecision(NUMPRECISION) << fixed << fEnergy_Grid - ergsgaussian << " kT" << endl;
            ergs = fEnergy_Grid - ergsgaussian;
        }

        /*
         * LinLi,Argo :if iGaussian==1, skip other energy terms
         */
        if (iGaussian == 0)
        {
            if (bGridEng && bAnalyEng)
            {
                #ifdef VERBOSE
                cout << " Difference energy, in kT, is  " << fEnergy_Grid-fEnergy_AnalyGrid << endl;
                cout << " Difference energy, in kcals, is  " << (fEnergy_Grid-fEnergy_AnalyGrid)*0.6 << endl;
                #endif
            }

            /*
             * For polarized membrane (not supported and not tested)
             */
            if (iBndyType == 5)
            {
                cout << " WARNING!!! Not completely tested routine for polarized membrane!!" << endl;
                exit (EXIT_FAILURE);

                if (bNonlinearEng || bAnalySurfEng)
                {
                    #ifdef VERBOSE
                    cout << " WARNING !!! This option is not yet working with fixed potential difference!" << endl;
                    #endif

                    fEnergy_AnalySurf = 0.0;
                    fEnergy_Nonlinear = 0.0;
                }
            }

            /*
             * -------------------- parallelize Reaction field energy calculation --------------------
             */
            if (bReactFieldlnFRC || bSolvEng || bNonlinearEng || bAnalySurfEng || bSurfEngOut || bSurfCrgOut)
            {
                // fEnergy_SolvToChgIn = interaction energy of the solvent and the fixed charges. //
                fEnergy_Solvation   = 0.0;
                fEnergy_AnalySurf   = 0.0;
                fEnergy_Nonlinear   = 0.0;
                fEnergy_SolvToChgIn = 0.0;

                if (bPotentiallnSite) iisitpot = 1;
                else                  iisitpot = 0;

                /*
                 * -------------------- parallel computing starts --------------------
                 */
                mpi_energy_react(fEnergy_Solvation, fEnergy_AnalySurf, iisitpot); // call reaction field function
            }
        }

        //if( bCoulombEng && ( !bIonsEng || !bNonlinearEng ) )
        if (bCoulombEng && (!bIonsEng || !bNonlinearEng) && !(iGaussian == 1 && inhomo == 1))
        {
            fEnergy_Coulombic = 0.0;

            if (bIonsEng)
            {
                fEnergy_SolvToChgIn = 0.0;

                /*
                 * -------------------- parallel computing starts --------------------
                 */
                mpi_energy_clbtotal(fEnergy_SolvToChgIn, fEnergy_Coulombic); // call clbtotal function.

                {
                    CIonicCalcUnderTest waring;
                }

                #ifdef VERBOSE
                cout << enerString << "Solvent contribution to fixed charges" << endl;

                cout << enerString << left << setw(MAXWIDTH) << "Respectively inside and outside the cube" << " : "
                << setw(NUMWIDTH) << right << fEnergy_SolvToChgIn << "  kT   " << setprecision(NUMPRECISION) << fixed << fEnergy_SolvToChgOut << "  kT" << endl; // where is ergestout??

                cout << enerString << left << setw(MAXWIDTH) << "Total ionic direct contribution" << " : " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << (fEnergy_SolvToChgIn+fEnergy_SolvToChgOut) << "  kT" << endl;
                #endif
            }
            else
            {
                /*
                 * -------------------- parallel computing starts --------------------
                 */
                if (iMediaNum == 1)
                {
                    mpi_energy_clb(fEnergy_Coulombic);  // call clb function.
                    fEnergy_Coulombic = fEnergy_Coulombic / fEpsin;
                }
                else
                {
                    mpi_energy_clbmedia(fEnergy_Coulombic); // call clbmedia function.
                }
            }

            cout << enerString << left << setw(MAXWIDTH) << "Coulombic energy" << " : " << setw(NUMWIDTH) << right
                << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT" << endl;

            if (bEngOut)
            {
                ofstream ofEnergyFile;
                ofEnergyFile.open(strEnergyFile, std::fstream::app);
                ofEnergyFile << "Total coulombic energy  :  " << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT \n";
                ofEnergyFile.close();
            }
        }

        if (bNonlinearEng)
        {
            energy_nonl(fEnergy_Nonlinear, iGridOutput);// call nonlinear function.

            if (bIonsEng)
            {
                fEnergy_Coulombic   = 0.0;
                fEnergy_SolvToChgIn = 0.0;

                /*
                 * -------------------- parallel computing starts --------------------
                 */
                MPI_Bcast(&iGridOutput, 1, MPI_INT, 0, MPI_COMM_WORLD);
                mpi_energy_clbnonl(fEnergy_Coulombic, fEnergy_SolvToChgIn, iGridOutput);

                cout << enerString << left << setw(MAXWIDTH) << "Direct ionic contribution inside the box" << " : "
                    << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_SolvToChgIn << " kT" << endl;

                cout << enerString << left << setw(MAXWIDTH) << "Coulombic energy" << " : " << setw(NUMWIDTH) << right
                    << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT" << endl;

                if (bEngOut)
                {
                    ofstream ofEnergyFile;
                    ofEnergyFile.open(strEnergyFile, std::fstream::app);
                    ofEnergyFile << "Total coulombic energy :  " << setw(NUMWIDTH) << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT \n";
                    ofEnergyFile.close();
                }
            }

            vector<SGridValue<delphi_real> >().swap(sout);
        }

        if (iGaussian == 0) //if iGaussian==1 skip these terms
        {
            if (bSolvEng && bIonsEng)
                cout << " Energy arising from solvent and boundary pol.  " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION)
                    << fixed << (fEnergy_Nonlinear + fEnergy_Solvation + fEnergy_SolvToChgIn + fEnergy_SolvToChgOut) << " kT" << "\n";

            if (bNonlinearEng && bGridEng)
                cout << enerString << left << setw(MAXWIDTH) << "Total non linear grid energy" << " : " << setw(NUMWIDTH)
                    << right << setprecision(NUMPRECISION) << fixed << (fEnergy_Grid + fEnergy_Nonlinear) << " kT \n";

            fEnergy_Total = fEnergy_Nonlinear + fEnergy_Coulombic + fEnergy_Solvation + fEnergy_SolvToChgIn + fEnergy_SolvToChgOut;

            if (bSolvEng || bCoulombEng)
            {
                cout << enerString << left << setw(MAXWIDTH) << "All required energy terms but grid energy    " << " : "
                    << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_Total << " kT" << endl;

                if (bEngOut)
                {
                    ofstream ofEnergyFile;
                    ofEnergyFile.open(strEnergyFile, std::fstream::app);
                    ofEnergyFile << "Total required energy (everything calculated but grid and self_reaction energies: "
                        << setw(NUMWIDTH) << setprecision(NUMPRECISION) << fixed << fEnergy_Total << " kT \n";
                    ofEnergyFile.close();
                }
            }

            if (bAnalySurfEng && bAnalyEng && bGridEng)
                cout << enerString << left << setw(MAXWIDTH) << "Excess grid energy" << " : " << setw(NUMWIDTH) << right
                    << setprecision(NUMPRECISION) << fixed << (fEnergy_Grid - fEnergy_AnalySurf - fEnergy_AnalyGrid) << endl;

            cout << endl;

            /*
             * Write into data container
             */
            ergg    = fEnergy_Grid;
            ergc    = fEnergy_Coulombic;
            ergs    = fEnergy_Solvation;
            ergions = fEnergy_SolvToChgIn + fEnergy_SolvToChgOut;
        }
        else if (iGaussian == 1)
        {
            ergg = fEnergy_Grid;
            ergc = fEnergy_Coulombic;
            ergs = fEnergy_Solvation;
        }
    }
    else /* slave processes */
    {
        /*
         * -------------------- parallelize Total grid energy calculations --------------------
         */
        if (bGridEng)
        {
            fEnergy_Grid = 0.0;
            lim_min = 2 + ieBuffz.nMin; //ieBuffz(pdc->getKey_constRef< SExtrema<int> >("bufz"))
            lim_max = iGrid - 1 - ieBuffz.nMax; //lim%min=2+bufz%min ; lim%max=igrid-1-bufz%max

            if ( iCrgedGridB <= mpi_num_workers )
            {
                if (mpi_rank <= iCrgedGridB)
                {
                    i = mpi_rank - 1;
                    ixyz = prgigGridCrgPose[i];
                    ido = 1;
                    if (optORLT<delphi_integer>(ixyz, lim_min) || optORGT<delphi_integer>(ixyz, lim_max)) ido = 0;
                    if (ido)
                    {
                        iw = ((ixyz.nZ - 1) * iGrid * iGrid) + (ixyz.nY - 1) * iGrid + (ixyz.nX - 1);  // i, j, k -> k, j, i (z, y, x)
                        fEnergy_Grid = fEnergy_Grid + prgfPhimap[iw] * prgfGridCrg[i];
                    }
                }
            }
            /*
             * a simple parallelization but does NOT produce high speedups
             */
//            else
//            {
//                int increment = mpi_num_workers;
//                for (i = mpi_rank - 1; i < iCrgedGridB; i += increment)
//                {
//                    ixyz = prgigGridCrgPose[i];
//                    ido = 1;
//                    if (optORLT<delphi_integer>(ixyz, lim_min) || optORGT<delphi_integer>(ixyz, lim_max)) ido = 0;
//                    if (ido)
//                    {
//                        iw = ((ixyz.nZ - 1) * iGrid * iGrid) + (ixyz.nY - 1) * iGrid + (ixyz.nX - 1);  // i, j, k -> k, j, i (z, y, x)
//                        fEnergy_Grid = fEnergy_Grid + prgfPhimap[iw] * prgfGridCrg[i];
//                    }
//                }
//            }
            /*
             * another parallelization to improve performance by using cache more effectively
             */
            else
            {
                int mpi_mystart, mpi_myend, mpi_split = iCrgedGridB % mpi_num_workers;

                if (mpi_rank <= mpi_split)
                {
                    mpi_mystart = (mpi_rank - 1) * (iCrgedGridB / mpi_num_workers + 1);
                    mpi_myend   = mpi_mystart + (iCrgedGridB / mpi_num_workers + 1) - 1;
                }
                else
                {
                    mpi_mystart = mpi_split * (iCrgedGridB / mpi_num_workers + 1) + (mpi_rank - mpi_split - 1) * (iCrgedGridB / mpi_num_workers);
                    mpi_myend   = mpi_mystart + (iCrgedGridB / mpi_num_workers) - 1;
                }

                if (mpi_myend >= iCrgedGridB - 1) mpi_myend = iCrgedGridB - 1;

                for (i = mpi_mystart; i <= mpi_myend; i++)
                {
                    ixyz = prgigGridCrgPose[i];
                    ido = 1;
                    if (optORLT<delphi_integer>(ixyz, lim_min) || optORGT<delphi_integer>(ixyz, lim_max)) ido = 0;
                    if (ido)
                    {
                        iw = ((ixyz.nZ - 1) * iGrid * iGrid) + (ixyz.nY - 1) * iGrid + (ixyz.nX - 1);  // i, j, k -> k, j, i (z, y, x)
                        fEnergy_Grid = fEnergy_Grid + prgfPhimap[iw] * prgfGridCrg[i];
                    }
                }
            }

            MPI_Reduce(&fEnergy_Grid, &fEnergy_Grid, 1, mpi_delphi_real, MPI_SUM, 0, MPI_COMM_WORLD);
        } /* end of if (bGridEng) for calculating otal grid energy */

        if (iGaussian == 0) //if iGaussian==1 skip these terms
        {
            /*
             * -------------------- parallelize Reaction field energy calculation --------------------
             */
            if (bReactFieldlnFRC || bSolvEng || bNonlinearEng || bAnalySurfEng || bSurfEngOut || bSurfCrgOut)
            {
                // fEnergy_SolvToChgIn = interaction energy of the solvent and the fixed charges. //
                fEnergy_Solvation   = 0.0;
                fEnergy_AnalySurf   = 0.0;
                fEnergy_Nonlinear   = 0.0;
                fEnergy_SolvToChgIn = 0.0;

                if (bPotentiallnSite) iisitpot = 1;
                else                  iisitpot = 0;

                mpi_energy_react(fEnergy_Solvation, fEnergy_AnalySurf, iisitpot); // call reaction field function
            }
        }

        if (bCoulombEng && (!bIonsEng || !bNonlinearEng) && !(iGaussian == 1 && inhomo == 1))
        {
            fEnergy_Coulombic = 0.0;

            if (bIonsEng)
            {
                fEnergy_SolvToChgIn = 0.0;
                mpi_energy_clbtotal(fEnergy_SolvToChgIn, fEnergy_Coulombic); // call clbtotal function.
            }
            else
            {
                if (iMediaNum == 1)
                {
                    mpi_energy_clb(fEnergy_Coulombic);  // call clb function.
                }
                else
                {
                    mpi_energy_clbmedia(fEnergy_Coulombic); // call clbmedia function.
                }
            }
        }

        if (bNonlinearEng)
        {
            if (bIonsEng)
            {
                fEnergy_Coulombic   = 0.0;
                fEnergy_SolvToChgIn = 0.0;

                MPI_Bcast(&iGridOutput, 1, MPI_INT, 0, MPI_COMM_WORLD);
                mpi_energy_clbnonl(fEnergy_Coulombic, fEnergy_SolvToChgIn, iGridOutput);
            }
        }


    } /* end of splitting master and slave processes */
}

#endif
