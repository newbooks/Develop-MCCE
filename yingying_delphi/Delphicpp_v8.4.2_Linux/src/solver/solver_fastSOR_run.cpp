/*
 * solver_fastSOR_run.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

void CDelphiFastSOR::run() 
{
    debug_solver = false;
    const char * infoString = " Info> ";
    const char * timeString = " Time> ";
    const char * bndcString = " Bndc> ";
    size_t MAXWIDTH = 45;

    if (inhomo == 1) fIonStrength = 0;

    if (iGaussian == 1 && inhomo == 0) //reset some vectors for 2nd run of Gaussian
    {
        if (prgfGridCrg.size() > 0) 
        {
            if (debug_solver) cout << "prgfGridCrg.size(): " << prgfGridCrg.size() << endl;
            vector<delphi_real>().swap(prgfGridCrg); //gchrg
        }
        if (prgigGridCrgPose.size() > 0) 
        {
            if (debug_solver) cout << "prgigGridCrgPose.size(): " << prgigGridCrgPose.size() << endl;
            vector<SGrid<delphi_integer> >().swap(prgigGridCrgPose); //gchrgp
        }

        //vector<delphi_integer> ibndx,ibndy,ibndz;
        if (ibndx.size() > 0) 
        {
            if (debug_solver) cout << "ibndx.size(): " << ibndx.size() << endl;
            vector<delphi_integer>().swap(ibndx); //ibndx
        }
        if (ibndy.size() > 0) 
        {
            if (debug_solver) cout << "ibndy.size(): " << ibndy.size() << endl;
            vector<delphi_integer>().swap(ibndy); //ibndy
        }
        if (ibndz.size() > 0) 
        {
            if (debug_solver) cout << "ibndz.size(): " << ibndz.size() << endl;
            vector<delphi_integer>().swap(ibndz); //ibndz
        }
        if (phimap1.size() > 0) 
        {
            if (debug_solver) cout << "phimap1.size(): " << phimap1.size() << endl;
            vector<delphi_real>().swap(phimap1); //phimap1
        }
        if (phimap2.size() > 0) 
        {
            if (debug_solver) cout << "phimap2.size(): " << phimap2.size() << endl;
            vector<delphi_real>().swap(phimap2); //phimap2
        }

        if (prgiBndyDielecIndex.size() > 0) 
        {
            if (debug_solver) cout << "prgiBndyDielecIndex.size(): " << prgiBndyDielecIndex.size() << endl;
            vector<delphi_integer>().swap(prgiBndyDielecIndex); //prgiBndyDielecIndex
        }
        if (prgfBndyDielec.size() > 0) 
        {
            if (debug_solver) cout << "prgfBndyDielec.size(): " << prgfBndyDielec.size() << endl;
            vector<vector<delphi_real> >().swap(prgfBndyDielec); //prgfBndyDielec
        }
        //vector< SDoubleGridValue >& prgdgvCrgBndyGrid;   // cgbp(ibc)
        if (prgdgvCrgBndyGrid.size() > 0) 
        {
            if (debug_solver) cout << "prgdgvCrgBndyGrid.size(): " << prgdgvCrgBndyGrid.size() << endl;
            vector<SDoubleGridValue>().swap(prgdgvCrgBndyGrid); //prgdgvCrgBndyGrid
        }
        if (debug_solver) cout << "prgfPhiMap.size(): " << prgfPhiMap.size() << endl;

        if (gaussianBoundaryDensity.size() > 0) 
        {
            if (debug_solver) cout << "gaussianBoundaryDensity.size(): " << gaussianBoundaryDensity.size() << endl;
            vector<delphi_real>().swap(gaussianBoundaryDensity); //prgdgvCrgBndyGrid
        }

        iDielecBndyOdd = 0;

    }
    
    if (debug_solver) cout << "iDielecBndyOdd: " << iDielecBndyOdd << endl;

    validateInput();
    setDielecBndySaltMap();
    setCrg();
    
    cout << timeString << "iepsmp to db, and charging done on ";
    pTimer->showTime();
    cout << endl;
    cout << left << setw(MAXWIDTH) << " Number of grid points assigned charge" << " : " << iCrgedGridSum << endl;

    if (bEpsOut && iGaussian == 0) //----- write dielectric map
    {
        unique_ptr < CIO > pio(new CIO()); // smart unique_ptr
        // pio->writeEpsMap(iAtomNum,iObjectNum,iGrid,fScale,fgBoxCenter,prgigEpsMap,prgbDielecMap,strEpsFile);
        pio->writeHomoEpsMap(iGrid, repsout, repsin, fScale, fgBoxCenter, prgbDielecMap, strEpsFile);
        pio.reset();
    }

    if (bEpsOut && iGaussian == 1) //----- write dielectric map
    {
        unique_ptr < CIO > pio(new CIO()); // smart unique_ptr
        // pio->writeEpsMap(iAtomNum,iObjectNum,iGrid,fScale,fgBoxCenter,prgigEpsMap,prgbDielecMap,strEpsFile);
        pio->writeGaussEpsMap(iGrid, repsout, fScale, fgBoxCenter, gepsmp2, strEpsFile);
        pio.reset();
    }

    //----- the order of calculateRelaxFactor and setBndy cannnot be reverted! prgfPhiMap is used
    //      as temporary array in calculateRelaxFactor...

    phimap1.assign((iGrid*iGrid*iGrid + 1) / 2, 0.0); //phimap1.assign(icgrid/2,0.0);
    phimap2.assign((iGrid*iGrid*iGrid + 1) / 2, 0.0); //phimap2.assign(icgrid/2,0.0);

    if (bSpectralRadius) 
    {
        fSpec = fSpectralRadius;
        cout << infoString << left << setw(MAXWIDTH) << "Using entered value for relaxation of" << " : " << fSpec << endl;
    } 
    else 
    {
        fSpec = calculateRelaxFactor();
    }

    int noit = (int) (7.8 / log(1.0 + sqrt(1 - fSpec)));

    cout << left << setw(MAXWIDTH) << " Estimated iterations to convergence" << " : " << noit << endl;

    if (bAutoConverge) iLinIterateNum = noit;

    /*
     * setBndy() requires using phimap. move allocating phimap here to save memory usage
     */
    prgfPhiMap.assign(iGrid * iGrid * iGrid, 0.0);

    setBndy();
    
    cout << endl;
    cout << timeString << "Setup time done on ";
    pTimer->showTime();
    cout << endl;
    
    if (0.3 < (delphi_real) iCrgBndyGridNum / iBndyGridNum && bAutoConverge) 
    {
        iLinIterateNum = iLinIterateNum * iCrgBndyGridNum / (0.3 * iBndyGridNum);
        cout << left << setw(MAXWIDTH) << " Re-estimated iterations now " << " : " << iLinIterateNum << endl;
    }

    cout << timeString << "Now iterating on ";
    pTimer->showTime();
    cout << endl;

    if (0 == iNonIterateNum || fZero > fIonStrength || inhomo == 1) 
    {
        if (0 >= iLinIterateNum) throw CZeorLinIterateNum(bAutoConverge, iLinIterateNum);
        itit();
    } 
    else 
    {
        if (50 < noit || 0 >= iLinIterateNum) 
        {
            iLinIterateNum = noit / 2;
            cout << left << setw(MAXWIDTH) << " Re-estimated iterations now " << " : " << iLinIterateNum << endl;
        }

        delphi_real qfact = abs(fNetCrg) * (delphi_real) iCrgBndyGridNum / iBndyGridNum;

        nitit(qfact);
    }
}

#ifdef PARALLEL_MPI

/*
 * implementation of pure abstract function mpi_run() in class IAbstractModule
 */
void CDelphiFastSOR::mpi_run()
{
    int noit;

    if (inhomo == 1) fIonStrength = 0;
	
   /*
    * master process carries out necessary preparations for iterations
    */
   if (0 == mpi_rank)
   {
      debug_solver = false;
      const char * infoString = " Info> ";
      const char * timeString = " Time> ";
      const char * bndcString = " Bndc> ";
      size_t MAXWIDTH = 45;
      
      if( iGaussian == 1 && inhomo == 0 ) //reset some vectors for 2nd run of Gaussian
      {
         if ( prgfGridCrg.size() > 0 )
         {
            if(debug_solver) cout << "prgfGridCrg.size(): " << prgfGridCrg.size() << endl;
            vector<delphi_real>().swap(prgfGridCrg) ;//gchrg
         }

         if ( prgigGridCrgPose.size() > 0 )
         {
            if(debug_solver) cout << "prgigGridCrgPose.size(): " << prgigGridCrgPose.size() << endl;
            vector< SGrid<delphi_integer> >().swap(prgigGridCrgPose) ; //gchrgp
         }

         //vector<delphi_integer> ibndx,ibndy,ibndz;
         if ( ibndx.size() > 0 )
         {
            if(debug_solver) cout << "ibndx.size(): " << ibndx.size() << endl;
            vector<delphi_integer>().swap(ibndx) ; //ibndx
         }

         if ( ibndy.size() > 0 )
         {
            if(debug_solver) cout << "ibndy.size(): " << ibndy.size() << endl;
            vector<delphi_integer>().swap(ibndy) ; //ibndy
         }

         if ( ibndz.size() > 0 )
         {
            if(debug_solver) cout << "ibndz.size(): " << ibndz.size() << endl;
            vector<delphi_integer>().swap(ibndz) ; //ibndz
         }

         if ( phimap1.size() > 0 )
         {
            if(debug_solver) cout << "phimap1.size(): " << phimap1.size() << endl;
            vector<delphi_real>().swap(phimap1) ; //phimap1
         }

         if ( phimap2.size() > 0 )
         {
            if(debug_solver) cout << "phimap2.size(): " << phimap2.size() << endl;
            vector<delphi_real>().swap(phimap2) ; //phimap2
         }

         if ( prgiBndyDielecIndex.size() > 0 )
         {
            if(debug_solver) cout << "prgiBndyDielecIndex.size(): " << prgiBndyDielecIndex.size() << endl;
            vector<delphi_integer>().swap(prgiBndyDielecIndex) ; //prgiBndyDielecIndex
         }

         if ( prgfBndyDielec.size() > 0 )
         {
            if(debug_solver) cout << "prgfBndyDielec.size(): " << prgfBndyDielec.size() << endl;
            vector< vector<delphi_real> >().swap(prgfBndyDielec) ; //prgfBndyDielec
         }

         //vector< SDoubleGridValue >& prgdgvCrgBndyGrid;   // cgbp(ibc)
         if ( prgdgvCrgBndyGrid.size() > 0 )
         {
            if(debug_solver) cout << "prgdgvCrgBndyGrid.size(): " << prgdgvCrgBndyGrid.size() << endl;
            vector< SDoubleGridValue >().swap(prgdgvCrgBndyGrid) ; //prgdgvCrgBndyGrid
         }

         if(debug_solver) cout << "prgfPhiMap.size(): " << prgfPhiMap.size() << endl;

         if (gaussianBoundaryDensity.size() > 0) 
         {
             if (debug_solver) cout << "gaussianBoundaryDensity.size(): " << gaussianBoundaryDensity.size() << endl;
             vector<delphi_real>().swap(gaussianBoundaryDensity); //prgdgvCrgBndyGrid
         }         
         
         iDielecBndyOdd=0;
      }

      if(debug_solver) cout << "iDielecBndyOdd: " << iDielecBndyOdd << endl;

      validateInput();

      /*
       * master syncs prgbDielecMap (a.k.a idebmap(igrid,igrid,igrid)) from global data container
       */
      pdc->readGlobalVector1D<char>("idebmap", 0, iGrid * iGrid * iGrid, prgbDielecMap);

      mpi_setDielecBndySaltMap();

      mpi_setCrg();

      /*
       * master frees prgbDielecMap (a.k.a idebmap(igrid,igrid,igrid)) in local data container
       * when running linear solver
       */
      if (0 == iNonIterateNum || fZero > fIonStrength || inhomo == 1)
          if (prgbDielecMap.size() > 0) vector<char>().swap(prgbDielecMap);

      cout << timeString << "iepsmp to db, and charging done on ";
      pTimer->showTime();
      cout << endl;
      cout << left << setw(MAXWIDTH) << " Number of grid points assigned charge" << " : " << iCrgedGridSum << endl;

      /*
       * MPI implementation does not output dielectric maps
       */
//      if (bEpsOut && iGaussian == 0) //----- write dielectric map
//      {
//         unique_ptr<CIO> pio(new CIO()); // smart unique_ptr
//         pio->writeHomoEpsMap(iGrid,repsout, repsin,fScale,fgBoxCenter,prgbDielecMap,strEpsFile);
//         pio.reset();
//      }
//
//      if (bEpsOut && iGaussian==1) //----- write dielectric map
//      {
//         unique_ptr<CIO> pio(new CIO()); // smart unique_ptr
//         pio->writeGaussEpsMap(iGrid, repsout, fScale, fgBoxCenter, gepsmp2, strEpsFile);
//         pio.reset();
//      }

      //----- the order of calculateRelaxFactor and setBndy cannot be reverted! prgfPhiMap is used
      //      as temporary array in calculateRelaxFactor...
      if (bSpectralRadius)
      {
         fSpec = fSpectralRadius;
         cout << infoString << left << setw(MAXWIDTH) << "Using entered value for relaxation of" << " : " << fSpec << endl;
      }
      else
      {
         fSpec = calculateRelaxFactor();
      }

      noit = (int)(7.8 / log(1.0 + sqrt(1 - fSpec)));

      cout << left << setw(MAXWIDTH) << " Estimated iterations to convergence" << " : " << noit << endl;

      if (bAutoConverge) iLinIterateNum = noit;

      /*
       * setBndy() requires using phimap. move allocating phimap here to save memory usage
       */
      prgfPhiMap.assign(iGrid * iGrid * iGrid, 0.0);

      /*
       * master syncs prggvCrgedAtom (a.k.a atmcrg(nqass)) from global data container if iBndyType = 4
       */
      if ( (7 == iBndyType && !(bIonsEng && 0 == iNonIterateNum && fZero < abs(fNetCrg)) ) || 4 == iBndyType )
          pdc->readGlobalVector1D< SGridValue<delphi_real> >("atmcrg", 0, iCrgGridNum, prggvCrgedAtom);

      setBndy();

      /*
       * master frees prggvCrgedAtom (a.k.a atmcrg(nqass)) in local data container to save memory usage
       */
      if (prggvCrgedAtom.size() > 0) vector< SGridValue<delphi_real> >().swap(prggvCrgedAtom);

      cout << endl;
      cout << timeString << "Setup done on ";
      pTimer->showTime();
      cout << endl;

      if (0.3 < (delphi_real)iCrgBndyGridNum / iBndyGridNum && bAutoConverge)
      {
         iLinIterateNum = iLinIterateNum * iCrgBndyGridNum / (0.3 * iBndyGridNum);
         cout << left << setw(MAXWIDTH) << " Re-estimated iterations now " << " : " << iLinIterateNum << endl;
      }

      cout << timeString << "Now iterating on ";
      pTimer->showTime();
      cout << endl;

      if (0 == iNonIterateNum || fZero > fIonStrength || inhomo == 1)
      {
    	  	  if (0 >= iLinIterateNum) throw CZeorLinIterateNum(bAutoConverge,iLinIterateNum);
      }
      else
      {
          if (50 < noit || 0 >= iLinIterateNum)
          {
             iLinIterateNum = noit / 2;
             cout << left << setw(MAXWIDTH) << " Re-estimated iterations now " << " : " << iLinIterateNum << endl;
          }
      }
   } // ----- end of if (0 == mpi_rank())

   /*
    * ---------- the following lines are performed on both master and salve processes ----------
    */

   if (iGaussian == 1) // Gaussian based runs
   {
       /*
        * Gaussian is not parallelized yet!
        */
       if (0 == mpi_rank)
       {
           if (0 == iNonIterateNum || fZero > fIonStrength || inhomo == 1)
               itit();
           else
           {
               delphi_real qfact = abs(fNetCrg) * (delphi_real) iCrgBndyGridNum / iBndyGridNum;
               nitit(qfact);
           }
       }
   }
   else // regular runs
   {
	   /*
	    * data container variables not to be changed
	    */
	   MPI_Bcast(&fIonStrength,     1, mpi_delphi_real, 0, MPI_COMM_WORLD); // rionst    (value from data container)
	   
	   /*
        * reference to read-and-write variables from data container
        */
	   MPI_Bcast(&iGrid,            1, MPI_INT,         0, MPI_COMM_WORLD); // igrid     (reference to data container)
	   MPI_Bcast(&iBndyType,        1, MPI_INT,         0, MPI_COMM_WORLD); // ibctyp    (reference to data container)
	   MPI_Bcast(&fNetCrg,          1, mpi_delphi_real, 0, MPI_COMM_WORLD); // qnet      (reference to data container)
	   MPI_Bcast(&iLinIterateNum,   1, MPI_INT,         0, MPI_COMM_WORLD); // nlit      (reference to data container)
	   MPI_Bcast(&iIterateInterval, 1, MPI_INT,         0, MPI_COMM_WORLD); // icon1     (reference to data container)
	   MPI_Bcast(&iConvergeFract,   1, MPI_INT,         0, MPI_COMM_WORLD); // icon2     (reference to data container)
	   MPI_Bcast(&fRelaxParam,      1, mpi_delphi_real, 0, MPI_COMM_WORLD); // relpar    (reference to data container)
	   
	   /*
	    * out into data container
	    */
	   MPI_Bcast(&iDielecBndyOdd,   1, MPI_INT,         0, MPI_COMM_WORLD); // icount2b  (reference to data container)
	   MPI_Bcast(&iCrgedGridSum,    1, MPI_INT,         0, MPI_COMM_WORLD); // icount1b  (reference to data container)
	   MPI_Bcast(&iCrgBndyGridNum,  1, MPI_INT,         0, MPI_COMM_WORLD); // ibc       (reference to data container)
	   
	   /*
	    * class variables
	    */
	   MPI_Bcast(&fSpec,            1, mpi_delphi_real, 0, MPI_COMM_WORLD); // fSpec     (local value)
	   MPI_Bcast(&iDielecBndyEven,  1, MPI_INT,         0, MPI_COMM_WORLD); // icount2a  (local value)
	   MPI_Bcast(&iCrgedGridEven,   1, MPI_INT,         0, MPI_COMM_WORLD); // icount1a  (local value)
	   
	   /*
	    * local function variables
	    */
	   MPI_Bcast(&noit,            1, MPI_INT,         0, MPI_COMM_WORLD); // noit      (local value)     
   
       if (0 == iNonIterateNum || fZero > fIonStrength || inhomo == 1)
       {
          mpi_itit(); 
       }     
       else
       {
          delphi_real qfact = abs(fNetCrg)*(delphi_real)iCrgBndyGridNum/iBndyGridNum;
          mpi_nitit(qfact);
       }
   }
}

#endif
