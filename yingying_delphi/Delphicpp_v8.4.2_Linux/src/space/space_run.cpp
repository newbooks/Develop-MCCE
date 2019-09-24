/*
 * Space_Run.cpp
 *
 */
#include <complex>
#include "space.h"
#include <boost/lexical_cast.hpp>

void CDelphiSpace::run()
{
    const char* infoString = " Info> ";

    delphi_integer i, j, k, ix, iy, iz, ic, ico, natom2;
    SGrid <delphi_integer> epstmp;
    SGrid <delphi_real> xl, xr;

    sgrid_temp_real.nX = 0.; sgrid_temp_real.nY = 0.; sgrid_temp_real.nZ = 0.;
    sgrid_temp_int.nX  = 0;   sgrid_temp_int.nY = 0;   sgrid_temp_int.nZ = 0;

    // initialize bDebMap
    get_pt3d <char>(idebmap, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1);
    for (ix = 1; ix <= iGrid.nX; ix++)
        for (iy = 1; iy <= iGrid.nY; iy++)
            for (iz = 1; iz <= iGrid.nZ; iz++)
                idebmap[ix][iy][iz] = true;

    // (ARGO) initialize zetaSurfMap
    get_pt3d <char>(zetaSurfMap, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1);
    for (ix = 1; ix <= iGrid.nX; ix++)
        for (iy = 1; iy <= iGrid.nY; iy++)
            for (iz = 1; iz <= iGrid.nZ; iz++)
                zetaSurfMap[ix][iy][iz] = true;

    if (iGaussian == 0)
        get_pt3d <SGrid <delphi_integer> >(iepsmp, global_iGrid.nX + 1, global_iGrid.nY + 1, global_iGrid.nZ + 1);
    else if (iGaussian == 1)
    {
        get_pt3d <SGrid <delphi_real> >(gepsmp, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1);
        get_pt3d <SGrid <delphi_real> >(gepsmp2, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1);
        get_pt3d <delphi_real>(gDensityMapOnGridPoint, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1);
    }

    if (iGaussian == 1 && inhomo == 0 && logs) //for 2nd Gaussian run
        for (i = 1; i <= iGrid.nX; i++)
            for (j = 1; j <= iGrid.nY; j++)
                for (k = 1; k <= iGrid.nZ; k++)
                    gepsmp[i][j][k] = fGepsMap_v[(i - 1)*iGrid.nY*iGrid.nZ + (j - 1)*iGrid.nY + k - 1];

    xn1 = &xn1_v[0]-1;
    xn2 = &xn2_v[0]-1;
    fRadPrb = &fRadPrb_v[0]-1;

    // for focusing
    if (ibctyp == 3)
    {
        SGrid <delphi_real> halfl;
        SGrid <delphi_real> edge_low, edge_high, xyz_temp;
        vector <CAtomPdb> delphipdb2; //delphipdb2 for focusing

        halfl     = (iGrid - delphi_integer(1)) / fScale / 2.0;
        edge_low  = acenter - halfl - 3.5;
        edge_high = acenter + halfl + 3.5;
        natom2    = 0;
		
		cout << "space_run> half1 = " << halfl << endl;
		cout << "space_run> acenter =  " << acenter.nX <<" " << acenter.nY << " " << acenter.nZ <<  endl;

        for (i = 0; i < iNatom; i++)
        {
            xyz_temp = delphipdb[i].getPose();
            if ( !(optORGT(xyz_temp, edge_high) || optORLT(xyz_temp, edge_low) ) )
            {
                delphipdb2.push_back(delphipdb[i]);
                natom2++;
            }
        }

        delphipdb.resize(natom2);

        for (i = 0; i < natom2; i++)
            delphipdb[i] = delphipdb2[i];

        iNatom = delphipdb.size();

        vector <CAtomPdb>().swap(delphipdb2);

        /*
         * write updated values back to global data container
         */
        pdc->getKey_Ref<delphi_integer>("natom")         = iNatom;
        pdc->getKey_Ref<vector <CAtomPdb> >("delphipdb") = delphipdb;
        myStart = {0,0,0};
        epsdim  = iNatom + iNObject + 2;
    }
    else
    {
        split();

        SGrid <delphi_real> xyz_temp;
        vector <CAtomPdb> delphipdb2; //temperary atom information storage

        natom2 = 0;

        for (i = 0; i < iNatom; i++)
        {
            xyz_temp = delphipdb[i].getPose();
            if ( !(optORGT(xyz_temp, myEndCoor) || optORLT(xyz_temp, myStartCoor) ) )
            {
                delphipdb2.push_back(delphipdb[i]);
                globalAtomIndex.push_back(i);
                natom2++;
            }
        }

        // update the atoms handled by current process
        delphipdb.swap(delphipdb2);
        iNatom = natom2;

        SExtrema<delphi_real> tmpExtrema = { {0,0,0},{0,0,0} };
        if (iNatom > 0)
            tmpExtrema = { delphipdb[0].getPose(), delphipdb[0].getPose() };

        for (int i = 0; i < iNatom; i++)
        {
            tmpExtrema.nMin = optMin<delphi_real>(tmpExtrema.nMin, delphipdb[i].getPose());
            tmpExtrema.nMax = optMax<delphi_real>(tmpExtrema.nMax, delphipdb[i].getPose());
        }
        sLimObject[0]=tmpExtrema;

    }

    sDelPhiPDB = new delphipdb_struc[iNatom + 1];

    for (int i = 0; i <= iNatom - 1; i++)
    {
        sDelPhiPDB[i + 1].radius = delphipdb[i].getRadius();
        sDelPhiPDB[i + 1].xyz = delphipdb[i].getPose();
        sDelPhiPDB[i + 1].charge = delphipdb[i].getCharge();
        sDelPhiPDB[i + 1].atom_resname = delphipdb[i].getAtResname(); //Argo
    }

    for (int i = 1; i <= iNatom; i++)
    {
        xn1[i] = sDelPhiPDB[i].xyz;
        if (ibctyp != 3)
            xn2[i] = (sDelPhiPDB[i].xyz - global_cOldMid)*fScale + fRMid - optCast<delphi_real,delphi_integer>(myStart);
        else
            xn2[i] = (sDelPhiPDB[i].xyz - global_cOldMid)*fScale + (iGrid+1)/2.0;
    }

    epsmak(); //Lin Li reset

    // Now start crgarr
    if (isolv)
    {
        // increased the starting dimension of crgatn and..
        extracrg = 0;
        if (ndistr > 0) extracrg = iGrid.nX*iGrid.nY*iGrid.nZ;

        /*
         * 2011-05-30 Allocation of the arrays below is moved to the body of crgarr void,
         * arrays are accessible via pointers module. Sizes of arrays are determined before
         * allocation inside the crgarr void
         */
        crgarr(); //Lin Li reset

        xl = myStartCoor;
        xr = myEndCoor;

        if (logs || lognl)
        {
            ico = 0;
            for (ic = 1; ic <= nqass; ic++)
            {
                if (optORLT(chgpos[ic], xl) || optORGT(chgpos[ic], xr))
                {
                    #ifdef VERBOSE
                    if (crgatn[ic] < 0)
                        cout << "!WARNING: distribution " << -crgatn[ic] << "outside the box" << endl;
                    else
                    {
                        if (crgatn[ic] > iNatom)
                            cout << "WARNING:crg " << ic << "object " << crgatn[ic] - iNatom << "outside the box " << chgpos[ic] << endl;
                        else
                            cout << "//!! WARNING : charge " << delphipdb[crgatn[ic] - 1].getAtInf() << "outside the box" << endl;
                    }
                    #endif

                    ico = 1;
                }
            }

            if (ico > 0 && ibctyp != 3)
            {
                cout << "CHARGES OUTSIDE THE BOX AND NOT DOING FOCUSSING << THEREFORE STOP" << endl;
                exit(0);
            }
        }
    }

    // initialize iepsmp
    bDebMap_v.assign(iGrid.nX*iGrid.nY*iGrid.nZ, true);
    zetaSurfMap_v.assign(iGrid.nX*iGrid.nY*iGrid.nZ, true);
    if (iGaussian == 0)

        iEpsMap_v.assign(iGrid.nX*iGrid.nY*iGrid.nZ, { repsout,repsout ,repsout });
    else if (iGaussian == 1)
    {
        fGepsMap_v.assign(iGrid.nX*iGrid.nY*iGrid.nZ, sgrid_temp_real);
        fGepsMap2_v.assign(iGrid.nX*iGrid.nY*iGrid.nZ, sgrid_temp_real);
    }

    // -------gepsmp file: --------
//    {
//        ofstream densfile;
//        densfile.open ("test_space_gaussian_gepsmp0.txt");
//
//        densfile << fixed << setprecision(7);
//        for(ix=1; ix<=iGrid.nX; ix++)
//            for(iy=1; iy<=iGrid.nY; iy++)
//                for(iz=1; iz<=iGrid.nZ; iz++)
//                    densfile << setw(6) << right << ix << " "
//                             << setw(6) << right << iy << " "
//                             << setw(6) << right << iz << " "
//                             << setw(8) << right << gepsmp[ix][iy][iz].nX << " "
//                             << setw(8) << right << gepsmp[ix][iy][iz].nY << " "
//                             << setw(8) << right << gepsmp[ix][iy][iz].nZ << endl;
//        densfile.close();
//    }

    if (iGaussian == 0)
        for (k = 1; k <= iGrid.nZ; k++)
            for (j = 1; j <= iGrid.nY; j++)
                for (i = 1; i <= iGrid.nX; i++)
                    iEpsMap_v[(k - 1)*iGrid.nY*iGrid.nX + (j - 1)*iGrid.nX + i - 1] = iepsmp[k][j][i];
    else if (iGaussian == 1)
        for (k = 1; k <= iGrid.nZ; k++)
            for (j = 1; j <= iGrid.nY; j++)
                for (i = 1; i <= iGrid.nX; i++)
                {
                    //fGepsMap2_v[(k - 1)*iGrid.nY*iGrid.nX + (j - 1)*iGrid.nX + i - 1] = gepsmp2[i][j][k];
                    //fGepsMap_v[(k - 1)*iGrid.nY*iGrid.nX + (j - 1)*iGrid.nX + i - 1] = gepsmp[i][j][k];
                    fGepsMap2_v[(i-1)*iGrid.nY*iGrid.nZ+(j-1)*iGrid.nZ+k-1]=gepsmp2[i][j][k];
                    fGepsMap_v[(i-1)*iGrid.nY*iGrid.nZ+(j-1)*iGrid.nZ+k-1]=gepsmp[i][j][k];
                }

    // Gaussian Density Map
    if ( !(iGaussian == 1 && inhomo == 0 && logs) )
    {
        fGDensityMap_v.assign(iGrid.nX*iGrid.nY*iGrid.nZ, 0.0);

        if (fIonStrenth > fZero)
        {
            if (iGaussian != 0)
                for (k = 1; k <= iGrid.nZ; k++)
                    for (j = 1; j <= iGrid.nY; j++)
                        for (i = 1; i <= iGrid.nX; i++)
                            fGDensityMap_v[(k - 1)*iGrid.nY*iGrid.nX + (j - 1)*iGrid.nX + i - 1] = gDensityMapOnGridPoint[i][j][k];
            else
                for (k = 1; k <= iGrid.nZ; k++)
                    for (j = 1; j <= iGrid.nY; j++)
                        for (i = 1; i <= iGrid.nX; i++)
                            if (!idebmap[i][j][k]) fGDensityMap_v[(k - 1)*iGrid.nY*iGrid.nX + (j - 1)*iGrid.nX + i - 1] = 1.0;

        }
    }

    for (k = 1; k <= iGrid.nZ; k++)
    {
        for (j = 1; j <= iGrid.nY; j++)
        {
            for (i = 1; i <= iGrid.nX; i++)
            {
                //ARGO commented the ifdef in original. Place ifdef-KJI later
                #ifdef IKJ
                bDebMap_v[(i - 1)*iGrid*iGrid + (j - 1)*iGrid + (k - 1)] = idebmap[k][j][i];
                if (zetaOn == 1) zetaSurfMap_v[(i - 1)*iGrid*iGrid + (j - 1)*iGrid + (k - 1)] = zetaSurfMap[k][j][i]; //ARGO
                #endif // IKJ

                bDebMap_v[(k - 1)*iGrid.nY*iGrid.nX + (j - 1)*iGrid.nX + i - 1] = idebmap[k][j][i];
                if (zetaOn == 1) zetaSurfMap_v[(k - 1)*iGrid.nY*iGrid.nX + (j - 1)*iGrid.nX + i - 1] = zetaSurfMap[i][j][k]; //ARGO
            }
        }
    }

    if (iepsmp != NULL) free_pt3d(iepsmp, global_iGrid.nX + 1, global_iGrid.nY + 1, global_iGrid.nZ + 1);
    if (idebmap != NULL) free_pt3d(idebmap, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1);
    if (gepsmp != NULL) free_pt3d(gepsmp, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1);
    if (gepsmp2 != NULL) free_pt3d(gepsmp2, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1);
    if (gDensityMapOnGridPoint != NULL) free_pt3d(gDensityMapOnGridPoint, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1);
    if (zetaSurfMap != NULL) free_pt3d(zetaSurfMap, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1); //ARGO
}

