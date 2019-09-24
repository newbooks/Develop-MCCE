/**
 * In:
 *    delphi_integer ibnum = iBndyGridNum;
 *    vector< SGrid<delphi_integer> > iepsmp = prgigEpsMap;
 *    vector<bool> idebmap = prgbDielecMap;
 *    vector< SGrid<delphi_integer> > ibgrd = prgigBndyGrid;
 *
 * return:
 *    delphi_integer                iDielecBndyEven     icount2a   used for realigning idpos and db,compressing to contingous space
 *    delphi_integer                iDielecBndyOdd      icount2b   used for realigning idpos and db,compressing to contingous space
 *    vector<delphi_integer>        prgiBndyDielecIndex idpos(nsp) indices of the dielectric boundary points used to recover the boundary values
 *    vector< vector<delphi_real> > prgfBndyDielec      db(6,nsp)  dielectric values on the boundary
 *    vector<delphi_real>           prgfSaltMap1        sf1(nhgp)  coefficients used in GS/SOR iterations
 *    vector<delphi_real>           prgfSaltMap2        sf2(nhgp)  coefficients used in GS/SOR iterations
 */

/*
 * sample code for creating multi-dimensional vector:
 * int num_of_col = 5;
 * int num_of_row = 9;
 * double init_value = 3.14;
 *
 * vector< vector<double> > matrix;
 * //now we have an empty 2D-matrix of size (0,0). Resizing it with one single command:
 * matrix.resize( num_of col , vector<double>( num_of_row , init_value ) );
 * // and we are good to go ...
 */

#include "solver_fastSOR.h"

//-----------------------------------------------------------------------//
void CDelphiFastSOR::setDielecBndySaltMap() 
{
    //++++++++++ INPUT:
    vector<SGrid<delphi_integer> >::const_iterator iepsmp = prgigEpsMap.cbegin();
    vector<char>::const_iterator idebmap = prgbDielecMap.cbegin();
    //vector< SGrid<delphi_integer> >::const_iterator ibgrd = prgigBndyGrid.cbegin();

    //++++++++++ LOCAL:
    vector<int> it(6, 0);                    // it(6) rgiDBIndex(6,0)
    vector < delphi_real > vecttemp(6, 0.0); // vecttemp rgfDBTemp(6,0.0)
    delphi_integer ieps;                     // iEps
    delphi_real temp;                        // fTempVal
    int deb;                                 // iDielec
    delphi_integer idbs = 0;                 // iDielecBndyGridNum
    vector < delphi_integer > iepsv;         // iepsv(nsp) prgiEpsv

    //++++++++++ OUTPUT:
    delphi_integer& icount2a = iDielecBndyEven;
    delphi_integer& icount2b = iDielecBndyOdd;
    vector < delphi_integer > &idpos = prgiBndyDielecIndex;
    vector < vector<delphi_real> > &db = prgfBndyDielec;
    vector < delphi_real > &sf1 = prgfSaltMap1;
    vector < delphi_real > &sf2 = prgfSaltMap2;

    //-------------------------------- dbsfd -----------------------------//
    /*
     *  create rgfDielecBndyValue(dbval) and rgfSaltFuncDiff(sfd) to store the values
     *  later to be assigned later in db(db) and prgSaltMap(sf1-2)
     */
    delphi_real rgfDielecBndyValue[2][7][2];       // dbval(0:1,0:6,0:1)
    delphi_real rgfSaltFuncDiff[6][2];             // sfd(5,0:1)
    delphi_real fDenom;                            // denom

    //fGepsMap2_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1]=gepsmp2[i][j][k];

    if (debug_solver) 
    {
        cout << "##### starting dbsf:  ###" << endl;
        cout << "icount2a: " << icount2a << endl;
        cout << "icount2b: " << icount2a << endl;
        pTimer->showTime();
    }
    
    if (debug_solver)
        cout << "iGaussian " << iGaussian << endl;
    if (iGaussian == 1)
        bDbOut = true; //for Gaussian

    if (0.0 < fIonStrength) 
    {
        for (int iz = 0; iz <= 1; iz += 1) 
        {
            for (int iy = 1; iy <= 3; iy += 1) 
            {
                fDenom = fSixEps + iy * fEpsDiff + iz * fDebFct;
                rgfDielecBndyValue[0][iy][iz] = 0.0;
                rgfDielecBndyValue[1][iy][iz] = fEpsDiff / fDenom;
                rgfSaltFuncDiff[iy][iz] = fEpsOut / fDenom;
            }
        }

        for (int iz = 0; iz <= 1; iz += 1) 
        {
            for (int iy = 4; iy <= 5; iy += 1) 
            {
                fDenom = fSixEps + iy * fEpsDiff + iz * fDebFct;
                rgfDielecBndyValue[0][iy][iz] = -fEpsDiff / fDenom;
                rgfDielecBndyValue[1][iy][iz] = 0.0;
                rgfSaltFuncDiff[iy][iz] = fEpsOut / fDenom;
            }
        }
    } 
    else 
    {
        for (int iz = 0; iz <= 1; iz += 1) 
        {
            for (int iy = 1; iy <= 5; iy += 1) 
            {
                fDenom = fSixEps + iy * fEpsDiff;
                rgfDielecBndyValue[0][iy][iz] = fEpsOut / fDenom - fSixth;
                rgfDielecBndyValue[1][iy][iz] = fEpsIn / fDenom - fSixth;
            }
        }
    }
    //---------------------------- End of dbsfd.f -----------------------//

    if (debug_solver) 
    {
        cout << "#####  dbsf:  line 120 ###" << endl;
        pTimer->showTime();
    }

    //-------------------------------- mkdbsf ----------------------------//
    int ix, iy, iz;
    delphi_integer iw;

    string strDbFile = "db.dat"; // db file name
    ofstream ofDbFileStream;

    if (bDbOut) //in fortran: idbwrt
    {
        ofDbFileStream.open(strDbFile.c_str());
        ofDbFileStream << fixed << setprecision(3);
        ofDbFileStream << "DELPHI DB FILE" << endl;
        ofDbFileStream << "FORMAT NUMBER=1" << endl;
        ofDbFileStream << "NUMBER OF BOUNDARY POINTS= " << iBndyGridNum << endl;
    }

    /*vector<delphi_integer>::iterator idposEven = idpos.end();
     vector< vector<delphi_real> >::iterator dbEven = db.end();
     vector<delphi_real>::iterator densityEven = gaussianBoundaryDensity.begin();
     vector<vector<delphi_real>>::iterator gdbEven = gaussianBoundaryDielec.begin();
     vector<delphi_integer>::iterator iepsvEven = iepsv.end();*/

    //predetermin number of even/odd points 
    for (vector<SGrid<delphi_integer> >::const_iterator const_itr = prgigBndyGrid.cbegin(); const_itr != prgigBndyGrid.cend(); ++const_itr) 
    {
        ix = const_itr->nX;
        iy = const_itr->nY;
        iz = const_itr->nZ;
        iw = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);

        if (0 == iw % 2) // iw is even
        {
            icount2a += 1;
        }

        icount2b += 1;
    }

    if (debug_solver) 
    {
        cout << "#####  dbsf:  line 178 ### icount2a=" << icount2a << "  icount2b=" << icount2b << endl;
        pTimer->showTime();
    }

    //preset array size
    idpos.resize(icount2b);
    db.resize(icount2b);
    gaussianBoundaryDensity.resize(icount2b);
    gaussianBoundaryDielec.resize(icount2b);
    iepsv.resize(icount2b);

    //------ Due to the nature of vector in C++, realignment of idpos, db, sf1 and sf2 is unnecessary.
    //preset array Offset for Even and Odd
    vector<delphi_integer>::iterator idposEven    = idpos.begin();
    vector<vector<delphi_real> >::iterator dbEven = db.begin();
    vector<delphi_real>::iterator densityEven     = gaussianBoundaryDensity.begin();
    vector<vector<delphi_real>>::iterator gdbEven = gaussianBoundaryDielec.begin();
    vector<delphi_integer>::iterator iepsvEven    = iepsv.begin();

    vector<delphi_integer>::iterator idposOdd     = idpos.begin() + icount2a;
    vector<vector<delphi_real> >::iterator dbOdd  = db.begin() + icount2a;
    vector<delphi_real>::iterator densityOdd      = gaussianBoundaryDensity.begin() + icount2a;
    vector<vector<delphi_real>>::iterator gdbOdd  = gaussianBoundaryDielec.begin() + icount2a;
    //vector<delphi_integer>::iterator iepsvOdd     = iepsv.end() + icount2a; //This line causes errors

    delphi_integer icountEven = 0;
    delphi_integer icountOdd  = 0;

    //Now loop over all boundary points
    for (vector<SGrid<delphi_integer> >::const_iterator const_itr = prgigBndyGrid.cbegin(); const_itr != prgigBndyGrid.cend(); ++const_itr) 
    {
        ix = const_itr->nX;
        iy = const_itr->nY;
        iz = const_itr->nZ;

        if (0 == iDirectEpsMap) 
        {
            it.assign(6, 0);

            iw = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);
            if (0 != (iepsmp + iw)->nX / iEpsDim) it[0] = 1;
            if (0 != (iepsmp + iw)->nY / iEpsDim) it[1] = 1;
            if (0 != (iepsmp + iw)->nZ / iEpsDim) it[2] = 1;

            iw = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1 - 1);
            if (0 != (iepsmp + iw)->nX / iEpsDim) it[3] = 1;

            iw = (iz - 1) * iGrid * iGrid + (iy - 1 - 1) * iGrid + (ix - 1);
            if (0 != (iepsmp + iw)->nY / iEpsDim) it[4] = 1;

            iw = (iz - 1 - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);
            if (0 != (iepsmp + iw)->nZ / iEpsDim) it[5] = 1;

            ieps = it[0] + it[1] + it[2] + it[3] + it[4] + it[5];
        } 
        else //Here is the key point for Gaussian
        {
            ieps = 0;
            temp = 0.0;

            delphi_integer iIndex;

            if (iGaussian == 0)
            {
                iw = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);
                iIndex      = (iepsmp + iw)->nX / iEpsDim;
                temp       += prgfMediaEps[iIndex];
                vecttemp[0] = prgfMediaEps[iIndex];

                iIndex      = (iepsmp + iw)->nY / iEpsDim;
                temp       += prgfMediaEps[iIndex];
                vecttemp[1] = prgfMediaEps[iIndex];

                iIndex      = (iepsmp + iw)->nZ / iEpsDim;
                temp       += prgfMediaEps[iIndex];
                vecttemp[2] = prgfMediaEps[iIndex];

                iw     = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1 - 1);
                iIndex = (iepsmp + iw)->nX / iEpsDim;
                temp  += prgfMediaEps[iIndex];
                vecttemp[3] = prgfMediaEps[iIndex];

                iw     = (iz - 1) * iGrid * iGrid + (iy - 1 - 1) * iGrid + (ix - 1);
                iIndex = (iepsmp + iw)->nY / iEpsDim;
                temp  += prgfMediaEps[iIndex];
                vecttemp[4] = prgfMediaEps[iIndex];

                iw     = (iz - 1 - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);
                iIndex = (iepsmp + iw)->nZ / iEpsDim;
                temp  += prgfMediaEps[iIndex];
                vecttemp[5] = prgfMediaEps[iIndex];
            }
            else if (iGaussian == 1) //Gaussian
            {
                iw    = (ix - 1) * iGrid * iGrid + (iy - 1) * iGrid + (iz - 1);
                temp += gepsmp2[iw].nX;
                vecttemp[0] = gepsmp2[iw].nX;

                temp += gepsmp2[iw].nY;
                vecttemp[1] = gepsmp2[iw].nY;

                temp += gepsmp2[iw].nZ;
                vecttemp[2] = gepsmp2[iw].nZ;

                //iw = (ix-1-1)*iGrid*iGrid + (iy-1)*iGrid + (iz-1);
                //iw = iw - iGrid*iGrid;
                delphi_integer iw3 = iw - iGrid * iGrid;
                temp += gepsmp2[iw3].nX;
                vecttemp[3] = gepsmp2[iw3].nX;

                //iw = (ix-1)*iGrid*iGrid + (iy-1-1)*iGrid + (iz-1);
                //iw = iw + iGrid*(iGrid - 1);
                delphi_integer iw4 = iw - iGrid;
                temp += gepsmp2[iw4].nY;
                vecttemp[4] = gepsmp2[iw4].nY;

                //iw = (ix-1)*iGrid*iGrid + (iy-1)*iGrid + (iz-1-1);
                //iw = iw + iGrid - 1;
                delphi_integer iw5 = iw - 1;
                temp += gepsmp2[iw5].nZ;
                vecttemp[5] = gepsmp2[iw5].nZ;
            }
        }

        deb = 0;

        //iw = (iz-1-1)*iGrid*iGrid + (iy-1-1)*iGrid + (ix-1);
        iw = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);
        //iw = (ix-1)*iGrid*iGrid + (iy-1)*iGrid + (iz-1);

        if (*(idebmap + iw)) 
        {
            deb = 1;
            idbs += 1;
        }

        delphi_real gridDensity;
        if (iGaussian==1)
        {
            gridDensity = gaussianDensityMap[iw];
        }

        //if (gridDensity < fvdwdens) gridDensity = 0;

        vector < delphi_real > dbrow;

        if (0 == iDirectEpsMap) 
        {
            dbrow.push_back(rgfDielecBndyValue[it[3]][ieps][deb]);
            dbrow.push_back(rgfDielecBndyValue[it[0]][ieps][deb]);
            dbrow.push_back(rgfDielecBndyValue[it[4]][ieps][deb]);
            dbrow.push_back(rgfDielecBndyValue[it[1]][ieps][deb]);
            dbrow.push_back(rgfDielecBndyValue[it[5]][ieps][deb]);
            dbrow.push_back(rgfDielecBndyValue[it[2]][ieps][deb]);
        } 
        else 
        {
            fDenom = temp + deb * fDebFct;

            if (0.0 == fIonStrength) 
            {
                dbrow.push_back(vecttemp[3] / fDenom - fSixth);
                dbrow.push_back(vecttemp[0] / fDenom - fSixth);
                dbrow.push_back(vecttemp[4] / fDenom - fSixth);
                dbrow.push_back(vecttemp[1] / fDenom - fSixth);
                dbrow.push_back(vecttemp[5] / fDenom - fSixth);
                dbrow.push_back(vecttemp[2] / fDenom - fSixth);
            } 
            else 
            {
                dbrow.push_back(vecttemp[3] / fDenom);
                dbrow.push_back(vecttemp[0] / fDenom);
                dbrow.push_back(vecttemp[4] / fDenom);
                dbrow.push_back(vecttemp[1] / fDenom);
                dbrow.push_back(vecttemp[5] / fDenom);
                dbrow.push_back(vecttemp[2] / fDenom);
            }
        }

        vector < delphi_real > dbrow_original;
        dbrow_original.push_back(vecttemp[3]);
        dbrow_original.push_back(vecttemp[0]);
        dbrow_original.push_back(vecttemp[4]);
        dbrow_original.push_back(vecttemp[1]);
        dbrow_original.push_back(vecttemp[5]);
        dbrow_original.push_back(vecttemp[2]);

        iw = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1); // iw=isgrid*(k-1) + igrid*(j-1) + i

        /*if (0 == iw%2) // iw is even
         {
			 idposEven = idpos.begin() + icount2a;
			 idposEven = idpos.insert(idposEven,iw/2+1);
	
			 iepsvEven = iepsv.begin() + icount2a;
			 iepsvEven = iepsv.insert(iepsvEven,ieps);
			 iepsvEven++;
	
			 dbEven    = db.begin() + icount2a;
			 dbEven    = db.insert(dbEven,dbrow);
			 dbEven++;
	
			 densityEven = gaussianBoundaryDensity.begin() + icount2a;
			 densityEven = gaussianBoundaryDensity.insert(densityEven, gridDensity);
	
			 gdbEven = gaussianBoundaryDielec.begin() + icount2a;
			 gdbEven = gaussianBoundaryDielec.insert(gdbEven, dbrow_original);
	
			 icount2a += 1;
         }
         else // iw is odd
         {
			 idpos.push_back( (iw+1)/2 );
			 iepsv.push_back(ieps);
			 db.push_back(dbrow);
			 gaussianBoundaryDensity.push_back(gridDensity);
			 gaussianBoundaryDielec.push_back(dbrow_original);
         }
         icount2b++;
         */

        //Faster method than the commented lines above
        if (0 == iw % 2) // iw is even 
                {
            /**idposEven = iw / 2 + 1;
             idposEven++;

             *iepsvEven = ieps;
             iepsvEven++;

             *dbEven = dbrow;
             dbEven++;

             *densityEven = gridDensity;
             densityEven++;

             *gdbEven = dbrow_original;
             gdbEven++;*/

            *(idpos.begin() + icountEven) = iw / 2 + 1;
            *(iepsv.begin() + icountEven) = ieps;
            *(db.begin() + icountEven) = dbrow;
            *(gaussianBoundaryDensity.begin() + icountEven) = gridDensity;
            *(gaussianBoundaryDielec.begin() + icountEven) = dbrow_original;
            icountEven++;
        } 
        else // iw is odd
        {
            /**idposOdd = (iw + 1) / 2;
             idposOdd++;

             *iepsvOdd = ieps;
             iepsvOdd++;

             *dbOdd = dbrow;
             dbOdd++;

             *densityOdd = gridDensity;
             densityOdd++;

             *gdbOdd = dbrow_original;
             gdbOdd++;*/

            *(idpos.begin() + icount2a + icountOdd) = (iw + 1) / 2;
            *(iepsv.begin() + icount2a + icountOdd) = ieps;
            *(db.begin() + icount2a + icountOdd) = dbrow;
            *(gaussianBoundaryDensity.begin() + icount2a + icountOdd) = gridDensity;
            *(gaussianBoundaryDielec.begin() + icount2a + icountOdd) = dbrow_original;
            icountOdd++;
        }

        if (bDbOut) 
        {
            ofDbFileStream << setw(3) << left << ix << " " << setw(3) << left << iy << " " << setw(3) << left << iz << " " << setw(8) << left << dbrow[0] << setw(8) << left << dbrow[1]
            << setw(8) << left << dbrow[2] << setw(8) << left << dbrow[3] << setw(8) << left << dbrow[4] << setw(8) << left << dbrow[5] << endl;
        }

        // dbs=dbs+db(1,ibnum3) is NOT used
    } //---------- end of for (delphi_integer ix = 0; ix < iBndyGridNum; ix += 1)

    if (debug_solver) 
    {
        cout << "#####  dbsf:  line 359 ###" << endl;
        pTimer->showTime();
    }

    if (bDbOut) ofDbFileStream.close();

#ifdef VERBOSE
    cout << " Info> Number of dielectric boundary points in salt = " << idbs << endl;
#endif

    //---------- realign idpos and db,compressing to contingous space

    //---------- set saltmaps 1 and 2, i.e., sf1 and sf2.
    if (0.0 < fIonStrength) 
    {
        const delphi_real fSixSalt = fSixth * (1.0 / (1.0 + fDebFct / fSixEps) - 1.0); // sixsalt

        iw = 0;

        for (vector<char>::const_iterator const_itr = prgbDielecMap.cbegin(); const_itr != prgbDielecMap.cend(); ++const_itr) 
        {
            if (0 == iw % 2) // even pts
            {
                if (*const_itr)
                    sf1.push_back(fSixth + fSixSalt);
                else
                    sf1.push_back(fSixth);
            } 
            else // odd pts
            {
                if (*const_itr)
                    sf2.push_back(fSixth + fSixSalt);
                else
                    sf2.push_back(fSixth);
            }

            iw += 1;
        }

        sf2.push_back(0.0); //sf[iGrid3]

#ifdef DEBUG_DELPHI_SOLVER_MKDBSF1
        {
            string strTestFile = "test_mkdbsf1.dat";
            ofstream ofTestStream(strTestFile.c_str());
            ofTestStream << boolalpha;
            ofTestStream << fixed << setprecision(7);

            int count = 1;
            for (vector<delphi_real>::const_iterator const_itr = sf1.cbegin(); const_itr != sf1.cend(); ++const_itr)
            {
                ofTestStream << "sf1(" << setw(6) << right << count << ") = " << setw(11) << right << *const_itr << endl;
                count += 1;
            }

            count = 1;
            for (vector<delphi_real>::const_iterator const_itr = sf2.cbegin(); const_itr != sf2.cend(); ++const_itr)
            {
                ofTestStream << "sf2(" << setw(6) << right << count << ") = " << setw(11) << right << *const_itr << endl;
                count += 1;
            }

            ofTestStream.close();
        }
#endif // DEBUG_DELPHI_SOLVER_MKDBSF1

        for (int i = 0; i < icount2a; i++) 
        {
            if (0 != iDirectEpsMap)
                sf1[idpos[i] - 1] = 0.0;
            else 
            {
                if (abs(sf1[idpos[i] - 1] - fSixth) < fZero)
                    sf1[idpos[i] - 1] = rgfSaltFuncDiff[iepsv[idpos[i] - 1]][0];
                else
                    sf1[idpos[i] - 1] = rgfSaltFuncDiff[iepsv[idpos[i] - 1]][1];
            }
        }

        for (int i = icount2a; i < icount2b; i++) 
        {
            if (0 != iDirectEpsMap)
                sf2[idpos[i] - 1] = 0.0;
            else 
            {
                if (abs(sf2[idpos[i] - 1] - fSixth) < fZero)
                    sf2[idpos[i] - 1] = rgfSaltFuncDiff[iepsv[idpos[i] - 1]][0];
                else
                    sf2[idpos[i] - 1] = rgfSaltFuncDiff[iepsv[idpos[i] - 1]][1];
            }
        }
    }

    if (debug_solver) 
    {
        cout << "#####  dbsf:  line 448 ###" << endl;
        pTimer->showTime();
    }

    //---------------------------- end of mkdbsf.f -----------------------//

#ifdef DEBUG_DELPHI_SOLVER_MKDBSF
    {
        string strTestFile = "test_mkdbsf_"+to_string(inhomo)+".dat";
        ofstream ofTestStream(strTestFile.c_str());
        ofTestStream << boolalpha;
        ofTestStream << fixed << setprecision(7);

        ofTestStream << "icount2a = " << setw(6) << right << icount2a << " icount2b = " << setw(6) << right << icount2b << endl;

        int count = 1;
        for (vector<delphi_integer>::const_iterator const_itr = idpos.cbegin(); const_itr != idpos.cend(); ++const_itr)
        {
            ofTestStream << "idpos(" << setw(6) << right << count << ") = " << setw(6) << right << *const_itr << endl;
            count += 1;
        }

        count = 1;
        for (vector< vector<delphi_real> >::const_iterator const_itr = db.cbegin(); const_itr != db.cend(); ++const_itr)
        {
            ofTestStream << "db(" << setw(6) << right << count << ") = " << setw(11) << right << const_itr->at(0) << setw(11) << right << const_itr->at(1)
            << setw(11) << right << const_itr->at(2) << setw(11) << right << const_itr->at(3) << setw(11) << right << const_itr->at(4) << setw(11) << right << const_itr->at(5) << endl;
            count += 1;
        }

        count = 1;
        for (vector<delphi_real>::const_iterator const_itr = sf1.cbegin(); const_itr != sf1.cend(); ++const_itr)
        {
            ofTestStream << "sf1(" << setw(6) << right << count << ") = " << setw(11) << right << *const_itr << endl;
            count += 1;
        }

        count = 1;
        for (vector<delphi_real>::const_iterator const_itr = sf2.cbegin(); const_itr != sf2.cend(); ++const_itr)
        {
            ofTestStream << "sf2(" << setw(6) << right << count << ") = " << setw(11) << right << *const_itr << endl;
            count += 1;
        }

        ofTestStream.close();
    }
#endif // DEBUG_DELPHI_SOLVER_MKDBSF

    if (debug_solver) 
    {
        cout << "### out of dbsf... ### " << endl;
        pTimer->showTime();
    }

}

#ifdef PARALLEL_MPI

//-----------------------------------------------------------------------//
void CDelphiFastSOR::mpi_setDielecBndySaltMap()
{
    //++++++++++ INPUT:
    //vector<SGrid<delphi_integer> >::const_iterator iepsmp = prgigEpsMap.cbegin();
    SGrid<delphi_integer> iepsmp_val;

    vector<char>::const_iterator idebmap = prgbDielecMap.cbegin();
    //vector< SGrid<delphi_integer> >::const_iterator ibgrd = prgigBndyGrid.cbegin();

    //++++++++++ LOCAL:
    vector<int> it(6, 0);                    // it(6) rgiDBIndex(6,0)
    vector < delphi_real > vecttemp(6, 0.0); // vecttemp rgfDBTemp(6,0.0)
    delphi_integer ieps;                     // iEps
    delphi_real temp;                        // fTempVal
    int deb;                                 // iDielec
    delphi_integer idbs = 0;                 // iDielecBndyGridNum
    vector < delphi_integer > iepsv;         // iepsv(nsp) prgiEpsv

    //++++++++++ OUTPUT:
    delphi_integer& icount2a = iDielecBndyEven;
    delphi_integer& icount2b = iDielecBndyOdd;
    vector < delphi_integer > &idpos = prgiBndyDielecIndex;
    vector < vector<delphi_real> > &db = prgfBndyDielec;
    vector < delphi_real > &sf1 = prgfSaltMap1;
    vector < delphi_real > &sf2 = prgfSaltMap2;

    //-------------------------------- dbsfd -----------------------------//
    /*
     *  create rgfDielecBndyValue(dbval) and rgfSaltFuncDiff(sfd) to store the values
     *  later to be assigned later in db(db) and prgSaltMap(sf1-2)
     */
    delphi_real rgfDielecBndyValue[2][7][2];       // dbval(0:1,0:6,0:1)
    delphi_real rgfSaltFuncDiff[6][2];             // sfd(5,0:1)
    delphi_real fDenom;                            // denom

    if (iGaussian == 1)
        bDbOut = true; //for Gaussian

    if (0.0 < fIonStrength)
    {
        for (int iz = 0; iz <= 1; iz += 1)
        {
            for (int iy = 1; iy <= 3; iy += 1)
            {
                fDenom = fSixEps + iy * fEpsDiff + iz * fDebFct;
                rgfDielecBndyValue[0][iy][iz] = 0.0;
                rgfDielecBndyValue[1][iy][iz] = fEpsDiff / fDenom;
                rgfSaltFuncDiff[iy][iz] = fEpsOut / fDenom;
            }
        }

        for (int iz = 0; iz <= 1; iz += 1)
        {
            for (int iy = 4; iy <= 5; iy += 1)
            {
                fDenom = fSixEps + iy * fEpsDiff + iz * fDebFct;
                rgfDielecBndyValue[0][iy][iz] = -fEpsDiff / fDenom;
                rgfDielecBndyValue[1][iy][iz] = 0.0;
                rgfSaltFuncDiff[iy][iz] = fEpsOut / fDenom;
            }
        }
    }
    else
    {
        for (int iz = 0; iz <= 1; iz += 1)
        {
            for (int iy = 1; iy <= 5; iy += 1)
            {
                fDenom = fSixEps + iy * fEpsDiff;
                rgfDielecBndyValue[0][iy][iz] = fEpsOut / fDenom - fSixth;
                rgfDielecBndyValue[1][iy][iz] = fEpsIn / fDenom - fSixth;
            }
        }
    }
    //---------------------------- End of dbsfd.f -----------------------//

    //-------------------------------- mkdbsf ----------------------------//
    int ix, iy, iz;
    delphi_integer iw;

    string strDbFile = "db.dat"; // db file name
    ofstream ofDbFileStream;

    if (bDbOut) //in fortran: idbwrt
    {
        ofDbFileStream.open(strDbFile.c_str());
        ofDbFileStream << fixed << setprecision(3);
        ofDbFileStream << "DELPHI DB FILE" << endl;
        ofDbFileStream << "FORMAT NUMBER=1" << endl;
        ofDbFileStream << "NUMBER OF BOUNDARY POINTS= " << iBndyGridNum << endl;
    }

    /*
     * predetermin number of even/odd points
     */
    size_t mpi_size = pdc->sizeofGlobal1D<SGrid<delphi_integer>>("ibgrd");
    SGrid<delphi_integer> mpi_ibgrd;

    for (size_t i = 0; i < mpi_size; i++)
    {
        mpi_ibgrd = pdc->readGlobalVector1D< SGrid<delphi_integer> >("ibgrd", i);

        ix = mpi_ibgrd.nX;
        iy = mpi_ibgrd.nY;
        iz = mpi_ibgrd.nZ;
        iw = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);

        if (0 == iw % 2) // iw is even
        {
            icount2a += 1;
        }

        icount2b += 1;
    }

    /*
     * preset array size
     */
    idpos.resize(icount2b);
    db.resize(icount2b);
    gaussianBoundaryDensity.resize(icount2b);
    gaussianBoundaryDielec.resize(icount2b);
    iepsv.resize(icount2b);

    /*
     * Due to the nature of vector in C++, realignment of idpos, db, sf1 and sf2 is unnecessary.
     *
     * preset array Offset for Even and Odd
     */
    vector<delphi_integer>::iterator idposEven    = idpos.begin();
    vector<vector<delphi_real> >::iterator dbEven = db.begin();
    vector<delphi_real>::iterator densityEven     = gaussianBoundaryDensity.begin();
    vector<vector<delphi_real>>::iterator gdbEven = gaussianBoundaryDielec.begin();
    vector<delphi_integer>::iterator iepsvEven    = iepsv.begin();

    vector<delphi_integer>::iterator idposOdd     = idpos.begin() + icount2a;
    vector<vector<delphi_real> >::iterator dbOdd  = db.begin() + icount2a;
    vector<delphi_real>::iterator densityOdd      = gaussianBoundaryDensity.begin() + icount2a;
    vector<vector<delphi_real>>::iterator gdbOdd  = gaussianBoundaryDielec.begin() + icount2a;
    //vector<delphi_integer>::iterator iepsvOdd     = iepsv.end() + icount2a; //This line causes errors

    delphi_integer icountEven = 0;
    delphi_integer icountOdd  = 0;

    /*
     * Now loop over all boundary points
     */
    for (size_t i = 0; i < mpi_size; i++)
    {
        mpi_ibgrd = pdc->readGlobalVector1D< SGrid<delphi_integer> >("ibgrd", i);

        ix = mpi_ibgrd.nX;
        iy = mpi_ibgrd.nY;
        iz = mpi_ibgrd.nZ;

        if (0 == iDirectEpsMap)
        {
            it.assign(6, 0);

            iw         = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);
            iepsmp_val = pdc->readGlobalVector1D< SGrid<delphi_integer> >("iepsmp", iw);
            if (0 != iepsmp_val.nX / iEpsDim) it[0] = 1;
            if (0 != iepsmp_val.nY / iEpsDim) it[1] = 1;
            if (0 != iepsmp_val.nZ / iEpsDim) it[2] = 1;

            iw         = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1 - 1);
            iepsmp_val = pdc->readGlobalVector1D< SGrid<delphi_integer> >("iepsmp", iw);
            if (0 != iepsmp_val.nX / iEpsDim) it[3] = 1;

            iw         = (iz - 1) * iGrid * iGrid + (iy - 1 - 1) * iGrid + (ix - 1);
            iepsmp_val = pdc->readGlobalVector1D< SGrid<delphi_integer> >("iepsmp", iw);
            if (0 != iepsmp_val.nY / iEpsDim) it[4] = 1;

            iw         = (iz - 1 - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);
            iepsmp_val = pdc->readGlobalVector1D< SGrid<delphi_integer> >("iepsmp", iw);
            if (0 != iepsmp_val.nZ / iEpsDim) it[5] = 1;

            ieps = it[0] + it[1] + it[2] + it[3] + it[4] + it[5];
        }
        else //Here is the key point for Gaussian
        {
            ieps = 0;
            temp = 0.0;
            delphi_real mpi_medeps;


            delphi_integer iIndex;

            if (iGaussian == 0)
            {
                iw         = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);
                iepsmp_val = pdc->readGlobalVector1D< SGrid<delphi_integer> >("iepsmp", iw);

                iIndex      = iepsmp_val.nX / iEpsDim;
                mpi_medeps  = pdc->readGlobalVector1D< delphi_real >("medeps", iIndex);
                temp       += mpi_medeps;
                vecttemp[0] = mpi_medeps;

                iIndex      = iepsmp_val.nY / iEpsDim;
                mpi_medeps  = pdc->readGlobalVector1D< delphi_real >("medeps", iIndex);
                temp       += mpi_medeps;
                vecttemp[1] = mpi_medeps;

                iIndex      = iepsmp_val.nZ / iEpsDim;
                mpi_medeps  = pdc->readGlobalVector1D< delphi_real >("medeps", iIndex);
                temp       += mpi_medeps;
                vecttemp[2] = mpi_medeps;

                iw          = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1 - 1);
                iepsmp_val  = pdc->readGlobalVector1D< SGrid<delphi_integer> >("iepsmp", iw);

                iIndex      = iepsmp_val.nX / iEpsDim;
                mpi_medeps  = pdc->readGlobalVector1D< delphi_real >("medeps", iIndex);
                temp       += mpi_medeps;
                vecttemp[3] = mpi_medeps;

                iw          = (iz - 1) * iGrid * iGrid + (iy - 1 - 1) * iGrid + (ix - 1);
                iepsmp_val  = pdc->readGlobalVector1D< SGrid<delphi_integer> >("iepsmp", iw);

                iIndex      = iepsmp_val.nY / iEpsDim;
                mpi_medeps  = pdc->readGlobalVector1D< delphi_real >("medeps", iIndex);
                temp       += mpi_medeps;
                vecttemp[4] = mpi_medeps;

                iw          = (iz - 1 - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);
                iepsmp_val  = pdc->readGlobalVector1D< SGrid<delphi_integer> >("iepsmp", iw);

                iIndex      = iepsmp_val.nZ / iEpsDim;
                mpi_medeps  = pdc->readGlobalVector1D< delphi_real >("medeps", iIndex);
                temp       += mpi_medeps;
                vecttemp[5] = mpi_medeps;
            }
            else if (iGaussian == 1) //Gaussian
            {
                SGrid<delphi_real> mpi_gepsmp2;

                iw          = (ix - 1) * iGrid * iGrid + (iy - 1) * iGrid + (iz - 1);
                mpi_gepsmp2 = pdc->readGlobalVector1D< SGrid<delphi_real> >("gepsmp2", iw);

                temp       += mpi_gepsmp2.nX;
                vecttemp[0] = mpi_gepsmp2.nX;

                temp       += mpi_gepsmp2.nY;
                vecttemp[1] = mpi_gepsmp2.nY;

                temp       += mpi_gepsmp2.nZ;
                vecttemp[2] = mpi_gepsmp2.nZ;

                delphi_integer iw3 = iw - iGrid * iGrid;
                mpi_gepsmp2 = pdc->readGlobalVector1D< SGrid<delphi_real> >("gepsmp2", iw3);

                temp       += mpi_gepsmp2.nX;
                vecttemp[3] = mpi_gepsmp2.nX;

                delphi_integer iw4 = iw - iGrid;
                mpi_gepsmp2 = pdc->readGlobalVector1D< SGrid<delphi_real> >("gepsmp2", iw4);

                temp       += mpi_gepsmp2.nY;
                vecttemp[4] = mpi_gepsmp2.nY;

                delphi_integer iw5 = iw - 1;
                mpi_gepsmp2 = pdc->readGlobalVector1D< SGrid<delphi_real> >("gepsmp2", iw5);

                temp       += mpi_gepsmp2.nZ;
                vecttemp[5] = mpi_gepsmp2.nZ;
            }
        }

        deb = 0;
        iw = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1);

        if (*(idebmap + iw))
        {
            deb = 1;
            idbs += 1;
        }

        delphi_real gridDensity;
        if (iGaussian==1)
        {
            //gridDensity = gaussianDensityMap[iw];
            gridDensity = pdc->readGlobalVector1D< delphi_real >("gdensity", iw);
        }

        vector < delphi_real > dbrow;

        if (0 == iDirectEpsMap)
        {
            dbrow.push_back(rgfDielecBndyValue[it[3]][ieps][deb]);
            dbrow.push_back(rgfDielecBndyValue[it[0]][ieps][deb]);
            dbrow.push_back(rgfDielecBndyValue[it[4]][ieps][deb]);
            dbrow.push_back(rgfDielecBndyValue[it[1]][ieps][deb]);
            dbrow.push_back(rgfDielecBndyValue[it[5]][ieps][deb]);
            dbrow.push_back(rgfDielecBndyValue[it[2]][ieps][deb]);
        }
        else
        {
            fDenom = temp + deb * fDebFct;

            if (0.0 == fIonStrength)
            {
                dbrow.push_back(vecttemp[3] / fDenom - fSixth);
                dbrow.push_back(vecttemp[0] / fDenom - fSixth);
                dbrow.push_back(vecttemp[4] / fDenom - fSixth);
                dbrow.push_back(vecttemp[1] / fDenom - fSixth);
                dbrow.push_back(vecttemp[5] / fDenom - fSixth);
                dbrow.push_back(vecttemp[2] / fDenom - fSixth);
            }
            else
            {
                dbrow.push_back(vecttemp[3] / fDenom);
                dbrow.push_back(vecttemp[0] / fDenom);
                dbrow.push_back(vecttemp[4] / fDenom);
                dbrow.push_back(vecttemp[1] / fDenom);
                dbrow.push_back(vecttemp[5] / fDenom);
                dbrow.push_back(vecttemp[2] / fDenom);
            }
        }

        vector < delphi_real > dbrow_original;
        dbrow_original.push_back(vecttemp[3]);
        dbrow_original.push_back(vecttemp[0]);
        dbrow_original.push_back(vecttemp[4]);
        dbrow_original.push_back(vecttemp[1]);
        dbrow_original.push_back(vecttemp[5]);
        dbrow_original.push_back(vecttemp[2]);

        iw = (iz - 1) * iGrid * iGrid + (iy - 1) * iGrid + (ix - 1); // iw=isgrid*(k-1) + igrid*(j-1) + i

        /*
         * Faster method than the commented lines above
         */
        if (0 == iw % 2) // iw is even
        {
            *(idpos.begin() + icountEven) = iw / 2 + 1;
            *(iepsv.begin() + icountEven) = ieps;
            *(db.begin() + icountEven) = dbrow;
            *(gaussianBoundaryDensity.begin() + icountEven) = gridDensity;
            *(gaussianBoundaryDielec.begin() + icountEven) = dbrow_original;
            icountEven++;
        }
        else // iw is odd
        {
            *(idpos.begin() + icount2a + icountOdd) = (iw + 1) / 2;
            *(iepsv.begin() + icount2a + icountOdd) = ieps;
            *(db.begin() + icount2a + icountOdd) = dbrow;
            *(gaussianBoundaryDensity.begin() + icount2a + icountOdd) = gridDensity;
            *(gaussianBoundaryDielec.begin() + icount2a + icountOdd) = dbrow_original;
            icountOdd++;
        }

        if (bDbOut)
        {
            ofDbFileStream << setw(3) << left << ix << " " << setw(3) << left << iy << " " << setw(3) << left << iz << " " << setw(8) << left << dbrow[0] << setw(8) << left << dbrow[1]
            << setw(8) << left << dbrow[2] << setw(8) << left << dbrow[3] << setw(8) << left << dbrow[4] << setw(8) << left << dbrow[5] << endl;
        }
    } //---------- end of for (delphi_integer ix = 0; ix < iBndyGridNum; ix += 1)

    if (bDbOut) ofDbFileStream.close();

    #ifdef VERBOSE
    cout << " Info> Number of dielectric boundary points in salt = " << idbs << endl;
    #endif

    /*
     * realign idpos and db,compressing to contingous space
     *
     * set saltmaps 1 and 2, i.e., sf1 and sf2.
     */
    if (0.0 < fIonStrength)
    {
        const delphi_real fSixSalt = fSixth * (1.0 / (1.0 + fDebFct / fSixEps) - 1.0); // sixsalt

        iw = 0;

        for (vector<char>::const_iterator const_itr = prgbDielecMap.cbegin(); const_itr != prgbDielecMap.cend(); ++const_itr)
        {
            if (0 == iw % 2) // even pts
            {
                if (*const_itr)
                    sf1.push_back(fSixth + fSixSalt);
                else
                    sf1.push_back(fSixth);
            }
            else // odd pts
            {
                if (*const_itr)
                    sf2.push_back(fSixth + fSixSalt);
                else
                    sf2.push_back(fSixth);
            }

            iw += 1;
        }

        sf2.push_back(0.0); //sf[iGrid3]

        for (int i = 0; i < icount2a; i++)
        {
            if (0 != iDirectEpsMap)
                sf1[idpos[i] - 1] = 0.0;
            else
            {
                if (abs(sf1[idpos[i] - 1] - fSixth) < fZero)
                    sf1[idpos[i] - 1] = rgfSaltFuncDiff[iepsv[idpos[i] - 1]][0];
                else
                    sf1[idpos[i] - 1] = rgfSaltFuncDiff[iepsv[idpos[i] - 1]][1];
            }
        }

        for (int i = icount2a; i < icount2b; i++)
        {
            if (0 != iDirectEpsMap)
                sf2[idpos[i] - 1] = 0.0;
            else
            {
                if (abs(sf2[idpos[i] - 1] - fSixth) < fZero)
                    sf2[idpos[i] - 1] = rgfSaltFuncDiff[iepsv[idpos[i] - 1]][0];
                else
                    sf2[idpos[i] - 1] = rgfSaltFuncDiff[iepsv[idpos[i] - 1]][1];
            }
        }
    }

    //---------------------------- end of mkdbsf.f -----------------------//
}

#endif
