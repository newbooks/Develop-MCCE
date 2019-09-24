/*
 * solver_fastSOR_nitit.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

void CDelphiFastSOR::nitit(const delphi_real& qfact) 
{
    delphi_integer ix, iy, iz;
    vector<delphi_real> rmsl(nxran, 0.0), rmaxl(nxran, 0.0);
    bool ichangeom, istop, inewfirst, inew;
    delphi_real relparprev, factor, fraction, derprec, der;
    int itr, itnum, ires, icountplus;
    delphi_real rmsch, rmxch;
    delphi_real conv[3] = { 0.0, 0.0, 0.0 };
    string nlstr;
    string strLine60 = " ----------------------------------------------------------------";

    if (0 == iConvergeFract) 
    {
        iIterateInterval = 10;
        iConvergeFract   = 1;
    }

    if (iLinIterateNum < iIterateInterval) iIterateInterval = iLinIterateNum;

    debmap1.assign(iHalfGridNum, 0.0);
    debmap2.assign(iHalfGridNum, 0.0);

    if (iGaussian != 0)
    {
        gaussianBoundaryNonlinear.assign(iDielecBndyOdd, 0.0);
        gaussianChargeNonlinear.assign(iCrgedGridSum, 0.0);
    }

    if (debug_solver) 
    {
        cout << "gaussianBoundaryDielec.size= "  << gaussianBoundaryDielec.size()  << endl;
        cout << "gaussianBoundaryDensity.size= " << gaussianBoundaryDensity.size() << endl;
        cout << "gaussianChargeDielec.size= "    << gaussianChargeDielec.size()    << endl;
        cout << "gaussianChargeDensity.size= "   << gaussianChargeDensity.size()   << endl;
    }

    for (ix = 0; ix < iHalfGridNum; ix++) 
    {
        iy = ix * 2;
        if (prgbDielecMap[iy])     debmap1[ix] = 1.0;
        if (prgbDielecMap[iy + 1]) debmap2[ix] = 1.0;
    }

    initOddEvenItr(2); // forWhom = 2

    cout << " Linear relaxation parameter" << " : " << om2 << endl;

    if (bManualRelaxParam) 
    {
        ichangeom = true;
        cout << " Non linear fixed relaxation parameter" << " : " << fRelaxParam << endl;
    } 
    else 
    {
        ichangeom = false;
        cout << " Non linear initial relaxation parameter " << " : " << fRelaxParam << endl;
        cout << " q factor " << " : " << qfact << endl;
    }

    cout << strLine60 << endl;
    //cout << "       rms-change     max change       #iterations" << endl;
    cout << "      " << " rms-change   max change       #iterations" << endl;
    cout << strLine60 << endl;

    #ifdef DEBUG_DELPHI_SOLVER
    {
        string strTestFile = "test_nitit.dat";
        ofstream ofTestStream(strTestFile.c_str());
        ofTestStream << boolalpha;
        ofTestStream << fixed << setprecision(7);

        ix = 1;
        for (vector<delphi_real>::iterator it = debmap1.begin(); it != debmap1.end(); ++it)
        {
            ofTestStream << "debmap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = debmap2.begin(); it != debmap2.end(); ++it)
        {
            ofTestStream << "debmap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = bndx1.begin(); it != bndx1.end(); ++it)
        {
            ofTestStream << "bndx1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = bndx2.begin(); it != bndx2.end(); ++it)
        {
            ofTestStream << "bndx2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = bndx3.begin(); it != bndx3.end(); ++it)
        {
            ofTestStream << "bndx3[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = bndx4.begin(); it != bndx4.end(); ++it)
        {
            ofTestStream << "bndx4[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

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

    bool testflag = true;
    #endif // DEBUG_DELPHI_SOLVER

    itr        = 1;
    itnum      = 0;
    ires       = 0;
    istop      = true;
    inewfirst  = false;
    inew       = true;
    factor     = 1.0;
    nlstr      = "                  ";
    icountplus = 0;
    relparprev = om2;
    fraction   = 0.0;

    do 
    {
        rmsch = 0.0;
        rmxch = 0.0;

        /*
         * iterate over odd points
         */
        itrOddPoints(2, itr); // forWhom = 2       
        
        #ifdef DEBUG_DELPHI_SOLVER
        if (false)
        {
            testflag = false;

            string strTestFile = "test_nitit.dat";
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
        #endif // DEBUG_DELPHI_SOLVER

        /*
         * Next update phimap2 using the new phimap1
         */
        itrEvenPoints(2, itr); // forWhom = 2

        #ifdef DEBUG_DELPHI_SOLVER
        if (false)
        {
            testflag = false;

            string strTestFile = "test_nitit.dat";
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
        #endif // DEBUG_DELPHI_SOLVER
        
        /*
         * we also save time by only checking convergence every 10 iterations, rather than every single iteration.
         * store phi2 in phi3 to compare against next iteration
         */
        if (iIterateInterval - 1 == itr % iIterateInterval) // itr = 9,19,29,...
        {
            for (ix = 1; ix < iHalfGridNum; ix += iConvergeFract)
                prgfPhiMap[ix] = phimap2[ix];
        }

        /*
         * check to see if accuracy is sufficient
         */
        if (0 == itr % iIterateInterval) 
        {
//            //if (2 == itr)
//            {
//                string strTestFile = "rank1_solver_nitit.dat";
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
//                ix = 0;
//                for (vector<delphi_real>::iterator it = prgfPhiMap.begin(); it != prgfPhiMap.end(); ++it)
//                {
//                    ofTestStream << "mpi_realbuff[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//                    ix++;
//                }
//
//                ofTestStream.close();
//                
//                return;
//            } 
            
            delphi_real rnorm2 = 0.0, temp2;

            for (ix = 1; ix < iHalfGridNum; ix += iConvergeFract) 
            {
                temp2   = prgfPhiMap[ix] - phimap2[ix];
                rnorm2 += temp2 * temp2;
                rmxch   = max(rmxch, abs(temp2));
            }

            conv[2] = conv[1];
            conv[1] = conv[0];
            conv[0] = rmxch;

            rmsch = sqrt((delphi_real) iConvergeFract * rnorm2 / ((iGrid - 2) * (iGrid - 2) * (iGrid - 2)));
            //rnormch = sqrt(rnorm2)/fRelaxParam;

            if ((fRmsc > rmsch || fMaxc > rmxch) && (22 < itnum)) ires = 1;

            if (0 == itnum) 
                //cout << scientific << rmsch << "  "  << rmxch << "  at  " << setw(5) << left << itr << " iterations\n";
                cout << "        " << scientific << rmsch << "  " << rmxch << "  at  " << setw(5) << left << itr << " iterations" << endl;

            istop = !(0.22 < rmxch);

            if (bLogGraph) 
            {
                int ibin;
                for (int j = itr - 9; j <= itr; j++) 
                {
                    ibin = (j - 1) * (60 - 1) / (iLinIterateNum - 1) + 1;
                    rmsl[ibin - 1]  = rmsch;
                    rmaxl[ibin - 1] = rmxch;
                }
            }

            //----- optimization of the first step
            inewfirst = (itr > iLinIterateNum - iLinIterateNum % iIterateInterval) && istop && (0 == itnum);

            if (!bManualRelaxParam && inewfirst && (3.2 < qfact)) 
            {
                factor    = exp(-qfact * 2.1) + fZero;
                ichangeom = true;
            }

            //----- nonlinear part
            inew = inew || inewfirst;

            if (0 != itnum) 
            {
                //cout << scientific << rmsch << " " << rmxch << " " << itnum << " it. " << nlstr << endl;
                cout << "        " << scientific << rmsch << "  " << rmxch << "  at  " << setw(5) << left << itnum << " iterations" 
                     << nlstr << endl;

                if (!bManualRelaxParam) 
                {
                    derprec = der;
                    der     = (conv[0] - conv[1]) / conv[1];

                    if (fZero > rmxch) 
                    {
                        factor = 1.2;
                        ichangeom = true;
                    } 
                    else 
                    {
                        if (0.0 < der && !inew) 
                        {
                            icountplus++;
                            factor = factor * pow(1.0 - der, 0.99);
                            
                            if (0.55 < der) 
                            {
                                ichangeom = true;
                                if (1.0 <= der) factor = 1.0e-5;
                            }
                            
                            if (0.35 < der && 0.1 < conv[0]) 
                            {
                                ichangeom = true;
                                factor = pow(factor * 0.05 / conv[0], 4.0);
                            }
                        }

                        if ((0.0 < der && 0.1 < conv[0]) && inew && !inewfirst) 
                        {
                            ichangeom = true;
                            factor    = min(factor * 0.05 / conv[0], factor * pow(1.0 - der, 0.86));
                        }

                        if (0.0 >= der) 
                        {
                            icountplus = 0;
                            factor = 1.0;
                            if (24 < itnum && itnum < 24 + 0.75 * (iNonIterateNum - 24) && 0.03 > rmxch && 0.0 >= derprec) 
                            {
                                if (-0.2 < der && -0.2 < derprec) 
                                {
                                    #ifdef VERBOSE
                                    cout << "Trying to speed up the convergence process\n";
                                    #endif
                                    
                                    factor = 1.1;
                                    ichangeom = true;
                                    if (0.2 > fRelaxParam && -0.05 < der && -0.05 < derprec) factor = 1 - 45.226 * (fRelaxParam - 0.2);
                                }
                            }
                        }

                        if (2 <= icountplus) ichangeom = true;
                    }

                    inewfirst = false;
                }
            }
        } //----- end of convergence check

        itr++;

        if ((iLinIterateNum >= itr || !istop) && fZero > abs(ires)) continue;

        if (0 < iNonIterateNum && 0 == itnum) 
        {
            cout << "\n  Now for the non-linear iterations" << endl;
            //cout << "\n       rms-change     max change         #iterations\n";
            cout << strLine60 << endl;
            cout << "      " << " rms-change   max change       #iterations" << endl;
            cout << strLine60 << endl;
        }

        iIterateInterval = 10; //----- icon1 = how many blocks each iteration convergence occurs
        iLinIterateNum   = 10; //----- nlit  = How many iterations in the block
        itnum++;

        if (iNonIterateNum < itnum || (1 == ires && !inew)) break;

        itr = 1;
        inew = ichangeom;
        if (ichangeom) 
        {
            delphi_real omcomp;
            fRelaxParam = fRelaxParam * factor;

            if (1.0e-4 > fRelaxParam) 
            {
                #ifdef VERBOSE
                cout << "estimation " << fRelaxParam << " 1E-4 preferred\n";
                #endif
                
                fRelaxParam = 1.0e-4;
            }

            factor = 1.0;

            #ifdef VERBOSE
            cout << " New relaxation parameter" << " : " << fRelaxParam << endl;
            #endif
            
            ichangeom  = false;
            icountplus = 0;
            omcomp     = fRelaxParam / relparprev;
            relparprev = fRelaxParam;
            om1        = 1.0 - fRelaxParam;

            for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
                *it = (*it) * omcomp;

            for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it)
                *it = (*it) * omcomp;

            for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)
                *it = (*it) * omcomp;

            for (delphi_integer iy = 0; iy < iDielecBndyOdd; iy++)
                for (delphi_integer ix = 0; ix < 6; ix++)
                    prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix] * omcomp;

            sixth = sixth * omcomp;
        }

        fraction += 0.05;

        if (1.0 < fraction) 
        {
            fraction = 1.0;
            nlstr    = " Full non-linearity";
        }

        /*
         * Here calculates the boundary and charge pure nonlinear for Gaussian based boundary
         */
        if (iGaussian == 1)
        {
            /*
             * The Even and Odd boundary points
             */
        	
            //for Odd boundary points
            for (int n = 0; n < iDielecBndyEven; n++) 
            {
                ix = prgiBndyDielecIndex[n];
                delphi_real myDensity        = gaussianBoundaryDensity[n];
                delphi_real myExpSolvE       = calcExpSolvE(myDensity);
                delphi_real myPhi            = phimap1[ix - 1];
                delphi_real myEpsSum         = gaussianBoundaryDielec[n][0] + gaussianBoundaryDielec[n][1] + gaussianBoundaryDielec[n][2] 
                                             + gaussianBoundaryDielec[n][3] + gaussianBoundaryDielec[n][4] + gaussianBoundaryDielec[n][5];
                delphi_real my_nonlinear     = calcPhiMinusSinh(myPhi) / (myEpsSum / (myExpSolvE * fDebFct) / fEPKT + 1);
                gaussianBoundaryNonlinear[n] = fraction * my_nonlinear;
            }

            //for Even boundary points
            for (int n = iDielecBndyEven; n < iDielecBndyOdd; n++) 
            {
                ix = prgiBndyDielecIndex[n];
                delphi_real myDensity        = gaussianBoundaryDensity[n];
                delphi_real myExpSolvE       = calcExpSolvE(myDensity);
                delphi_real myPhi            = phimap2[ix - 1];
                delphi_real myEpsSum         = gaussianBoundaryDielec[n][0] + gaussianBoundaryDielec[n][1] + gaussianBoundaryDielec[n][2] 
                                             + gaussianBoundaryDielec[n][3] + gaussianBoundaryDielec[n][4] + gaussianBoundaryDielec[n][5];
                delphi_real my_nonlinear     = calcPhiMinusSinh(myPhi) / (myEpsSum / (myExpSolvE * fDebFct) / fEPKT + 1);
                gaussianBoundaryNonlinear[n] = fraction * my_nonlinear;
            }

            /*
             * The pure nonlinear part for even and odd charge points
             */

            //for Odd charged points
            for (int n = 0; n < iCrgedGridEven; n++) 
            {
                ix = prgiCrgPose[n];
                delphi_real myDensity      = gaussianChargeDensity[n];
                delphi_real myExpSolvE     = calcExpSolvE(myDensity);
                delphi_real myPhi          = phimap1[ix - 1];
                delphi_real myEpsSum       = gaussianChargeDielec[n][0] + gaussianChargeDielec[n][1] + gaussianChargeDielec[n][2] 
                                           + gaussianChargeDielec[n][3] + gaussianChargeDielec[n][4] + gaussianChargeDielec[n][5];
                delphi_real my_nonlinear   = calcPhiMinusSinh(myPhi) / (myEpsSum / (myExpSolvE * fDebFct) / fEPKT + 1);
                gaussianChargeNonlinear[n] = fraction * my_nonlinear;
            }

            //for Even charged points
            for (int n = iCrgedGridEven; n < iCrgedGridSum; n++) 
            {
                ix = prgiCrgPose[n];
                delphi_real myDensity      = gaussianChargeDensity[n];
                delphi_real myExpSolvE     = calcExpSolvE(myDensity);
                delphi_real myPhi          = phimap2[ix - 1];
                delphi_real myEpsSum       = gaussianChargeDielec[n][0] + gaussianChargeDielec[n][1] + gaussianChargeDielec[n][2] 
                                           + gaussianChargeDielec[n][3] + gaussianChargeDielec[n][4] + gaussianChargeDielec[n][5];
                delphi_real my_nonlinear   = calcPhiMinusSinh(myPhi) / (myEpsSum / (myExpSolvE * fDebFct) / fEPKT + 1);
                gaussianChargeNonlinear[n] = fraction * my_nonlinear;
            }
        } // End of if Gaussian
        else // if not Gaussian
        {
            delphi_real temp1, temp2;
            delphi_real fac1 = fraction * fDebFct / (2.0 * fIonStrength * fEpsOut);

            for (ix = 0; ix < iHalfGridNum; ix++) 
            {
                temp1     = phimap1[ix] * debmap1[ix];
                temp2     = phimap2[ix] * debmap2[ix];
                qmap1[ix] = fac1 * calcPhiMinusSinh(temp1);
                qmap2[ix] = fac1 * calcPhiMinusSinh(temp2);
            }
        }
    } while (true);

    if (0.05 > fRelaxParam) CSmallRelaxParam waring;

    postItr(rmaxl, rmsl);

    #ifdef DEBUG_DELPHI_SOLVER
    if (false)
    {
        string strTestFile = "test_nitit.dat";
        ofstream ofTestStream(strTestFile.c_str());
        ofTestStream << boolalpha;
        ofTestStream << fixed << setprecision(7);

        const delphi_real *** phimap = pdc->getKey_constPtr<delphi_real>("phimap",iGrid,iGrid,iGrid); // const pointer to 3D phimap
        for (iz = 0; iz < iGrid; iz++)
        {
            for (iy = 0; iy < iGrid; iy++)
            {
                for (ix = 0; ix < iGrid; ix++)
                {
                    ofTestStream << "phimap[" << setw(6) << right << iz+1 << "," << setw(6) << right << iy+1 
                                 << "," << setw(6) << right << ix+1 << "] = " << setw(11) << right << phimap[iz][iy][ix] << endl;
                }
            }
        }

        ofTestStream.close();
    }
    #endif // DEBUG_DELPHI_SOLVER
}

/*
 This function calculates the exponential of solvation energy at a given grid point with gaussian solute density gdens.
 The solvation energy is caused by the difference between dielectric constant at this grid point and at bulk solvent.
 */
delphi_real CDelphiFastSOR::calcExpSolvE(delphi_real gdens) 
{
    delphi_real result = 0;
    delphi_real repsdens;
    delphi_real solvationEnergy;
    delphi_real halfSternRadiusInverse;

    repsdens               = gdens * repsin + (1 - gdens) * repsout;
    halfSternRadiusInverse = 1 / 2.0 * 0.5;
    solvationEnergy        = -fEPKT * halfSternRadiusInverse * (1 / repsdens - 1 / repsout);

    result = exp(solvationEnergy);

    return result;
}

/*
 This function returns sinh(x) using Taylor expansion. 
 A cut off is applied for input x to prevent extremely large results.
 */
delphi_real CDelphiFastSOR::calcSinh(delphi_real x) 
{
    //cut off for faster converge
    delphi_real cutOff = 20;

    if (fabs(x) > cutOff) 
    {
        x = (x > 0 ? cutOff : -cutOff);
    }

    //-------Taylor Series 5--------
    delphi_real x2 = x * x;
    delphi_real result = x * (1 + x2 * (1 / 6));  //Taylor 3
    //delphi_real result = x *(1 + x2 *(1/6 + phi2 / 120)); // Taylor 5

    // sinh
    //delphi_real result = sinh(x);

    return result;
}

/*
 This function calculates the non-linear sinh function used in DelPhi.
 It returns (1-sinh(x))/x.
 This calculation is usually done by a 5th order Taylor expansion.
 But it can alternatively done by the original sinh function or 3rd order Taylor expansion.
 A cut-off value of 20 is set for input x values to prevent the divergence caused by extremely large results.
 */
delphi_real CDelphiFastSOR::calcPhiMinusSinh(delphi_real x) 
{
    //cut off for faster converge
    delphi_real cutOff = 20;

    if (fabs(x) > cutOff) 
    {
        x = (x > 0 ? cutOff : -cutOff);
    }

    //-------Taylor Series 5--------
    delphi_real x2 = x * x;

    //Taylor 3
    //delphi_real result = x2 *(fTaylorCoeff2 + x * fTaylorCoeff3 ); 

    // Taylor 5
    delphi_real result = x2 * (fTaylorCoeff2 + x * (fTaylorCoeff3 + x * (fTaylorCoeff4 + x * fTaylorCoeff5)));

    // sinh
    //delphi_real result = (1 - sinh(x))/x;

    return result;
}

/*
 This function calculates the non-linear exp function used in DelPhi.
 It returns exp(x).
 This calculation is usually done by a 5th order Taylor expansion.
 But it can alternatively done by the original exp function or 3rd order Taylor expansion.
 */
delphi_real CDelphiFastSOR::calcExp(delphi_real x) 
{
    delphi_real x2     = x * x;
    delphi_real x3     = x * x2;
    delphi_real x4     = x2 * x2;
    delphi_real x5     = x2 * x3;
    delphi_real result = 1 + x + (x2 / 2) + (x3 / 6) + (x4 / 24) + (x5 / 120);
    
    return result;
}

#ifdef PARALLEL_MPI

/*
 * subroutine mpi_nitit(qfact)
 */
void CDelphiFastSOR::mpi_nitit(const delphi_real& qfact)
{
    delphi_integer ix, iy, iz;
    vector<delphi_real> rmsl(nxran, 0.0), rmaxl(nxran, 0.0);
    bool ichangeom, istop, inewfirst, inew;
    delphi_real relparprev, factor, fraction, derprec, der;
    int itr, itnum, ires, icountplus;
    delphi_real rmsch, rmxch;
    delphi_real conv[3] = { 0.0, 0.0, 0.0 };
    string nlstr;
    string strLine60 = " ----------------------------------------------------------------";
    
    vector<delphi_real> realbuff;
    delphi_real *mpi_realbuff;
       
    if (0 == iConvergeFract)
    {
        iIterateInterval = 10;
        iConvergeFract   = 1;
    }

    if (iLinIterateNum < iIterateInterval) iIterateInterval = iLinIterateNum;
    
    /* 
     * on master: debmap1-2 are of size iHalfGridNum; on slaves: debmap1-2 are of the same size of phimap1-2
     */
    if (0 == mpi_rank) /* master process */
    {
        debmap1.assign(iHalfGridNum, 0.0);
        debmap2.assign(iHalfGridNum, 0.0);
        
        for (ix = 0; ix < iHalfGridNum; ix++) 
        {
            iy = ix * 2;
            if (prgbDielecMap[iy])     debmap1[ix] = 1.0;
            if (prgbDielecMap[iy + 1]) debmap2[ix] = 1.0;
        }

        /*
         * master frees prgbDielecMap (a.k.a idebmap(igrid,igrid,igrid)) in local data container
         */
        if (prgbDielecMap.size() > 0) vector<char>().swap(prgbDielecMap);
    }
    
    mpi_initOddEvenItr(2); // forWhom = 2, 
                           // distributing debmap1-2 and debmap2 to all slave processes is moved here

    /*
     * Gaussian based runs are not parallelized yet!
     */
//    if (iGaussian != 0)
//    {
//        gaussianBoundaryNonlinear.assign(iDielecBndyOdd, 0.0);
//        gaussianChargeNonlinear.assign(iCrgedGridSum, 0.0);
//    }
    
    
    //MPI_Bcast(&bManualRelaxParam, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    if (bManualRelaxParam)
        ichangeom = true;
    else
        ichangeom = false;

    if (0 == mpi_rank) /* master process */
    {
        cout << " Linear relaxation parameter" << " : " << om2 << endl;

        if (bManualRelaxParam)
        {
            cout << " Non linear fixed relaxation parameter" << " : " << fRelaxParam << endl;
        }
        else
        {
            cout << " Non linear initial relaxation parameter " << " : " << fRelaxParam << endl;
            cout << " q factor " << " : " << qfact << endl;
        }

        cout << strLine60 << endl;
        cout << "      " << " rms-change   max change       #iterations" << endl;
        cout << strLine60 << endl;
    }

    itr        = 1;
    itnum      = 0;
    ires       = 0;
    istop      = true;
    inewfirst  = false;
    inew       = true;
    factor     = 1.0;
    nlstr      = "                  ";
    icountplus = 0;
    relparprev = om2;
    fraction   = 0.0;

    do
    {
        rmsch = 0.0;
        rmxch = 0.0;

        /*
         * iterate over odd points
         */
        mpi_itrOddPoints(2, itr); // forWhom = 2
        
        /*
         * Synchronization before getting into the loop
         */
        MPI_Barrier (MPI_COMM_WORLD);

        /*
         * Next update phimap2 using the new phimap1
         */
        mpi_itrEvenPoints(2, itr); // forWhom = 2
        
        /*
         * if not periodic bdy
         */
        if ( !(rgbPeriodicBndy[0] || rgbPeriodicBndy[1] || rgbPeriodicBndy[2]) )
        {
            delphi_real mpi_sumall;

            if (0 == mpi_rank) /* master process */
            {
                /*
                 * check to see if accuracy is sufficient
                 */
                if (0 == itr % iIterateInterval)
                {
//                    //if (2 == itr)
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
//                        if (0 == mpi_rank) strTestFile = "rank0_solver_nitit.dat";
//                        if (1 == mpi_rank) strTestFile = "rank1_solver_nitit.dat";
//                        if (2 == mpi_rank) strTestFile = "rank2_solver_nitit.dat";
//                        if (3 == mpi_rank) strTestFile = "rank3_solver_nitit.dat";
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
                    
                    delphi_real rnorm2 = 0.0, temp2;

                    MPI_Allreduce(&rnorm2, &mpi_sumall, 1, mpi_delphi_real, MPI_SUM, MPI_COMM_WORLD);
                    rnorm2 = mpi_sumall;

                    MPI_Allreduce(&rmxch, &mpi_sumall, 1, mpi_delphi_real, MPI_MAX, MPI_COMM_WORLD);
                    rmxch = mpi_sumall;

                    conv[2] = conv[1];
                    conv[1] = conv[0];
                    conv[0] = rmxch;

                    rmsch = sqrt( (delphi_real) iConvergeFract * rnorm2 / ((iGrid - 2) * (iGrid - 2) * (iGrid - 2)) );

                    if ((fRmsc > rmsch || fMaxc > rmxch) && (22 < itnum)) ires = 1;

                    if (0 == itnum)
                        cout << "        " << scientific << rmsch << "  " << rmxch << "  at  " << setw(5) << left << itr 
                             << " iterations" << endl;

                    istop = !(0.22 < rmxch);

                    if (bLogGraph)
                    {
                        int ibin;
                        for (int j = itr - 9; j <= itr; j++)
                        {
                            ibin = (j - 1) * (60 - 1) / (iLinIterateNum - 1) + 1;
                            rmsl[ibin - 1]  = rmsch;
                            rmaxl[ibin - 1] = rmxch;
                        }
                    }

                    //----- optimization of the first step
                    inewfirst = (itr > iLinIterateNum - iLinIterateNum % iIterateInterval) && istop && (0 == itnum);

                    if (!bManualRelaxParam && inewfirst && (3.2 < qfact))
                    {
                        factor = exp(-qfact * 2.1) + fZero;
                        ichangeom = true;
                    }

                    //----- nonlinear part
                    inew = inew || inewfirst;

                    if (0 != itnum)
                    {
                        cout << "        " << scientific << rmsch << "  " << rmxch << "  at  " << setw(5) << left << itnum 
                             << " iterations" << nlstr << endl;;

                        if (!bManualRelaxParam)
                        {
                            derprec = der;
                            der = (conv[0] - conv[1]) / conv[1];

                            if (fZero > rmxch)
                            {
                                factor = 1.2;
                                ichangeom = true;
                            }
                            else
                            {
                                if (0.0 < der && !inew)
                                {
                                    icountplus++;
                                    factor = factor * pow(1.0 - der, 0.99);
                                    if (0.55 < der)
                                    {
                                        ichangeom = true;
                                        if (1.0 <= der) factor = 1.0e-5;
                                    }
                                    if (0.35 < der && 0.1 < conv[0])
                                    {
                                        ichangeom = true;
                                        factor = pow(factor * 0.05 / conv[0], 4.0);
                                    }
                                }

                                if ((0.0 < der && 0.1 < conv[0]) && inew && !inewfirst)
                                {
                                    ichangeom = true;
                                    factor = min(factor * 0.05 / conv[0], factor * pow(1.0 - der, 0.86));
                                }

                                if (0.0 >= der)
                                {
                                    icountplus = 0;
                                    factor = 1.0;
                                    if (24 < itnum && itnum < 24 + 0.75 * (iNonIterateNum - 24) && 0.03 > rmxch && 0.0 >= derprec)
                                    {
                                        if (-0.2 < der && -0.2 < derprec)
                                        {
                                            #ifdef VERBOSE
                                            cout << "Trying to speed up the convergence process\n";
                                            #endif
                                            factor = 1.1;
                                            ichangeom = true;
                                            if (0.2 > fRelaxParam && -0.05 < der && -0.05 < derprec) factor = 1 - 45.226 * (fRelaxParam - 0.2);
                                        }
                                    }
                                }

                                if (2 <= icountplus) ichangeom = true;

                            } /* end of if (fZero > rmxch) */

                            inewfirst = false;

                        } /* end of if (!bManualRelaxParam) */

                    } /* end of if (0 != itnum) */

                } /* end of if (0 == itr % iIterateInterval) */

            } /* end of running on master process */

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

                /*
                 * check to see if accuracy is sufficient
                 */
                if (0 == itr % iIterateInterval)
                {
//                    //if (2 == itr)
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
//                        if (0 == mpi_rank) strTestFile = "rank0_solver_nitit.dat";
//                        if (1 == mpi_rank) strTestFile = "rank1_solver_nitit.dat";
//                        if (2 == mpi_rank) strTestFile = "rank2_solver_nitit.dat";
//                        if (3 == mpi_rank) strTestFile = "rank3_solver_nitit.dat";
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
//                        ix = mpi_wrphimap2_start;
//                        for (vector<delphi_real>::iterator it = realbuff.begin(); it != realbuff.end(); ++it)
//                        {
//                            ofTestStream << "mpi_realbuff[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
//                            ix++;
//                        }
//
//                        ofTestStream.close();
//                        
//                        return;
//                    }        
                    
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

                    conv[2] = conv[1];
                    conv[1] = conv[0];
                    conv[0] = rmxch;

                    rmsch = sqrt( (delphi_real) iConvergeFract * rnorm2 / ((iGrid - 2) * (iGrid - 2) * (iGrid - 2)));

                    if ((fRmsc > rmsch || fMaxc > rmxch) && (22 < itnum)) ires = 1;

                    istop = !(0.22 < rmxch);

                    //----- optimization of the first step
                    inewfirst = (itr > iLinIterateNum - iLinIterateNum % iIterateInterval) && istop && (0 == itnum);

                    if (!bManualRelaxParam && inewfirst && (3.2 < qfact))
                    {
                        factor = exp(-qfact * 2.1) + fZero;
                        ichangeom = true;
                    }

                    //----- nonlinear part
                    inew = inew || inewfirst;

                    if (0 != itnum)
                    {
                        if (!bManualRelaxParam)
                        {
                            derprec = der;
                            der = (conv[0] - conv[1]) / conv[1];

                            if (fZero > rmxch)
                            {
                                factor = 1.2;
                                ichangeom = true;
                            }
                            else
                            {
                                if (0.0 < der && !inew)
                                {
                                    icountplus++;
                                    factor = factor * pow(1.0 - der, 0.99);
                                    if (0.55 < der)
                                    {
                                        ichangeom = true;
                                        if (1.0 <= der) factor = 1.0e-5;
                                    }
                                    if (0.35 < der && 0.1 < conv[0])
                                    {
                                        ichangeom = true;
                                        factor = pow(factor * 0.05 / conv[0], 4.0);
                                    }
                                }

                                if ((0.0 < der && 0.1 < conv[0]) && inew && !inewfirst)
                                {
                                    ichangeom = true;
                                    factor = min(factor * 0.05 / conv[0], factor * pow(1.0 - der, 0.86));
                                }

                                if (0.0 >= der)
                                {
                                    icountplus = 0;
                                    factor = 1.0;
                                    if (24 < itnum && itnum < 24 + 0.75 * (iNonIterateNum - 24) && 0.03 > rmxch && 0.0 >= derprec)
                                    {
                                        if (-0.2 < der && -0.2 < derprec)
                                        {
                                            factor = 1.1;
                                            ichangeom = true;
                                            if (0.2 > fRelaxParam && -0.05 < der && -0.05 < derprec) factor = 1 - 45.226 * (fRelaxParam - 0.2);
                                        }
                                    }
                                }

                                if (2 <= icountplus) ichangeom = true;

                            } /* end of if (fZero > rmxch) */

                            inewfirst = false;

                        } /* end of if (!bManualRelaxParam) */

                    } /* end of if (0 != itnum) */

                } /* end of if (0 == itr % iIterateInterval) */

            } /* end of running on slave processes */

            itr++;

            if ((iLinIterateNum >= itr || !istop) && fZero > abs(ires)) continue;

        } /* end of  non-periodic bdy */

        if (0 == mpi_rank) /* master process */
        {
            if (0 < iNonIterateNum && 0 == itnum)
            {
                cout << "\n  Now for the non-linear iterations" << endl;
                cout << strLine60 << endl;
                cout << "      " << " rms-change   max change       #iterations" << endl;
                cout << strLine60 << endl;
            }
        }

        iIterateInterval = 10; //----- icon1 = how many blocks each iteration convergence occurs
        iLinIterateNum   = 10; //----- nlit  = How many iterations in the block
        itnum++;

        if (iNonIterateNum < itnum || (1 == ires && !inew)) break;

        itr = 1;
        inew = ichangeom;

        if (ichangeom)
        {
            delphi_real omcomp;
            fRelaxParam = fRelaxParam * factor;

            if (1.0e-4 > fRelaxParam)
            {
                if (0 == mpi_rank) /* master process */
                {
                    #ifdef VERBOSE
                    cout << "estimation " << fRelaxParam << " 1E-4 preferred\n";
                    #endif
                }
                fRelaxParam = 1.0e-4;
            }

            factor = 1.0;

            if(0 == mpi_rank) /* master process */
            {
                #ifdef VERBOSE
                cout << " New relaxation parameter" << " : " << fRelaxParam << endl;
                #endif
            }

            ichangeom = false;
            icountplus = 0;
            omcomp = fRelaxParam / relparprev;
            relparprev = fRelaxParam;
            om1 = 1.0 - fRelaxParam;

            for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
                *it = (*it) * omcomp;

            for (vector<delphi_real>::iterator it = prgfSaltMap2.begin();it != prgfSaltMap2.end(); ++it)
                *it = (*it) * omcomp;

            for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)
                *it = (*it) * omcomp;
            
            for (delphi_integer iy = 0; iy < prgfBndyDielec.size(); iy++)
                for (delphi_integer ix = 0; ix < 6; ix++)
                    prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix] * omcomp;

            sixth = sixth * omcomp;
        } /* end of if (ichangeom) */

        fraction += 0.05;
        
        if (1.0 < fraction)
        {
            fraction = 1.0;
            nlstr = " Full non-linearity";
        }

        /*
         * Here calculates the boundary and charge pure nonlinear for Gaussian based boundary
         * 
         * Gaussian based runs are not parallelized yet!
         */
        if (iGaussian == 1)
        {
//            /*
//             * The Even and Odd boundary points
//             */
//        	
//        	    //for Odd boundary points
//            for (int n = 0; n < iDielecBndyEven; n++) 
//            {
//                ix = prgiBndyDielecIndex[n];
//                delphi_real myDensity  = gaussianBoundaryDensity[n];
//                delphi_real myExpSolvE = calcExpSolvE(myDensity);
//                delphi_real myPhi      = phimap1[ix - 1];
//                delphi_real myEpsSum   = gaussianBoundaryDielec[n][0] + gaussianBoundaryDielec[n][1] + gaussianBoundaryDielec[n][2] + gaussianBoundaryDielec[n][3] + gaussianBoundaryDielec[n][4] + gaussianBoundaryDielec[n][5];
//                delphi_real my_nonlinear = calcPhiMinusSinh(myPhi) / (myEpsSum / (myExpSolvE * fDebFct) / fEPKT + 1);
//                gaussianBoundaryNonlinear[n] = fraction * my_nonlinear;
//            }
//
//            //for Even boundary points
//            for (int n = iDielecBndyEven; n < iDielecBndyOdd; n++) 
//            {
//                ix = prgiBndyDielecIndex[n];
//                delphi_real myDensity = gaussianBoundaryDensity[n];
//                delphi_real myExpSolvE = calcExpSolvE(myDensity);
//                delphi_real myPhi = phimap2[ix - 1];
//                delphi_real myEpsSum = gaussianBoundaryDielec[n][0] + gaussianBoundaryDielec[n][1] + gaussianBoundaryDielec[n][2] + gaussianBoundaryDielec[n][3] + gaussianBoundaryDielec[n][4] + gaussianBoundaryDielec[n][5];
//                delphi_real my_nonlinear = calcPhiMinusSinh(myPhi) / (myEpsSum / (myExpSolvE * fDebFct) / fEPKT + 1);
//                gaussianBoundaryNonlinear[n] = fraction * my_nonlinear;
//            }
//
//            /*
//             * The pure nonlinear part for even and odd charge points
//             */
//
//            //for Odd charged points
//            for (int n = 0; n < iCrgedGridEven; n++) 
//            {
//                ix = prgiCrgPose[n];
//                delphi_real myDensity = gaussianChargeDensity[n];
//                delphi_real myExpSolvE = calcExpSolvE(myDensity);
//                delphi_real myPhi = phimap1[ix - 1];
//                delphi_real myEpsSum = gaussianChargeDielec[n][0] + gaussianChargeDielec[n][1] + gaussianChargeDielec[n][2] + gaussianChargeDielec[n][3] + gaussianChargeDielec[n][4] + gaussianChargeDielec[n][5];
//                delphi_real my_nonlinear = calcPhiMinusSinh(myPhi) / (myEpsSum / (myExpSolvE * fDebFct) / fEPKT + 1);
//                gaussianChargeNonlinear[n] = fraction * my_nonlinear;
//            }
//
//            //for Even charged points
//            for (int n = iCrgedGridEven; n < iCrgedGridSum; n++) 
//            {
//                ix = prgiCrgPose[n];
//                delphi_real myDensity = gaussianChargeDensity[n];
//                delphi_real myExpSolvE = calcExpSolvE(myDensity);
//                delphi_real myPhi = phimap2[ix - 1];
//                delphi_real myEpsSum = gaussianChargeDielec[n][0] + gaussianChargeDielec[n][1] + gaussianChargeDielec[n][2] + gaussianChargeDielec[n][3] + gaussianChargeDielec[n][4] + gaussianChargeDielec[n][5];
//                delphi_real my_nonlinear = calcPhiMinusSinh(myPhi) / (myEpsSum / (myExpSolvE * fDebFct) / fEPKT + 1);
//                gaussianChargeNonlinear[n] = fraction * my_nonlinear;
//            }        	
        } // End of if Gaussian
        else // if not Gaussian
        {
            if (0 != mpi_rank) /* slave processes */
            {
                delphi_real temp1, temp2;
                delphi_real fac1 = fraction * fDebFct / (2.0 * fIonStrength * fEpsOut);

                for (ix = mpi_wrstar1[mpi_rank] - 1; ix < mpi_wrfinl1[mpi_rank]; ix++)
                {
                    temp1         = mpi_phimap1[ix] * mpi_debmap1[ix];             
                    mpi_qmap1[ix] = fac1 * calcPhiMinusSinh(temp1);
                }

                for (ix = mpi_wrstar2[mpi_rank] - 1; ix < mpi_wrfinl2[mpi_rank]; ix++)
                {
                    temp2         = mpi_phimap2[ix] * mpi_debmap2[ix];
                    mpi_qmap2[ix] = fac1 * calcPhiMinusSinh(temp2);
                }
            }
        } //  end of calculating the boundary and charge pure nonlinear for Gaussian based boundary
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

    if (0 == mpi_rank) /* master process */
    {
        if (0.05 > fRelaxParam)  CSmallRelaxParam waring;
        postItr(rmaxl, rmsl);
    }
}

#endif
