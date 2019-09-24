#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"
#include "space.h"

using namespace std;

void CDelphiSpace::crgarr()
{
    SGrid <delphi_real> rxyz;
    SGrid <delphi_integer> jxyz;

    //2011-05-30 Declarations added due to IMPLICIT NONE
    delphi_integer ic1,ic2,i,ix,ii,imed,jx,jy,jz,ridx,cidx,delta=0;
    delphi_real chrg,rgrid;
    bool ipassed;
    SGridValue<delphi_real> temp;
    delphi_real radpolext=1.0;

    // STRING FORMATTING VARIABLES
    const char * infoString = " Info> ";
    size_t MAXWIDTH = 45;

    //VARIABLES FOR SURFPOT-HIGHER ORDER MULTIPOLES
    SGrid <delphi_real> r;
    delphi_real xyz[3] = { 0., 0., 0. };
    delphi_real r2=0;

    /*
     * 2011-05-30 First determined ic1 = nqass - number of assigned charges
     * in order to avoid re-sizing of arrays atmcrg,chgpos and crgatn.
     */

    // Part1 (easy) : from molecules
    temp.nGrid.nX=0.; temp.nGrid.nY=0.; temp.nGrid.nZ=0.; temp.nValue=0; //temp.nGrid={0.,0.,0.};

    qMoment.assign(9,0.0);
    dMoment.assign(3,0.0);

    ic1=0;
    for(i=1; i<=iNatom; i++)
    {
        // (Parallel) check if atom is in the net grid
        if (checkNotBufferCoor(sDelPhiPDB[i].xyz) || ibctyp == 3)
            if(abs(sDelPhiPDB[i].charge)>1.e-6)
                ic1=ic1+1;
    }

    /*
     * 2011-05-30 Part 2 (more difficult) : from charge distributions
     * running distrTOpointeger without assigning values to  arrays not allocated yet
     */
    #ifdef VERBOSE
    if (ndistr>0) cout << "Warning: call distrTOpoFloat2Int(ic1,false)..." << endl;
    #endif

    nqass=ic1;

    /*
     * atmcrg contains grid positions of all charges AND the charge in the 4th field atmeps6 contains
     * 6*epsilon/epkt as a function of ic2-th charge internal atmeps6 to the grid - NO LONGER USED
     * nqgrdtonqass maps ic2 to ic1 atmeps contains epsilon/epkt as a funcion of ic1-th general
     * charge
     */
    atmcrg_v.assign(nqass,temp);
    chgpos_v.assign(nqass,sgrid_temp_real);

    crgatn_v.assign(nqass,0);
    atmeps_v.assign(nqass,0.);
    atmcrg=&atmcrg_v[0]-1;
    chgpos=&chgpos_v[0]-1;
    crgatn=&crgatn_v[0]-1;
    atmeps=&atmeps_v[0]-1;

    //epsdim=iNatom+iNObject+2; replaced in space module semi global

    /**
     * find charge moments for dipole approximation
     */
    qnet=0.0;
    qplus=0.0;
    qmin=0.0;
    cqplus.nX=0.; cqplus.nY=0.; cqplus.nZ=0.; //cqplus={0.,0.,0.};
    cqmin.nX=0.;  cqmin.nY=0.;  cqmin.nZ=0.;  //cqmin={0.,0.,0.};

    ic1=0;
    for(ix=1; ix<=iNatom; ix++)
    {
        // (Parallel) check if atom is in the net grid
        if(checkNotBufferCoor(sDelPhiPDB[ix].xyz) || ibctyp == 3)
        {
            if (abs(sDelPhiPDB[ix].charge)>1.e-6)
            {
                ic1=ic1+1;
                chrg=sDelPhiPDB[ix].charge;
                atmcrg[ic1].nGrid=xn2[ix];
                atmcrg[ic1].nValue=chrg;
                chgpos[ic1]=xn1[ix];

                //modified for Parallel
                crgatn[ic1]= ix;
                qnet=qnet + chrg;

                if (chrg>0.)
                {
                    qplus=qplus + chrg;

                    // 2011-05-30 Using operations on coord type variables defined in
                    // module operators_on_coordinates
                    cqplus=cqplus+(chrg*atmcrg[ic1].nGrid);
                }
                else
                {
                    qmin=qmin + chrg;
                    cqmin=cqmin+(chrg*atmcrg[ic1].nGrid);
                }

                xyz[0] = (xn1[ix].nX - fgBoxCenter.nX);
                xyz[1] = (xn1[ix].nY - fgBoxCenter.nY);
                xyz[2] = (xn1[ix].nZ - fgBoxCenter.nZ);

                for (ridx = 0; ridx < 3; ridx++)
                    dMoment[ridx] += chrg*xyz[ridx];

                r = xn1[ix] - fgBoxCenter;
                r2 = optDot(r,r);

                for ( ridx = 0; ridx < 3; ridx++)
                {
                    for (cidx = 0; cidx < 3; cidx++)
                    {
                        if ( ridx == cidx)
                            delta = 1;
                        else
                            delta = 0;

                        qMoment[(3*ridx) + cidx] += chrg * ((3*xyz[ridx]*xyz[cidx] - r2*delta));
                    }
                }

            }
        }
    }

    #ifdef VERBOSE
    if(!(iGaussian==1) && inhomo==1 && logs)
        cout <<"number of charges coming from molecules " << ic1 << endl;

    if (ndistr>0) //insert charges from charge distributions
        cout << "Warning: call distrTOpoFloat2Int(ic1,false)..." << endl;
    #endif // VERBOSE

    //assign charges for boundary conditions: ic1 = number of charges

    if (    qplus > 1.e-6) cqplus = cqplus/qplus; //divide by charge totals
    if (abs(qmin) > 1.e-6) cqmin  = cqmin/qmin;

    // select those charges which will be charging the grid
    // Arrays of correct size are already allocated
    // rgrid=iGrid;
    ic2=0;

    // 2011-06-02 First determine correct size of arrays chgrv2 and nqgrdtongass
    for(ix=1; ix<=nqass; ix++)
    {
        // (Parallel) check if atom is in the net grid
        if ( ibctyp != 3)
        {
            if (checkNotBufferGrid(optCast<delphi_integer, delphi_real>(atmcrg[ix].nGrid)))
                ic2 = ic2 + 1;
        }
        else // ibctyp == 3
        {
            rgrid = iGrid.nX;
            if ( optANDGT(atmcrg[ix].nGrid,1.0) && optANDLT(atmcrg[ix].nGrid,rgrid) )
                ic2=ic2+1;
        }
    }

    nqgrd=ic2;
    ic2=0;

    chrgv2_v.assign(nqgrd,temp);
    nqgrdtonqass_v.assign(nqgrd,0);
    chrgv2=&chrgv2_v[0]-1;
    nqgrdtonqass=&nqgrdtonqass_v[0]-1;

    for(ix=1; ix<=nqass; ix++)
    {
        //crgatn[crg number]=atom number or iNatom+objectnumber or - distr.number
        ii=crgatn[ix];

        if (ii<0)
        {
            /*
             * now we have to consider charge distributions. In this case the distribution
             * is not linked to any object
             *    [jx,jy,jz]=coordinates of closest grid pointeger to charge
             *    (rx,ry,rz)=coordinates of the charge relatives to the current grid point
             */

            // 2011-06-02 Using operations on coord and int_coord type
            // variables defined in module operators_on_coordinates
            jxyz = optCast<delphi_integer,delphi_real>(atmcrg[ix].nGrid+0.5);
            rxyz = atmcrg[ix].nGrid -optCast<delphi_real,delphi_integer>(jxyz);
            jx   = jxyz.nX;
            jy   = jxyz.nY;
            jz   = jxyz.nZ;

            if (rxyz.nZ>rxyz.nX)
            {
                if (rxyz.nZ>-rxyz.nX)
                {
                    if (rxyz.nZ>rxyz.nY)
                    {
                        if (rxyz.nZ>-rxyz.nY)
                            imed=iepsmp[jz][jy][jx].nZ;
                        else
                            imed=iepsmp[jz][jy-1][jx].nY;
                    }
                    else
                        imed=iepsmp[jz][jy][jx].nY;
                }
                else
                {
                    if (rxyz.nY>rxyz.nX)
                    {
                        if (rxyz.nY>-rxyz.nX)
                            imed=iepsmp[jz][jy][jx].nY;
                        else
                            imed=iepsmp[jz][jy][jx-1].nX;
                    }
                    else
                        imed=iepsmp[jz][jy-1][jx].nY;
                }
            }
            else
            {
                if (rxyz.nZ>-rxyz.nX)
                {
                    if (rxyz.nY>rxyz.nX)
                        imed=iepsmp[jz][jy][jx].nY;
                    else
                    {
                        if (rxyz.nY>-rxyz.nX)
                            imed=iepsmp[jz][jy][jx].nX;
                        else
                            imed=iepsmp[jz][jy-1][jx].nY;
                    }
                }
                else
                {
                    if (rxyz.nZ>rxyz.nY)
                        imed=iepsmp[jz][jy-1][jx].nY;
                    else
                    {
                        if (rxyz.nZ>-rxyz.nY)
                            imed=iepsmp[jz][jy][jx].nY;
                        else
                            imed=iepsmp[jz-1][jy][jx].nZ;
                    }
                }
            }

            imed=imed/epsdim;
        }
        else
            imed=iAtomMed[ii];

        atmeps[ix]=medeps[imed];

        // (Parallel) check if atom is in the net grid
        if (ibctyp != 3)
        {
            if (checkNotBufferGrid(optCast<delphi_integer, delphi_real>(atmcrg[ix].nGrid)))
            {
                ic2=ic2+1;
                chrgv2[ic2].nGrid=atmcrg[ix].nGrid;
                chrgv2[ic2].nValue=atmcrg[ix].nValue;
                nqgrdtonqass[ic2]=ix;
            }
        }
        else // ibctyp == 3
        {
            if ( optANDGT(atmcrg[ix].nGrid,1.0) && optANDLT(atmcrg[ix].nGrid,rgrid) )
            {
                ic2=ic2+1;
                chrgv2[ic2].nGrid=atmcrg[ix].nGrid;
                chrgv2[ic2].nValue=atmcrg[ix].nValue;
                nqgrdtonqass[ic2]=ix;
            }// if
        }
    } // end of for(ix=1; ix<=nqass; ix++)

    ipassed=false;
    int zero_radius_cnt = 0;

    for(i=1; i<=nqass; i++)
    {
        ii=crgatn[i];

        if (ii>0&&ii<=iNatom)
        {
            if (sDelPhiPDB[ii].radius<=0.)
            {
                ipassed=true;
                zero_radius_cnt++;

                #ifdef VERBOSE
                cout << ii << " " << delphipdb[ii-1].getAtInf() << " is charged but its radius is absurd. Chaging radius to 0. " << endl;
                #endif

                sDelPhiPDB[ii].radius=0;  // leave it at zero instead of changing it to radpolext
            }
        }
    }

    if(ipassed) {CZeroChargeRadius war(zero_radius_cnt);};

    /*
     * argo -- surfpot module
     * Charge distribution moments
     */
    if ( zetaOn )
    {
        cout << infoString << left << setw(MAXWIDTH)  << " Charge Monopole "      << " : " << qnet << endl;
        cout << infoString << left << setw(MAXWIDTH)  << " Charge Dipole Moment " << " : " << "[ " << dMoment[0] << " " << dMoment[1] << " " << dMoment[2] << " ]" << endl;
        cout << infoString << left << setw(MAXWIDTH)  << " Charge Quadrupole Moment " << " : " << endl ;
        for ( ridx = 0; ridx < 3; ridx++)
        {
            cout << infoString << left << "    " << "|    " ;
            for (cidx = 0; cidx < 3; cidx++)
                cout << qMoment[(3*ridx) + cidx] << "    ";
            cout << "    |" << endl;
        }
    }

}
