//ARGO: 2016-FEB-09 --> Including headers for I/O for VDW_ssurface grids
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/assign/std/vector.hpp>

#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "../space/space.h"

#ifdef PARALLEL_OMP
#include <omp.h>
#endif

using namespace std;
using namespace boost::assign;

/*-------------------------ZETA ZETA ZETA ----------------------------
 A whole lot of changes have been made for the zetaSurface purposes.
 Look for if (zetaON == 1 ) conditions and statements when you
 need to put up a different file for zeta-calculation.
---------------------------------------------------------------------*/
void CDelphiSpace::setout()
{
    SGrid <delphi_real>  sqtemp[31],rad2aavtemp[31];
    SGrid <delphi_real> * sq=sqtemp+15;
    SGrid <delphi_real> * rad2aav=rad2aavtemp+15;
    SGrid <delphi_integer>* ioff;   ioff = NULL;

    bool itobig,itest2;
    string strtmp,strtmp1;

    //2011-05-12 Non-standard float variable, thus necessary
    //delphi_real modul,modul2, mod2,modx,mody;

    //2011-05-12 Non-standard type of variable, thus necessary
    //SGrid <delphi_real> xa,xb,xc,xd,xp;
    SGrid <delphi_real> ddist,dxyz,ddxyz,xn;
    //SGrid <delphi_real> tmpvect1,tmpvect2;
    SGrid <delphi_real> rad2av,fxn,vtemp;
    SGrid <delphi_integer> ismin,ismax,idist,ixyz,itest,ixn,i123;

    /*
     * here fRadPrb is not zero only if one wants to map the ext}//ed
     * surf. that has to be done with fRadPrb(1) if solute-solvent
     * interface is concerned imedia = medium number for a object
     * a non-zero entry in iEpsMap indicates an atom # plus 1 (to
     * properly treat re-entrant mid-points later (15 Aug 93)
     */
    delphi_integer iac,ibox,igrdc,i,iv,ix,iy,iz;
    delphi_integer limmax,lim;
    delphi_real dis2min2,dis2min1,distsq,dist,rad,rad2;
    delphi_real rad2a,radmax2,radp2,radtest,radprobe;
    delphi_real radp,temp;

    //argo : atom residue name
    string atom_residue;
    bool isProtein = false;
    vector<string> AA;
    AA += "ALA","ARG","ASN","ASP","ASX","CYS","GLU","GLN","GLX","GLY","HIS","HSE","HSD","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL";

    delphi_real extendedRad,extendedRad2;

    #ifdef VERBOSE
    cout << " Creating Van der Waals Epsilon Map: " << endl;
    #endif

    radprobe=0; //this radprobe seems not useful.

    //Parallel, we change this to:
    //epsdim=iNatom+iNObject+2; replaced in space module semi global

    radmax2=0.0;
    itest2=false;

    //ARGO 15-FEB,2016 -> Taken directly from site_writeSite_cube.cpp
    delphi_real coeff = 0.5291772108, stepsize = 1.0/fScale;
    SGrid<delphi_real> origin = (fgBoxCenter-stepsize*optCast<delphi_real,delphi_integer>(iGrid- delphi_integer(1))/2.0)/coeff;

    if (iNatom>0)
    {
        for(ix=1; ix<=iNatom; ix++)
            radmax2=max(radmax2,sDelPhiPDB[ix].radius);

        /*
         * this is probably the best way to do it,depending on which
         * surf. is desired
         */
        temp=max(radprobe,fExternRadius);

        //ARGO
        if ( zetaOn == 1 ) temp=max(temp,zetaDistance);

        #ifdef VERBOSE
        if ( temp == zetaDistance && zetaOn == 1)
        {
            cout << " The box size of individual atoms increased to Surface Distance" << endl;
            cout << " ... Expect Delays ... " << endl;
        }
        #endif

        radmax2=fScale*(radmax2+temp);

        lim=1+radmax2;
        limmax = 12;
        itobig=false;

        if(lim>limmax) itobig=true;
        igrdc=pow((2*lim+1),3);
        ioff = new SGrid <delphi_integer> [igrdc];

        if (!itobig)
        {
            radtest= pow( (radmax2 + 0.5*sqrt(3.0)),2 );
            ibox=-1;

            //2011-05-12 Strange statement. May allocate or may not
            //allocate array that used later in the program
            //irrespectively of itobig value, thus moved array
            //allocation before if condition

            #ifdef PARALLEL_OMP
            #pragma omp for schedule(auto)
            #endif

            for(ix=-lim; ix<=lim; ix++)
            {
                for(iy=-lim; iy<=lim; iy++)
                {
                    for(iz=-lim; iz<=lim; iz++)
                    {
                        //2011-05-13 Replaced by faster operation
                        //2011-05-13 Using operations on coord and
                        //int_coord type variables defined in module
                        //operators_on_coordinates
                        idist=int_coord(ix,iy,iz);
                        dist=float ( optDot(idist,idist) );
                        ddist = dist + float(0.25) + optCast <delphi_real,delphi_integer> (idist);

                        if ((dist<radtest)|| optORLT( ddist,radtest ))
                        {
                            ibox++;
                            ioff[ibox]=idist;
                        }
                    }
                }
            }
        }
    }

    /*
     * set interiors in MOLECULES
     */
    //ARGO: 2016-FEB-09: modified lines here
    //ORIGINIAL: 	if(itest2||itobig) cout <<"setout method 1 " << itest2 << " " << itobig << endl;

    #ifdef VERBOSE
    if(itest2||itobig) cout <<" setout method 1 : itest2 = " << itest2 << " and itobig = " << itobig << endl;
    #endif

    //DoATOMS:
    for( iv=1; iv<=iNatom; iv++)
    {
        /**
         * restore values
         */
        rad= sDelPhiPDB[iv].radius;

        //ARGO: also fetch the resname of the atoms and check if its a protein residue
        atom_residue = sDelPhiPDB[iv].atom_resname;
        isProtein = std::find(AA.begin(), AA.end(), atom_residue) != AA.end();
        isProtein = true; //for the present version 7.0 | Will remove it later

        xn=xn2[iv];

        //2011-05-13 Removed GOTO statement
        if (rad<1.e-6)
            continue;

        //fScale radius to grid
        rad=rad*fScale;
        //ARGO:
        extendedRad = rad + zetaDistance*fScale; //ARGO -> *fscale with zetaDistance to go with other radial distance parameters (e.g. rad,radp,etc.)

        radp=rad+fExternRadius*fScale;
        rad=rad+radprobe*fScale;

        //ARGO
        extendedRad2 = extendedRad*extendedRad;
        rad2=rad*rad;
        radp2=radp*radp;

        /*
         * set dielectric map
         * check if sphere sits within limits of box
         */
        itest2=false;

        ismin=optCast <delphi_integer,delphi_real> (xn-radmax2-1.0);
        ismax=optCast <delphi_integer,delphi_real> (xn+radmax2+1.0);
        itest=ismin;
        ismin=optMin(ismin,iGrid);
        ismin=optMax(ismin, delphi_integer(1));
        if(itest!=ismin) itest2=true;
        itest=ismax;
        ismax=optMin(ismax,iGrid);
        ismax=optMax(ismax, delphi_integer(1));
        if(itest != ismax) itest2=true;

        if (itest2||itobig)   //slow method;
        {
            rad2a = rad2 - 0.25;

            for(iz=ismin.nZ; iz<=ismax.nZ; iz++)
            {
                for(iy=ismin.nY; iy<=ismax.nY; iy++)
                {
                    for(ix=ismin.nX; ix<=ismax.nX; ix++)
                    {
                        ixyz=int_coord(ix,iy,iz);
                        dxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
                        distsq=optDot(dxyz,dxyz);
                        dxyz=dxyz+distsq;

                        if (dxyz.nX<rad2a)
                            iepsmp[iz][iy][ix].nX=iv+1+iAtomMed[iv]*epsdim;

                        if (dxyz.nY<rad2a)
                            iepsmp[iz][iy][ix].nY=iv+1+iAtomMed[iv]*epsdim;

                        if (dxyz.nZ<rad2a)
                            iepsmp[iz][iy][ix].nZ=iv+1+iAtomMed[iv]*epsdim;

                        if(distsq<radp2) idebmap[iz][iy][ix] =false;

                        //ARGO: 2016-FEB-09 : if zeta is 'on' then check for that too
                        if (zetaOn == 1 && isProtein)
                            if ( distsq<extendedRad2 ) zetaSurfMap[iz][iy][ix] =false;
                    }
                }
            }
        }
        else  /**faster method;*/
        {
            //IT HAS PROBLEMS!!!! Walter (be careful before using also in multidilectric case!!!&&!isitmd
            //cout << "####faster method:" << endl;
            rad2a=rad2-0.25;

            ixn=optRound(xn);

            fxn=optCast <delphi_real,delphi_integer> (ixn)-xn;
            rad2av=rad2a-fxn;

            for(ix=-lim; ix<=lim; ix++)
            {
                vtemp= double(ix)+fxn;
                sqtemp[ix+15]=vtemp*vtemp;
                rad2aavtemp[ix+15]=rad2a-vtemp;

            }

            if (iNMedia>1&&bOnlyMol)
            {
                for(i=0; i<=ibox; i++)
                {
                    i123=ioff[i];
                    ixyz=ixn+i123;
                    ix=ixyz.nX;
                    iy=ixyz.nY;
                    iz=ixyz.nZ;
                    distsq = sqtemp[i123.nX+15].nX +sqtemp[i123.nY+15].nY + sqtemp[i123.nZ+15].nZ;

                    if (distsq<rad2aavtemp[i123.nX+15].nX)
                    {
                        iac=(iepsmp[ix][iy][iz].nX % epsdim)-1;

                        if (iac==-1||iac>iNatom)
                        {
                            iepsmp[ix][iy][iz].nX=iv+1+iAtomMed[iv]*epsdim;
                        }
                        else
                        {
                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
                            ddxyz.nX= ddxyz.nX+0.5;
                            dis2min1=optDot(ddxyz,ddxyz)-rad2;

                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn2[iac];
                            ddxyz.nX= ddxyz.nX+0.5;
                            dis2min2=optDot(ddxyz,ddxyz)-pow((sDelPhiPDB[iac].radius*fScale),2);

                            if (dis2min2>dis2min1) iac=iv;
                            iepsmp[ix][iy][iz].nX=iac+1+iAtomMed[iac]*epsdim;
                        }
                    }

                    if (distsq<rad2aavtemp[i123.nY+15].nY)
                    {
                        iac=(iepsmp[ix][iy][iz].nY % epsdim)-1;
                        if (iac==-1||iac>iNatom)
                        {
                            iepsmp[ix][iy][iz].nY=iv+1+iAtomMed[iv]*epsdim;
                        }
                        else
                        {
                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
                            ddxyz.nY= ddxyz.nY+0.5;
                            dis2min1=optDot(ddxyz,ddxyz)-rad2;

                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn2[iac];
                            ddxyz.nY= ddxyz.nY+0.5;
                            dis2min2=optDot(ddxyz,ddxyz)-pow((sDelPhiPDB[iac].radius*fScale),2);

                            if (dis2min2>dis2min1) iac=iv;
                            iepsmp[ix][iy][iz].nY=iac+1+iAtomMed[iac]*epsdim;
                        }
                    }

                    if (distsq<rad2aavtemp[i123.nZ+15].nZ)
                    {
                        iac=(iepsmp[ix][iy][iz].nZ%epsdim)-1;
                        if (iac==-1||iac>iNatom)
                        {
                            iepsmp[ix][iy][iz].nZ=iv+1+iAtomMed[iv]*epsdim;
                        }
                        else
                        {

                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
                            ddxyz.nZ= ddxyz.nZ+0.5;
                            dis2min1=optDot(ddxyz,ddxyz)-rad2;

                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn2[iac];
                            ddxyz.nZ=ddxyz.nZ+0.5;
                            dis2min2=optDot(ddxyz,ddxyz)-pow((sDelPhiPDB[iac].radius*fScale),2);

                            if (dis2min2>dis2min1) iac=iv;
                            iepsmp[ix][iy][iz].nZ=iac+1+iAtomMed[iac]*epsdim;
                        }
                    }

                    if(distsq<radp2) idebmap[ix][iy][iz]=false;

                    //ARGO: 2016-FEB-09 : if zeta is 'on' then check for that too
                    if (zetaOn == 1 && isProtein)
                        if ( distsq<extendedRad2 ) zetaSurfMap[ix][iy][iz] =false;
                }
            }
            else
            {

                #ifdef PARALLEL_OMP
                #pragma omp for schedule(auto)
                #endif

                for(i=0; i<=ibox; i++)
                {
                    i123=ioff[i];
                    ixyz=ixn+i123;
                    ix=ixyz.nX;
                    iy=ixyz.nY;
                    iz=ixyz.nZ;
                    distsq = sq[i123.nX].nX +sq[i123.nY].nY + sq[i123.nZ].nZ;

                    if (distsq<rad2aav[i123.nX].nX)
                        //Lin Li: all indexes -1, because of C++ arrays' indexes start from 0;
                        //Lin Li: because iv start from 0, so +2:
                        iepsmp[iz][iy][ix].nX=iv+1+iAtomMed[iv]*epsdim;

                    if (distsq<rad2aav[i123.nY].nY)
                        iepsmp[iz][iy][ix].nY=iv+1+iAtomMed[iv]*epsdim;

                    if (distsq<rad2aav[i123.nZ].nZ)
                        iepsmp[iz][iy][ix].nZ=iv+1+iAtomMed[iv]*epsdim;

                    if (distsq<radp2) idebmap[iz][iy][ix]=false;

                    //ARGO: 2016-FEB-09 : if zeta is 'on' then check for that too
                    if (zetaOn == 1 && isProtein)
                        if ( distsq<extendedRad2 )
                            zetaSurfMap[iz][iy][ix]=false;
                }
            }
        }
    }

    //ARGO
    if (ioff != NULL) delete [] ioff;

    #ifdef VERBOSE
    cout << " Not Testing iepsmp and idebmap xxx " << endl;
    cout <<" Ending creating Van der Waals Epsilon Map " << endl;
    #endif

    return;

}// void setout;
