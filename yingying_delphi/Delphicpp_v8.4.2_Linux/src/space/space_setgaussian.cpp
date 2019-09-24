/**
 * Tips for Gaussian:
 * the new statement can be added by modifying following files:
 *
 * delphi_constants.h : increase iStatementNum
 * delphi_datamarshal.h
 * delphi_data_setMap
 * delphi_data_setDefault
 * delphi_data_getstatement
 *
 * new variables and arrays increased in Gaussian & MEMPOT:
 * cutoff,sigma,srfcut,radipz,inhomo,gaussian,ergsgaussian
 * gepsmp(:,:,:),gepsmp2(:,:,:)
 */

#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"
#include "space.h"

#ifdef PARALLEL_OMP
#include <omp.h>
#endif

using namespace std;

void CDelphiSpace::setGaussian()
{
    SGrid <delphi_real>  sqtemp[301],rad2aavtemp[301],vtemp2temp[301];
    delphi_real eps_min,eps_max,eps_diff,bnd_buff;
    SGrid <delphi_real> * sq=sqtemp+150;
    SGrid <delphi_real> * rad2aav=rad2aavtemp+150;
    SGrid <delphi_real> * vtemp2=vtemp2temp+150;
    SGrid <delphi_integer>* ioff;

    ioff=NULL;
    bool itobig,itest2,ipore,bOnlyMolbDebug;
    string strtmp,strtmp1;

    // 2011-05-12 Non-standard int variable, thus necessary
    //int epsdim, objecttype,iflag;
    //delphi_integer epsdim;

    // 2011-05-12 Non-standard float variable, thus necessary
    float modul,modul2, mod2,modx,mody,dentemp;

    // 2011-05-12 Non-standard type of variable, thus necessary
    SGrid <delphi_real> xa,xb,xc,xd,xp,ddist,xyz2,dxyz,ddxyz,xn,vecdist;
    SGrid <delphi_real> tmpvect1,tmpvect2,origin;
    SGrid <delphi_real> vectx,vecty,vectz,rad2av,fxn,vtemp;
    SGrid <delphi_integer> ismin,ismax,idist,idist1,ixyz,itest,ixn,i123;
    delphi_real coeff,stepsize;
    delphi_integer longint;

    int iboxt,iac,ibox,ii,igrdc,i,j,k,imedia,iv,ix,iy,iz,kind;
    delphi_integer limmax,lim,n;
    //integer,dimension(1:6) ::inwater //mid point;
    //logical,dimension(1:6) ::ifinwater //mid point;
    delphi_real alpha,dis2min2,dis2min1,dot,distsq,exrd,dist,rad,rad2;
    delphi_real rad2a,rad4,radius,radmax2,rad5,radp2,radtest,radprobe;
    delphi_real radp,tan2,temp,tmp,tmp1,tmp2,epstemp;
    delphi_real peak,pi,den,distance2,sigmatime,radsq; //sigma=sigma*radius.
                                                       //peak=4/(3*pi**0.5*sigma**3)

    eps_diff=1.05;
    //coeff=0.5291772108;
    //stepsize=1.0/fScale;
    //origin=(cOldMid-stepsize*(iGrid-1)/2)/coeff;

    if(!(iGaussian==1&&inhomo==0&&logs))
    {
        pi=3.14159265359;
        sigmatime=3.0;
        for(i=1; i<=iGrid.nX; i++)
        {
            for(j=1; j<=iGrid.nY; j++)
            {
                for(k=1; k<=iGrid.nZ; k++)
                {
                    gepsmp2[i][j][k].nX=repsout;
                    gepsmp2[i][j][k].nY=repsout;
                    gepsmp2[i][j][k].nZ=repsout;

                    gepsmp[i][j][k].nX=0.0;
                    gepsmp[i][j][k].nY=0.0;
                    gepsmp[i][j][k].nZ=0.0;

                    //Initialize Gaussian Salt density array to zeros
                    gDensityMapOnGridPoint[i][j][k] = 0.0;
                }
            }
        }

        // -------gepsmp file: --------
//        {
//            ofstream densfile;
//            densfile.open ("test_space_gaussian_gepsmp1.txt");
//
//            densfile << fixed << setprecision(7);
//            for(ix=1; ix<=iGrid.nX; ix++)
//                for(iy=1; iy<=iGrid.nY; iy++)
//                    for(iz=1; iz<=iGrid.nZ; iz++)
//                        densfile << setw(6) << right << ix << " "
//                                 << setw(6) << right << iy << " "
//                                 << setw(6) << right << iz << " "
//                                 << setw(8) << right << gepsmp[ix][iy][iz].nX << " "
//                                 << setw(8) << right << gepsmp[ix][iy][iz].nY << " "
//                                 << setw(8) << right << gepsmp[ix][iy][iz].nZ << endl;
//            densfile.close();
//        }

        #ifdef VERBOSE
        cout << "Starting creating Van der Waals Epsilon Map: " << endl;
        #endif

        radprobe=0; // this radprobe seems not useful.
        //epsdim=iNatom+iNObject+2; replaced in space module semi global
        //iboxt=0; //not useful
        radmax2=0.0;
        itest2=false;

        if (iNatom>0)
        {
            for(ix=1; ix<=iNatom; ix++)
            {
                //2011-05-13 Changed to derive-type array sDelPhiPDB (pointers module)
                radmax2=max(radmax2,sDelPhiPDB[ix].radius);

            }

            //this is probably the best way to do it,dep}//ing on which surf. is desired
            temp=max(radprobe,fExternRadius);

            //radmax2=fScale*(radmax2+temp);
            radmax2=sigmatime*fScale*(radmax2*sigma+temp); //Gaussian: 3 sigma plus temp. sigma=sigma*radius.Now radmax is 3 sigma.
            lim=1+radmax2;
            limmax = 1200; //Gaussian:original value:12

            itobig=false;
            if(lim>limmax) itobig=true;
            igrdc=pow((2*lim+1),3);
            ioff = new SGrid <delphi_integer> [igrdc];

            if (!itobig)
            {
                radtest= pow( (radmax2 + 0.5*sqrt(3.0)),2 );
                ibox=-1;

                /*
                 * 2011-05-12 Strange statement. May allocate or may not allocate
                 * array that used later in the program irrespectively of itobig value,
                 * thus moved array allocation before if condition
                 */

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
                            //int_coord type variables defined in module operators_on_coordinates
                            idist=int_coord(ix,iy,iz);
                            dist=float ( optDot(idist,idist) );
                            ddist = dist + float(0.25) + optCast <delphi_real,delphi_integer> (idist);

                            if ((dist<radtest)|| optORLT( ddist,radtest ))
                            {
                                ibox++;
                                ioff[ibox]=idist;
                            }
                        }// iz
                    }// iy
                }// ix
            }// itobig
        }// iNatom > 0

        //set interiors in OBJECTS

        /*
         * set interiors in MOLECULES
         */
        #ifdef VERBOSE
        if(itest2||itobig) cout <<"setout method 1 " << itest2 << " " << itobig << endl;
        #endif

        //DoATOMS:
        for( iv=1; iv<=iNatom; iv++)
        {
            /*
             * restore values
             */
            rad= sDelPhiPDB[iv].radius;
            rad= sigmatime*sigma*rad; // 3 sigma for Gaussian

            xn=xn2[iv];

            if (rad<1.e-6)
                continue;

            rad=rad*fScale; //fScale radius to grid
            radp=rad+fExternRadius*fScale; //rad5=pow( (rad+0.5),2);
            rad=rad+radprobe*fScale;
            rad2=rad*rad; //rad4=pow( (rad+0.5),2); // not used
            radp2=radp*radp;

            radsq=rad2/(sigmatime*sigmatime); // for Gaussian

            /*
             * set dielectric map
             *check if sphere sits within limits of box
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

            //if (itest2||itobig)   //slow method;
            if (itobig)   //slow method;
            {
                //2011-05-13 Seems redundant statement
                rad2a = rad2 - 0.25;

                #ifdef PARALLEL_OMP //Lin:2017.04.18
                #pragma omp for schedule(auto)
                #endif

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

                            //Gaussian Salt for slow method is under developing....

                        }
                    }
                }
            }
            else  /**faster method;*/
            {
                //IT HAS PROBLEMS!!!! Walter (be careful before using also in multidilectric case!!!&&!isitmd
                rad2a=rad2-0.25;

                ixn=optRound(xn);

                fxn=optCast <delphi_real,delphi_integer> (ixn)-xn;
                rad2av=rad2a-fxn;

                for(ix=-lim; ix<=lim; ix++)
                {
                    vtemp= double(ix)+fxn;
                    sqtemp[ix+150]=vtemp*vtemp;
                    rad2aavtemp[ix+150]=rad2a-vtemp;
                    vtemp2[ix]=vtemp;
                }

                //adjust inter-atom, different epsilon bgps+++04/2004 Walter
                if (iNMedia>1&&bOnlyMol)
                {
                    for(i=0; i<=ibox; i++)
                    {
                        i123=ioff[i];
                        ixyz=ixn+i123;
                        ix=ixyz.nX;
                        iy=ixyz.nY;
                        iz=ixyz.nZ;
                        distsq = sqtemp[i123.nX+150].nX +sqtemp[i123.nY+150].nY + sqtemp[i123.nZ+150].nZ;
                        if (distsq<rad2aavtemp[i123.nX+150].nX)
                        {
                            iac=(iepsmp[ix][iy][iz].nX % epsdim)-1;

                            if (iac==-1||iac>iNatom)
                            {
                                iepsmp[ix][iy][iz].nX=iv+1+iAtomMed[iv]*epsdim;
                            }
                            else
                            {
                                //2011-05-14 Using operations on coord and int_coord type variables defined in module operators_on_coordinates
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

                        if (distsq<rad2aavtemp[i123.nY+150].nY)
                        {
                            iac=(iepsmp[ix][iy][iz].nY % epsdim)-1;
                            if (iac==-1||iac>iNatom)
                            {
                                iepsmp[ix][iy][iz].nY=iv+1+iAtomMed[iv]*epsdim;
                            }
                            else
                            {
                                //2011-05-14 Using operations on coord and int_coord type variables defined in module operators_on_coordinates
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

                        if (distsq<rad2aavtemp[i123.nZ+150].nZ)
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

                        //Gaussian salt concentration
                        if (distsq < radp2)
                        {
                            //there may or may not be salt inside the atoms, lol
                        }
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

                        if(ix<iGrid.nX&&ix>1 && iy<iGrid.nY&&iy>1 && iz<iGrid.nZ&&iz>1)
                            distsq = sq[i123.nX].nX +sq[i123.nY].nY + sq[i123.nZ].nZ;
                        else // for grid outside the box
                            continue;

                        //---------------Lin Li: key section for Gaussian:
                        if (distsq<rad2aav[i123.nX].nX)
                        {
                            distance2=distsq+0.25+vtemp2[i123.nX].nX;
                            den=exp(-(distance2/(sigma*sigma*radsq)));
                            gepsmp[ix][iy][iz].nX=1-(1-gepsmp[ix][iy][iz].nX)*(1-den);
                        }

                        if (distsq<rad2aav[i123.nY].nY)
                        {
                            distance2=distsq+0.25+vtemp2[i123.nY].nY;
                            den=exp(-(distance2/(sigma*sigma*radsq)));
                            gepsmp[ix][iy][iz].nY=1-(1-gepsmp[ix][iy][iz].nY)*(1-den);
                        }

                        if (distsq<rad2aav[i123.nZ].nZ)
                        {
                            distance2=distsq+0.25+vtemp2[i123.nZ].nZ;
                            den=exp(-(distance2/(sigma*sigma*radsq))); //Gaussian: make density at center is 1.
                            gepsmp[ix][iy][iz].nZ=1-(1-gepsmp[ix][iy][iz].nZ)*(1-den);

                        }

                        //Section for Gaussian salt concentration
                        //The concentrations are recorded for each grid point
                        if (distsq < radtest)
                        {
                            den = exp(-(distsq / (sigma*sigma*radsq)));
                            gDensityMapOnGridPoint[ix][iy][iz] = 1 - (1 - gDensityMapOnGridPoint[ix][iy][iz])*(1 - den);
                        }
                    }
                }
            }
        }

        // -------gepsmp file: --------
//        {
//            ofstream densfile;
//            densfile.open ("test_space_gaussian_gepsmp2.txt");
//
//            densfile << fixed << setprecision(7);
//            for(ix=1; ix<=iGrid.nX; ix++)
//                for(iy=1; iy<=iGrid.nY; iy++)
//                    for(iz=1; iz<=iGrid.nZ; iz++)
//                        densfile << setw(6) << right << ix << " "
//                                 << setw(6) << right << iy << " "
//                                 << setw(6) << right << iz << " "
//                                 << setw(8) << right << gepsmp[ix][iy][iz].nX << " "
//                                 << setw(8) << right << gepsmp[ix][iy][iz].nY << " "
//                                 << setw(8) << right << gepsmp[ix][iy][iz].nZ << endl;
//            densfile.close();
//        }

    } //if(!(iGaussian==1&&inhomo==0&&logs))1010 continue;

    delphi_real temp_fDebFct= fEpsOut / (fDebyeLength*fScale*fDebyeLength*fScale);

    //Calculate Gaussian eps and fDebFct used for Gaussian Salt

    #ifdef PARALLEL_OMP //Lin:2017.04.18
    #pragma omp for schedule(auto)
    #endif

    for (ix = 1; ix <= iGrid.nX; ix++)
    {
        for (iy = 1; iy <= iGrid.nY; iy++)
        {
            for (iz = 1; iz <= iGrid.nZ; iz++)
            {
                gepsmp2[ix][iy][iz].nX=gepsmp[ix][iy][iz].nX*repsin+(1-gepsmp[ix][iy][iz].nX)*repsout;
                gepsmp2[ix][iy][iz].nY=gepsmp[ix][iy][iz].nY*repsin+(1-gepsmp[ix][iy][iz].nY)*repsout;
                gepsmp2[ix][iy][iz].nZ=gepsmp[ix][iy][iz].nZ*repsin+(1-gepsmp[ix][iy][iz].nZ)*repsout;

                //################### for set epsout in protein larger than in water:##########
                if(gepsmp[ix][iy][iz].nX<0.02)
                    gepsmp2[ix][iy][iz].nX=80.0;

                if(gepsmp[ix][iy][iz].nY<0.02)
                    gepsmp2[ix][iy][iz].nY=80.0;

                if(gepsmp[ix][iy][iz].nZ<0.02)
                    gepsmp2[ix][iy][iz].nZ=80.0;

                //Gaussian salt concentration
                if (gDensityMapOnGridPoint[ix][iy][iz]<0.02)
                    //density cut off 0.02
                    gDensityMapOnGridPoint[ix][iy][iz] = 0;

                //Gaussian salt density + protein density = 1
                gDensityMapOnGridPoint[ix][iy][iz] = 1 - gDensityMapOnGridPoint[ix][iy][iz];

                //Ion concentration = density * bulk conc
                //gDensityMapOnGridPoint[ix][iy][iz] *= temp_fDebFct;

                //################### }// for this epsout in protein different than in water######
            }
        }
    }

    if(inhomo==1)  //reduce epsilon out side protein;
    {
        epstemp=srfcut;

        if(dencut > 1)
        {
            cout << "dencut is greater than 1." << endl;
            exit(0);
        }
        if( dencut > 0)
        {
            epstemp=dencut;
        }

        if(dencut <0 && epstemp<repsin)
        {
            cout << "srfcut is lower than epsin." << endl;
            exit(0);
        }

        //Lin Li 2016 April: 3rd epsilon. The cutoff can be based on either epsilon or density:
        // epsilon cut off: srfcut;   density cut off: dencut
        // density cut off has higher priority than epsilon cut off.
        if(dencut > 0) // density cut off:
        {
            #ifdef PARALLEL_OMP //Lin:2017.04.18
            #pragma omp for schedule(auto)
            #endif
            for(ix=1; ix<=iGrid.nX; ix++)
            {
                for(iy=1; iy<=iGrid.nY; iy++)
                {
                    for(iz=1; iz<=iGrid.nZ; iz++)
                    {
                        if(gepsmp[ix][iy][iz].nX<epstemp) gepsmp2[ix][iy][iz].nX=repsout2;
                        if(gepsmp[ix][iy][iz].nY<epstemp) gepsmp2[ix][iy][iz].nY=repsout2;
                        if(gepsmp[ix][iy][iz].nZ<epstemp) gepsmp2[ix][iy][iz].nZ=repsout2;
                    }
                }
            }
        }
        else // epsilon cut off:
        {
            #ifdef PARALLEL_OMP //Lin:2017.04.18
            #pragma omp for schedule(auto)
            #endif
            for (ix = 1; ix <= iGrid.nX; ix++)
            {
                for (iy = 1; iy <= iGrid.nY; iy++)
                {
                    for (iz = 1; iz <= iGrid.nZ; iz++)
                    {
                        if(gepsmp2[ix][iy][iz].nX>epstemp) gepsmp2[ix][iy][iz].nX=repsout2;
                        if(gepsmp2[ix][iy][iz].nY>epstemp) gepsmp2[ix][iy][iz].nY=repsout2;
                        if(gepsmp2[ix][iy][iz].nZ>epstemp) gepsmp2[ix][iy][iz].nZ=repsout2;
                    }
                }
            }
        }
    }

    iBoundNum=0;
    longint=0;
    dentemp=0.1;


    #ifdef PARALLEL_OMP //Lin:2017.04.18
    #pragma omp for schedule(auto)
    #endif

    for(i=2; i<iGrid.nX; i++)
    {
        for(j=2; j<iGrid.nY; j++)
        {
            for(k=2; k<iGrid.nZ; k++)
            {
                if(gepsmp[i][j][k].nX>dentemp||gepsmp[i][j][k].nY>dentemp ||gepsmp[i][j][k].nZ>dentemp||gepsmp[i-1][j][k].nX>dentemp ||gepsmp[i][j-1][k].nY>dentemp||gepsmp[i][j][k-1].nZ>dentemp)
                {
                    idebmap[k][j][i]=false; //iGaussian change method of generating bDebMap;

                    // see if the epslon distribution is flat or not
                    eps_min=gepsmp[i][j][k].nX;
                    eps_max=gepsmp[i][j][k].nX;

                    eps_min=min(eps_min,gepsmp[i][j][k].nY);
                    eps_min=min(eps_min,gepsmp[i][j][k].nZ);
                    eps_min=min(eps_min,gepsmp[i-1][j][k].nX);
                    eps_min=min(eps_min,gepsmp[i][j-1][k].nY);
                    eps_min=min(eps_min,gepsmp[i][j][k-1].nZ);

                    eps_max=max(eps_max,gepsmp[i][j][k].nY);
                    eps_max=max(eps_max,gepsmp[i][j][k].nZ);
                    eps_max=max(eps_max,gepsmp[i-1][j][k].nX);
                    eps_max=max(eps_max,gepsmp[i][j-1][k].nY);
                    eps_max=max(eps_max,gepsmp[i][j][k-1].nZ);

                    if(eps_max/eps_min > eps_diff) longint=longint+1;
                }
            }
        }
    }

    iBoundNum=longint;
    ibgrd_v.assign(iBoundNum, sgrid_temp_int);
    ibgrd=&ibgrd_v[-1];

    n=0;

    #ifdef PARALLEL_OMP //Lin:2017.04.18
    #pragma omp for schedule(auto)
    #endif

    for(i=2; i<iGrid.nX; i++)
    {
        for(j=2; j<iGrid.nY; j++)
        {
            for(k=2; k<iGrid.nZ; k++)
            {
                if( gepsmp[i][j][k].nX>dentemp || gepsmp[i][j][k].nY>dentemp ||
                    gepsmp[i][j][k].nZ>dentemp||gepsmp[i-1][j][k].nX>dentemp ||
                    gepsmp[i][j-1][k].nY>dentemp||gepsmp[i][j][k-1].nZ>dentemp )
                {
                    eps_min=gepsmp[i][j][k].nX;
                    eps_max=gepsmp[i][j][k].nX;

                    eps_min=min(eps_min,gepsmp[i][j][k].nY);
                    eps_min=min(eps_min,gepsmp[i][j][k].nZ);
                    eps_min=min(eps_min,gepsmp[i-1][j][k].nX);
                    eps_min=min(eps_min,gepsmp[i][j-1][k].nY);
                    eps_min=min(eps_min,gepsmp[i][j][k-1].nZ);

                    eps_max=max(eps_max,gepsmp[i][j][k].nY);
                    eps_max=max(eps_max,gepsmp[i][j][k].nZ);
                    eps_max=max(eps_max,gepsmp[i-1][j][k].nX);
                    eps_max=max(eps_max,gepsmp[i][j-1][k].nY);
                    eps_max=max(eps_max,gepsmp[i][j][k-1].nZ);

                    if(eps_max/eps_min > eps_diff)
                    {
                        n=n+1;
                        ibgrd[n].nX=i;
                        ibgrd[n].nY=j;
                        ibgrd[n].nZ=k;
                    }
                }
            }
        }
    }

    if (ioff != NULL)
    {
        delete [] ioff;
        ioff=NULL;
    }

    #ifdef VERBOSE
    cout <<"Ending creating Van der Waals Epsilon Map " << endl;
    #endif

    return;

}// void setiGaussian;

