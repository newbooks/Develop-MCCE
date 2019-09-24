#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "../space/space.h"
#include <fstream>; //ARGO 12-FEB,2016

/* ARGO: 2016-FEB-09
 * Snippet of code meant to find the find boudary points on the zeta-surface Will use
 * the zetaSurfMap-3D array of booleans to gather the boundary points along the lines
 * of how Delphi already does it for VDW surface. --> Couple of new varibales have been defined
 */

using namespace std;

void CDelphiSpace::VdwToMs()
{
    string line;
    delphi_integer nbra[SPACE_NBRA_SIZE];
    delphi_integer eps[7],nt;
    delphi_integer itmp[7],dim1,cont;

    //ARGO
    bool zeta_tmp[7];
    int zeta = 1;
    delphi_integer zbgp,zext;
    delphi_integer lx,ly,lz;
    delphi_integer px,py,pz; // for reading the points at which BND is detected from zetaSurfMap
    delphi_real stepsize = 1.0/fScale, loc_fRMid=(iGrid.nX+1)/2.0;
    SGrid<delphi_real> cornerOrigin = (fgBoxCenter-stepsize*((iGrid-1)/2.0)), xyz;
    delphi_integer zleft = 0, zright = 0, ztop = 0, zbottom = 0, zfront = 0, zback = 0;
    SGrid <delphi_integer> lxyz;

    delphi_integer iaprec,dim,isign;
    delphi_integer imap[5][7]= {0};
    delphi_integer kind,eps2[7];
    bool remov;
    string strtmp;
    SGrid <delphi_real> xq,dxyz,dr123,dx123,u123;
    SGrid <delphi_real> goff[7]= {0.,0.,0.},xg,x1xyz,s123,xxyyzz;
    SGrid <delphi_integer> ixyz,it123,jxyz;

    delphi_real coordx,coordy,coordz; //for int_coord operations on zetaSurface

    bool exists,flag,cycle_flag=false;
    bool nbe[7]= {false};
    delphi_real cbln,cba,del,dis,dmn,dist=0,ds2,off,dsr,r0a;
    SGrid<delphi_real> offf;
    delphi_real x1,prbrd12,prbrd22,rsm,prbrd2;
    delphi_integer ia,i,iac,ibgp,ii,iext,iarv,iacv,imedia;
    delphi_integer iiord,ix2,iy2,iz2,ix,iy,iz;
    delphi_integer iii,iord,j,jj,jjj,jx,jy,jz,limu,liml,kk,m,k,n;
    delphi_integer mr,mpr,n2,ndv,ncav,nnn,nnbr;
    delphi_integer nn,ntt,n1,Nar;
    delphi_integer **** bndeps;  //bndeps is a local variable in C++
    delphi_integer ibmx=1000000; //ibmx is a local variable in C++

    if(optORGT(iGrid, delphi_integer(300))) ibmx=50000000;

    #ifdef VERBOSE
    cout << " Info> Drawing MS from vdW surface" << endl;
    #endif

    iarv=0; // Lin Li
    SGrid <delphi_integer> * ibnd = new SGrid <delphi_integer> [ibmx];

    r0.assign(iNatom+1,0.0);
    r02.assign(iNatom+1,0.0);
    rs2.assign(iNatom+1,0.0);
    ast.assign(iNatom+1,0);

    get_pt4d <delphi_integer> (bndeps,iGrid.nX+1,iGrid.nY+1,iGrid.nZ+1,3);

    offf=(iGrid+ delphi_integer(1))/2.;

    //epsdim=iNatom+iNObject+2; replaced in space module semi global

    /*
     * imap maps from midpoint position to iepsmap entry positions
     */
    imap[1][4]=-1;
    imap[2][5]=-1;
    imap[3][6]=-1;
    imap[4][1]=1;
    imap[4][2]=2;
    imap[4][3]=3;
    imap[4][4]=1;
    imap[4][5]=2;
    imap[4][6]=3;

    nbe[1]=true;
    nbe[2]=true;
    nbe[3]=true;
    nbe[4]=true;
    nbe[5]=true;

    off=0.5/fScale;

    goff[1].nX=off;
    goff[2].nY=off;
    goff[3].nZ=off;
    goff[4].nX=-off;
    goff[5].nY=-off;
    goff[6].nZ=-off;

    radpmax=max(fRadPrb[1],fRadPrb[2]); //transferred via qlog module;

    // conversion from grid to float coordinates(can also use routine gtoc)
    // 2011-05-17 Using operations on coord and int_coord type
    // variables defined in module operators_on_coordinates
    x1=1.0/fScale;
    if (ibctyp == 3)
        x1xyz = cOldMid - double( 0.5*x1*(iGrid.nX + 1) );
    else
        x1xyz = myStartCoor - 0.5*x1;

    // find extrema
    cMin.nX=6000.;
    cMin.nY=6000.;
    cMin.nZ=6000.;

    cMax.nX=-6000.;
    cMax.nY=-6000.;
    cMax.nZ=-6000.;

    for(ii=0; ii<=iNObject-1; ii++)
    {
        cMin=optMin(cMin,sLimObject[ii].nMin);
        cMax=optMax(cMax,sLimObject[ii].nMax);

    }

    // find vanderwaals boundary
    n=0;
    nn=0;

    // NB change limits to those of the molecule. set for iepsmp NOT equal to unity

    // Lin Li: determine cube 2 grids smaller than limeps;
    for(k=LimEps.nMin.nZ+1; k<=LimEps.nMax.nZ-1; k++)
    {
        for(j=LimEps.nMin.nY+1; j<=LimEps.nMax.nY-1; j++)
        {
            for(i=LimEps.nMin.nX+1; i<=LimEps.nMax.nX-1; i++)
            {
                // one distinguishes between internal,external,internal bgp and external bgp
                iext=0;
                ibgp=0;

                //Lin Li: six neighbors: itmp is iatmmed or 0;
                itmp[1]=abs(iepsmp[k][j][i].nX)/epsdim;   //right
                itmp[2]=abs(iepsmp[k][j][i].nY)/epsdim;   //front
                itmp[3]=abs(iepsmp[k][j][i].nZ)/epsdim;   //up
                itmp[4]=abs(iepsmp[k][j][i-1].nX)/epsdim; //left
                itmp[5]=abs(iepsmp[k][j-1][i].nY)/epsdim; //back
                itmp[6]=abs(iepsmp[k-1][j][i].nZ)/epsdim; //down

                if(itmp[1]==0) iext=1; //external point
                if(itmp[1]!=itmp[6]) ibgp=1; //bnd point

                for(cont=2; cont<=6; cont++)
                {
                    if(itmp[cont]==0) iext=1; //external point
                    if(itmp[cont]!=itmp[cont-1]) ibgp=1; //bnd point
                }

                /*
                 * assignment of right values to bndeps according to the point nature from now
                 * iBoundNum is the total number of internal and external boundary grid points
                 */
                if (ibgp>0)
                {
                    n=n+1; //boundary bnd point
                    bndeps[i][j][k][1]=n;
                    bndeps[i][j][k][2]=iext;

                    if (iext>0) nn=nn+1; //bnd point and external point: boundary + surface
                    ibnd[n]=int_coord(i,j,k);
                }
                else
                {
                    bndeps[i][j][k][1]=0;
                    bndeps[i][j][k][2]=0;
                }
            }
        }
    }

    /* ARGO:
     * Finding the boundary points on the zeta-surface
     * In this case, it writes a file "strZetaPhiFile + ".zTemp""
     * It resets zetaSurfMap to {true} and then assigns bnd at points
     * the file lists.
     * requires the file to be read else where
     */
    if (zetaOn == 1)
    {
        for (lz=2; lz < iGrid.nZ; lz++)
        {
            for (ly=2; ly < iGrid.nY; ly++)
            {
                for (lx=2; lx < iGrid.nX; lx++)
                {
                    zext=0;
                    zbgp=0;

                    zeta_tmp[0] = zetaSurfMap[lz][ly][lx];
                    zeta_tmp[1] = zetaSurfMap[lz+1][ly][lx]; // six neighbors: IN THE NEWS+/- FORMAT
                    zeta_tmp[2] = zetaSurfMap[lz-1][ly][lx];
                    zeta_tmp[3] = zetaSurfMap[lz][ly+1][lx];
                    zeta_tmp[4] = zetaSurfMap[lz][ly-1][lx];
                    zeta_tmp[5] = zetaSurfMap[lz][ly][lx+1];
                    zeta_tmp[6] = zetaSurfMap[lz][ly][lx-1];

                    zright=abs(iepsmp[lz][ly][lx].nX)/epsdim;    //right
                    zfront=abs(iepsmp[lz][ly][lx].nY)/epsdim;    //front
                    ztop=abs(iepsmp[lz][ly][lx].nZ)/epsdim;      //up
                    zleft=abs(iepsmp[lz][ly][lx-1].nX)/epsdim;   //left
                    zback=abs(iepsmp[lz][ly-1][lx].nY)/epsdim;   //back
                    zbottom=abs(iepsmp[lz-1][ly][lx].nZ)/epsdim; //down

                    //external point
                    if(zeta_tmp[0] && zright == 0 && zleft == 0 && ztop == 0 && zbottom == 0 && zfront == 0 && zback == 0)
                        zext=1;

                    for(cont=1; cont<=6; cont++)
                    {
                        // if(zeta_tmp[cont]) zext=1; //external point
                        if(zeta_tmp[cont]!=zeta_tmp[cont-1]) zbgp=1; //bnd point
                    }// do

                    if(zeta_tmp[6]!=zeta_tmp[1]) zbgp=1; //bnd point

                    if (zbgp>0 && zext==1)
                    {
                        lxyz = int_coord(lx, ly, lz);
                        xyz = fgBoxCenter + (optCast <delphi_real,delphi_integer> (lxyz) - loc_fRMid)/fScale;
                        coordz = xyz.nZ;
                        coordy = xyz.nY;
                        coordx = xyz.nX;

                        surf_grid_coords_v.push_back(coordx);
                        surf_grid_coords_v.push_back(coordy);
                        surf_grid_coords_v.push_back(coordz);

                        surf_grid_index_v.push_back(lx);
                        surf_grid_index_v.push_back(ly);
                        surf_grid_index_v.push_back(lz);
                    }
                }
            }
        }
    }//IF-ZETA IS TRUE

    iBoundNum=n;
    iBoundNumsurf=nn;
    nn=0;

    #ifdef VERBOSE
    cout <<" VdMS> Boundary points facing continuum solvent= " << iBoundNumsurf << endl;
    cout <<" VdMS> Total number of boundary points before elab.= " << iBoundNum << endl;
    #endif //VERBOSE

    if (iBoundNum>ibmx)
    {
        cout <<" WARNING !!! iBoundNum= " << iBoundNum << " is greater than ibmx = " << ibmx << endl;
        cout <<" WARNING !!! Increase ibmx in vwtms.f" << endl;
        exit(0);
    }

    //2011-05-17 Arrays allocated by ordinary F95 allocate statement

    if (radpmax<1.e-6) //prob radius = 0
    {
        ibgrd_v.assign(iBoundNum+1, sgrid_temp_int);
        ibgrd=&ibgrd_v[0]-1;

        //2011-05-17 Changed to array operations, but keeping some assignment in a cycle due to array size mismatch
        //ast=0; already initialized.
        for(i=1; i<=iBoundNum; i++)
            ibgrd[i]=ibnd[i];
    }
    else
    {
        if (!bOnlyMol && fRadPrb[1] != fRadPrb[2])
        {
            // !bOnlyMol: means object exist, so delete
        }
        else
        {
            for(i=1; i<=iNatom; i++)
            {
                r0a=sDelPhiPDB[i].radius+fRadPrb[1];
                r0[i]=r0a; //atom radius+ probe radius
                r02[i]=r0a*r0a;
                rs2[i]=0.99999*r02[i];
            }
        }

        /*
         * make a list of accessible points..,expos. all scaling of grid points will be done to
         * make these points..
         */
        prbrd12=fRadPrb[1]*fRadPrb[1]; //fRadPrb[1]: probe radius;fRadPrb[2]=0, no idea.
        prbrd22=fRadPrb[2]*fRadPrb[2];

        /*
         * calculate an identity for this conformation
         */
        rsm=0.0;
        if (bOnlyMol)
            for(i=1; i<=iNatom; i++)
                rsm=rsm+sDelPhiPDB[i].radius*optSum(optABS(xn1[i])); //Lin: is it a mistake? rsm? If rms is just for ARCDAT....

        flag=true;

        sas();

        del=1./fScale;
        del=max(del,radpmax);
        cbln=fRMax+del;

        cubedata(2.0,cbln);

        dim=(lcb+1)*(mcb+1)*(ncb+1);

        //2011-05-17 Array allocation is done by ordinary F95 ALLOCATE
        //statement Array allocated as 3D as in the cube void
        cbn1_v.assign(dim+1,0);
        cbn2_v.assign(dim+1,0);
        dim1=27;
        if ((iNObject-numbmol)>0) dim1=max(dim, delphi_integer(27));

        cbal.assign(dim1*(iNatom+iNObject-numbmol)+1,0);

        cube();

        //link the accessible points into iab1 and iab2
        indverdata(radpmax);

        cba=1./grdi;
        nnn=(lcb1+1)*(mcb1+1)*(ncb1+1);
        icume.assign(extot+1,0);

        #ifdef VERBOSE
        cout << " VdMS> grid for indexing accessible points =  " << cba << endl;
        #endif // VERBOSE
        indver(extot);

        //now start the expansion
        //m1= the number of boundary points removed from list
        ncav=0;
        n1=1;
        n2=iBoundNum;

        /*
         * m= number of new boundary elements..
         */
        mpr=100000;
        ndv=0;

        while(true)
        {
            m=0;
            mr=0;

            for(i=n1; i<=n2; i++)
            {
                ixyz=ibnd[i];
                ix=ixyz.nX;
                iy=ixyz.nY;
                iz=ixyz.nZ;

                //considering both internal and external b.g.p.
                if (bndeps[ix][iy][iz][1]!=0)
                {

                    /*
                     * still has to be considered what is external and what internal!!!!!WWW
                     * remov is true if it is an internal midpoint close to an interface where a molecule is present
                     * (expansion has to take place also in objects)
                     */
                    remov=false;

                    eps[1]=(iepsmp[iz][iy][ix].nX%epsdim);
                    eps[2]=(iepsmp[iz][iy][ix].nY%epsdim);
                    eps[3]=(iepsmp[iz][iy][ix].nZ%epsdim);
                    eps[4]=(iepsmp[iz][iy][ix-1].nX%epsdim);
                    eps[5]=(iepsmp[iz][iy-1][ix].nY%epsdim);
                    eps[6]=(iepsmp[iz-1][iy][ix].nZ%epsdim);

                    remov=((eps[1]>1&&eps[1]<=iNatom+1)||(eps[2]>1&&eps[2]<=iNatom+1));
                    remov=((eps[3]>1&&eps[3]<=iNatom+1)||(eps[4]>1&&eps[4]<=iNatom+1))||remov;
                    remov=((eps[5]>1&&eps[5]<=iNatom+1)||(eps[6]>1&&eps[6]<=iNatom+1))||remov;

                    eps2[1]=(iepsmp[iz][iy][ix].nX/epsdim);
                    eps2[2]=(iepsmp[iz][iy][ix].nY/epsdim);
                    eps2[3]=(iepsmp[iz][iy][ix].nZ/epsdim);
                    eps2[4]=(iepsmp[iz][iy][ix-1].nX/epsdim);
                    eps2[5]=(iepsmp[iz][iy-1][ix].nY/epsdim);
                    eps2[6]=(iepsmp[iz-1][iy][ix].nZ/epsdim);

                    /*
                     * cWWW there is still an issue in case there are both molecules and objects: since parent object of
                     * reentrant points is only known in sclbp, filling reentrant regions due to molecules in objects might fail
                     */
                    remov=remov&&(numbmol>0);

                    xg=(x1*optCast<delphi_real,delphi_integer>(ixyz))+x1xyz;

                    for (j=1; j<=6; j++)
                    {
                        if (eps[j]==0||(remov&&eps[j]>iNatom+1)||(eps2[j]==0&&eps[j]>0))
                        {
                            prbrd2=prbrd22;
                            if (eps[j]==0||eps2[j]==0) prbrd2=prbrd12;

                            //add midpoint offset to grid point..
                            s123=xg+goff[j];
                            //determine if this virgin midpoint is in or out
                            //2011-05-18 mn(x,y,z) and grdi were assigned values in INDVER void now coord type variable mnxyz is declared
                            //in pointers module and float grdi declared and thus accessible in qlog module
                            xxyyzz=(s123-mnxyz)*grdi;
                            jxyz=optCast<delphi_integer,delphi_real>(xxyyzz);
                            jx=jxyz.nX;
                            jy=jxyz.nY;
                            jz=jxyz.nZ;
                            //2011-05-18 Indexes lcb1, mcb1, ncb1 are transfred via qlog module and are set in INDVER void
                            if (optORLE(jxyz, delphi_integer(0))||optORGE(jxyz,lmncb1))
                            {
                                cout <<" VdMS> midpoint out of cube" << endl;
                                cout <<iepsmp[ iz][ iy ][ix ].nX << endl;
                                cout <<iepsmp[ iz][ iy ][ix ].nY << endl;
                                cout <<iepsmp[ iz][ iy ][ix ].nZ << endl;
                                cout <<iepsmp[ iz][ iy ][ix-1 ].nX << endl;
                                cout <<iepsmp[ iz][ iy-1 ][ix ].nY << endl;
                                cout <<iepsmp[ iz-1][ iy ][ix ].nZ << endl;
                            }// if
                            dmn=1000.;
                            iacv=0;

                            //2011-05-18 Repeating piece of the code now is in the separate file
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,0,0
                            jx=jx-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,0,0
                            jx=jx+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,-1,0
                            jx=jx-1;
                            jy=jy-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,1,0
                            jy=jy+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,0,-1
                            jy=jy-1;
                            jz=jz-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,0,1
                            jz=jz+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //nn=2
                            //1,0,1
                            jx=jx+1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,0,1
                            jx=jx-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,1,1
                            jx=jx+1;
                            jy=jy+1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,-1,1
                            jy=jy-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,-1,0
                            jz=jz-1;
                            jx=jx-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,-1,0
                            jx=jx+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,1,0
                            jy=jy+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,1,0
                            jx=jx-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,0,-1
                            jz=jz-1;
                            jy=jy-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,0,-1
                            jx=jx+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,1,-1
                            jx=jx-1;
                            jy=jy+1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,-1,-1
                            jy=jy-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //nn=3
                            //-1,-1,-1
                            jx=jx-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,-1,-1
                            jx=jx+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,1,-1
                            jy=jy+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,1,-1
                            jx=jx-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,1,1
                            jz=jz+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,1,1
                            jx=jx+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,-1,1
                            jy=jy-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,-1,1
                            jx=jx-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //it might be in the contact region findthe closest atom surface

                            it123=optCast<delphi_integer,delphi_real>((s123-xyzo)*cbai);
                            dmn=100.;
                            iac=0;
                            nnbr=0;
                            lmncb=int_coord(lcb,mcb,ncb);

                            if(optORLT(it123, delphi_integer(0))||optORGT(it123,lmncb))
                            {
                                //if the bgp is outside the cube,probably it is due to some object
                                for(ii=iNObject; ii<=1; ii--)
                                {
                                    strtmp=dataobject_v[(ii-1)*2];
                                    kind = atoi(strtmp.substr(15,3).c_str());
                                    if (strtmp.substr(0,4)!="is a" && kind!=2)
                                    {
                                        if ( optANDLE(s123,(sLimObject[ii].nMax+x1) ) && optANDGT(s123,(sLimObject[ii].nMin-x1)) )
                                        {
                                            nnbr=nnbr+1;
                                            nbra[nnbr]=ii+iNatom;
                                            liml=0;
                                            limu=0;
                                        }// if
                                    }// if
                                }// do

                                if(liml!=0||limu!=0)
                                {
                                    #ifdef VERBOSE
                                    cout <<" VdMS> a bgp close to nothing" << endl;
                                    #endif
                                }
                            }
                            else
                            {
                                //2011-05-19 Changed 1d array to 3d array as /in cube void
                                liml=cbn1[it123.nX][it123.nY][it123.nZ];
                                limu=cbn2[it123.nX][it123.nY][it123.nZ];
                            }

                            iaprec=0;

                            for( kk=liml; kk<=limu; kk++)
                            {
                                ia=cbal[kk];

                                #ifdef VERBOSE
                                if (ia==0) cout <<" VdMS> problems with cube" << endl;
                                #endif

                                if (ia<=iNatom&&ia>0)
                                {
                                    if (ast[ia]==0)
                                    {
                                        nnbr=nnbr+1;
                                        nbra[nnbr]=ia;
                                    }// if
                                }
                                else
                                {
                                    if (ia!=iaprec&&eps[j]==0)
                                    {
                                        //different from sclbp, I have to consider the object only if the
                                        //midpoint is external O SE C'E' UN PORO E VEDERE IN SEGUITO
                                        iaprec=ia;

                                        //assuming any object is not buried
                                        nnbr=nnbr+1;
                                        nbra[nnbr]=ia;
                                    }
                                }
                            }

                            for( ii=1; ii<=nnbr; ii++)
                            {
                                if(ii >= SPACE_NBRA_SIZE)
				cout << "space_vwtoms>> index beyond size of nbra: ii= "<< ii << endl;
				ia=nbra[ii];

                                if (ia>iNatom)
                                {
                                    iii=ia-iNatom;

                                    xq=s123;

                                    //an object can compete with an atom for a midpoint only if this midpoint is out of the object itself
                                    if (dist>=0.&&dist<dmn)
                                    {
                                        dmn=dist;
                                        iac=ia;
                                        dr123=(-dxyz)*(fRadPrb[1]-dist);
                                    }// if
                                }
                                else
                                {
                                    dx123=s123-xn1[ia];
                                    ds2=optDot(dx123,dx123);
                                    dis=sqrt(ds2)-sDelPhiPDB[ia].radius;

                                    //dis= distance to atom surface
                                    if (dis<dmn)
                                    {
                                        dmn=dis;
                                        iac=ia;
                                    }// if
                                }// if
                            }// do DOII;


                            if (iac==0)
                            {
                                ncav=ncav+1;
                                //possibly a cavity point
                            }
                            else
                            {
                                /**
                                 *check to see if it is in the contact region of that atom or object by
                                 *projecting it to the atom's acc surface and checking against the acc volumes of nearby atoms
                                 */
                                if (iac<=iNatom)
                                {
                                    dr123=s123-xn1[iac];
                                    dsr=sqrt(optDot(dr123,dr123));
                                    u123=xn1[iac]+((r0[iac]*dr123)/dsr);
                                }
                                else
                                {
                                    u123=s123-dr123;
                                }// if

                                it123=optCast<delphi_integer,delphi_real>((u123-xyzo)*cbai);

                                liml=cbn1[it123.nX][it123.nY][it123.nZ];
                                limu=cbn2[it123.nX][it123.nY][it123.nZ];

                                //DLIM:
                                for( kk=liml; kk<=limu; kk++)
                                {
                                    ia=cbal[kk];
                                    if (ia<=iNatom)
                                    {
                                        dx123=u123-xn1[ia];
                                        ds2=optDot(dx123,dx123);
                                        flag=true;

                                        if (ds2<rs2[ia])
                                        {
                                            flag=false;
                                            break;
                                        }// if
                                    }
                                    else
                                    {
                                        if (ia!=iac&&eps[j]==0)
                                        {
                                            xq=u123;

                                            flag=true;
                                            if (dist<-1.e-6)
                                            {
                                                flag=false;
                                                break;
                                            }

                                            //oriented distance from ext}//ed
                                            //object surface if negative => reentrant region
                                        }
                                    }
                                }

                                /*
                                 * it is in the contact region. flag the midpoint so it is not checked again iac is atom number...NOT increased by 1
                                 */

                                //2011-05-18 To get rid of goto 201
                                //statements above
                                if (flag)
                                {
                                    eps[j]=-iac;
                                    eps2[j]=-iac;
                                    continue;
                                }
                            }

                            eps[j]=1; //eps = 1 means cavity or reentrant;

                            //remap iepsmp
                            if (iac==0)
                            {
                                /**
                                 *this is an assumption, still to deeply understand meaning of cavity here and to improve this choice!WWW
                                 */
                                if (ia>0)
                                {
                                    imedia=iAtomMed[ia];
                                }
                                else
                                {
                                    #ifdef VERBOSE
                                    cout <<" VdMS> assigning arbitrary epsilon in cavity" << endl;
                                    #endif
                                    imedia=iAtomMed[1];
                                }
                            }
                            else
                            {
                                imedia=iAtomMed[iac];
                            }// if

                            switch (imap[4][j])
                            {
                                case 1:
                                    iepsmp[iz+imap[3][j]][iy+imap[2][j]][ix+imap[1][j]].nX=eps[j]+imedia*epsdim;
                                    //cout << "i,j,ix+imap[1][j],iy+imap[2][j],iz+imap[3][j],iepsmp: " << i << " " << j << " " << ix+imap[1][j] << " " <<iy+imap[2][j] << " " << iz+imap[3][j]<< " "  << iepsmp[iz+imap[3][j]][iy+imap[2][j]][ix+imap[1][j]] << endl;
                                    break;
                                case 2:
                                    iepsmp[iz+imap[3][j]][iy+imap[2][j]][ix+imap[1][j]].nY=eps[j]+imedia*epsdim;
                                    break;
                                case 3:
                                    iepsmp[iz+imap[3][j]][iy+imap[2][j]][ix+imap[1][j]].nZ=eps[j]+imedia*epsdim;
                                    break;
                                default:
                                    #ifdef VERBOSE
                                    cout <<" VdMS> ????? flag1" << endl;
                                    #endif
                                ;
                            }

                            eps2[j]=imedia;

                            /**
                             * not assigning the owner but only the medium, the former will be assigned in the fScale routine
                             * check to see if the nearest neighbour status has been changed..
                             */
                            ix2=ix;
                            iy2=iy;
                            iz2=iz;

                            /**
                             * if the nearest neighbour is a box boundary point { skip this since box boundary
                             * points can not also be dielctric boundary points
                             * 2011-05-18 Multiple IFs replaced by SELECT CASE
                             */
                            switch (j)
                            {
                                case 1:
                                    ix2=ix+1;
                                    if(ix2==iGrid.nX) continue;
                                    break;
                                case 2:
                                    iy2=iy+1;
                                    if(iy2==iGrid.nY) continue;
                                    break;
                                case 3:
                                    iz2=iz+1;
                                    if(iz2==iGrid.nZ) continue;
                                    break;
                                case 4:
                                    ix2=ix-1;
                                    if(ix2==1) continue;
                                    break;
                                case 5:
                                    iy2=iy-1;
                                    if(iy2==1) continue;
                                    break;
                                case 6:
                                    iz2=iz-1;
                                    if(iz2==1) continue;
                            }

                            /*
                             * once again one distinguishes between internal,external,internal bgp and external bgp
                             */
                            iext=0;
                            ibgp=0;

                            //2011-05-18 Changed to i nt_coord derived type
                            itmp[1]=abs(iepsmp[iz2][iy2][ix2].nX)/epsdim;
                            itmp[2]=abs(iepsmp[iz2][iy2][ix2].nY)/epsdim;
                            itmp[3]=abs(iepsmp[iz2][iy2][ix2].nZ)/epsdim;
                            itmp[4]=abs(iepsmp[iz2][iy2][ix2-1].nX)/epsdim;
                            itmp[5]=abs(iepsmp[iz2][iy2-1][ix2].nY)/epsdim;
                            itmp[6]=abs(iepsmp[iz2-1][iy2][ix2].nZ)/epsdim;

                            if(itmp[1]==0) iext=1;
                            if(itmp[1]!=itmp[6]) ibgp=1;

                            for(cont=2; cont<=6; cont++)
                            {
                                if(itmp[cont]==0) iext=1;
                                if(itmp[cont]!=itmp[cont-1]) ibgp=1;
                            }

                            if ((ibgp==0)&&(bndeps[ix2][iy2][iz2][1]!=0))
                            {
                                //reset bndeps for that point (i.e. remove
                                //bgp flag).
                                //a bgp become internal
                                iBoundNumsurf=iBoundNumsurf-bndeps[ix2][iy2][iz2][2];
                                bndeps[ix2][iy2][iz2][1]=0;
                                bndeps[ix2][iy2][iz2][2]=0;
                                mr=mr+1;
                            }
                            else
                            {
                                if (ibgp==1&&iext==0&&bndeps[ix2][iy2][iz2][2]==1)
                                {
                                    //an ext bgp is turned into an internal bgp
                                    iBoundNumsurf=iBoundNumsurf-1;
                                    bndeps[ix2][iy2][iz2][2]=0;
                                }
                            }

                            if (ibgp==1&&bndeps[ix2][iy2][iz2][1]==0)
                            {
                                //create a new boundary point..
                                m=m+1;
                                bndeps[ix2][iy2][iz2][1]=n2+m;
                                if(n2+m > ibmx)
                                {
                                    cout << " WARNING !!! This case is too big, ibmx need to be encreased." << endl; // Lin Li
                                    exit(0);
                                }
                                ibnd[n2+m]=int_coord(ix2,iy2,iz2);
                                bndeps[ix2][iy2][iz2][2]=iext;
                                iBoundNumsurf=iBoundNumsurf+bndeps[ix2][iy2][iz2][2];
                            }
                        }
                        //now jump to the next midpoint of the same grid point
                    }

                    //remap iepsmp in case there have been changes.. (that is some 0's became -1's)
                    //in other words: midpoint must remain external to objects
                    for(jj=1; jj<=6; jj++)

                    {
                        //in this way I can deal with eps[jj]<0
                        isign=1;

                        //iord=owner of the midpoint jj before change or after eps=1
                        iiord=optComp(iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]],imap[4][jj]);
                        iord=(iiord%epsdim);

                        //the last changed for sure has not iord<0 there can be iord<0 due to nearest neighbors already changed
                        if (iord<0) continue;

                        //if it has changed at previous step, dont change anymore
                        if (eps[jj]<0)
                        {
                            isign=-1;
                            if (iord==0) iord=1;
                        }

                        jjj=abs(iiord)/epsdim;

                        switch (imap[4][jj])
                        {
                            case 1:
                                iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nX=isign*(iord+jjj*epsdim);
                                break;
                            case 2:
                                iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nY=isign*(iord+jjj*epsdim);
                                break;
                            case 3:
                                iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nZ=isign*(iord+jjj*epsdim);
                                break;
                            default:
                                break;
                        }

                        //left iord definition with mod since if it is <> 0, it keeps its identity
                    }

                    //at this point one still can trace what changes have been made check to see if this is still a boundary point
                    //once again one distinguishes between internal,external,internal bgp and external bgp
                    iext=0;
                    ibgp=0;

                    itmp[1]=abs(iepsmp[iz][iy][ix].nX)/epsdim;
                    itmp[2]=abs(iepsmp[iz][iy][ix].nY)/epsdim;
                    itmp[3]=abs(iepsmp[iz][iy][ix].nZ)/epsdim;
                    itmp[4]=abs(iepsmp[iz][iy][ix-1].nX)/epsdim;
                    itmp[5]=abs(iepsmp[iz][iy-1][ix].nY)/epsdim;
                    itmp[6]=abs(iepsmp[iz-1][iy][ix].nZ)/epsdim;

                    if(itmp[1]==0) iext=1;
                    if(itmp[1]!=itmp[6]) ibgp=1;

                    for(cont=2; cont<=6; cont++)
                    {
                        if(itmp[cont]==0) iext=1;
                        if(itmp[cont]!=itmp[cont-1]) ibgp=1;
                    }

                    //if not now a boundary element change bndeps
                    if ((iext==0)||(ibgp==0))
                    {
                        iBoundNumsurf=iBoundNumsurf-bndeps[ix][iy][iz][2];
                        if(ibgp==1) bndeps[ix][iy][iz][2]=iext;

                        if (ibgp==0)
                        {
                            bndeps[ix][iy][iz][1]=0;
                            bndeps[ix][iy][iz][2]=0;
                            mr=mr+1;
                            if(iext==1)cout <<" //!!born a new external point!!!" << endl;
                        }
                    }

                }//if end for whether bndeps is nonzero;

            }//next boundary point FINISH;

            n1=n2+1;
            n2=n2+m;

            #ifdef VERBOSE
            cout <<" VdMS> bgp added m=" << m << " bgp removed mr =" << mr << endl;
            cout <<" VdMS> bgp added m=" << m << " bgp removed mr =" << mr << endl;
            #endif // VERBOSE

            if (m>mpr)
            {
                ndv=ndv+1;
                if (ndv>20)   //Lin Li: the value used to be 2,;
                {
                    // sometimes not enough
                    cout <<" WARNING !!! Surface iteration did not converge" << endl;
                    exit (0);
                }
            }
            else
            {
                ndv=0;
            }

            //2011-05-18 Replaced goto 100 statement
            if(m<=0)
            {
                break;
            }
        }

        if (n2>ibmx)
        {
            cout <<" WARNING !!! ibnd upper bound " << n2 << " exceeds ibmx" << endl;
            exit (0);
        }

        if(cbn1_v.size()>0) vector <delphi_integer>().swap(cbn1_v);
        if(cbn2_v.size()>0) vector <delphi_integer>().swap(cbn2_v);
        if(cbal.size()>0) vector <delphi_integer>().swap(cbal);

        if(cbn1 != NULL) free_pt3d<delphi_integer>(cbn1,lcb+1,mcb+1,ncb+1);
        if(cbn2 != NULL) free_pt3d<delphi_integer>(cbn2,lcb+1,mcb+1,ncb+1);


        #ifdef VERBOSE
        cout <<" VdMS> Number of cavity mid-points inaccessible to solvent = " << ncav << endl;
        #endif // VERBOSE

        //consolidate the list, removing dead boundary points, adding new ones..
        j=0;

        //2011-05-19 Array is re-sized keeping old values in the memory.
        if(ibgrd_v.size() > 0)
        {
            Nar=sizeof(ibgrd_v);
            if (Nar<ibmx)
            {
                ibgrd_v.resize(ibmx);
                ibgrd=&ibgrd_v[0]-1;
            }
        }
        else
        {
            ibgrd_v.assign(ibmx+1,sgrid_temp_int);
            ibgrd=&ibgrd_v[0]-1;
        }

        for(i=1; i<=n2; i++)
        {
            ixyz=ibnd[i];
            ix=ixyz.nX;
            iy=ixyz.nY;
            iz=ixyz.nZ;

            if (bndeps[ix][iy][iz][1]!=0)
            {
                j=j+1;
                bndeps[ix][iy][iz][1]=j;

                //2011-05-19 Precaution not to exceed array size (see above comment)
                if (j<=ibmx)
                {
                    ibgrd[j]=ixyz;
                }
                else
                {
                    cout << " VdMS> j=" << j << " is larger than ibmx= " << ibmx << " << thus stopped..." << endl;
                    exit (0);
                }
            }

            for(jj=1; jj<=6; jj++)
            {
                ntt=optComp(iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]],imap[4][jj]);
                nt=(ntt%epsdim);

                if (nt<0)
                {
                    ntt=-ntt;

                    switch (imap[4][jj])
                    {
                        case 1:
                            iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nX=ntt;
                            break;
                        case 2:
                            iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nY=ntt;
                            break;
                        case 3:
                            iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nZ=ntt;
                            break;
                        default:
                            break;
                    }

                    if (nt==-1)
                    {
                        ntt=ntt-1;

                        switch (imap[4][jj])
                        {
                            case 1:
                                iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nX=ntt;
                                break;
                            case 2:
                                iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nY=ntt;
                                break;
                            case 3:
                                iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nZ=ntt;
                                break;
                            default:
                                break;
                        }
                    }
                }
            }
        }

        if (j>ibmx)
        {
            cout << " WARNING !!! Number of  MS points exceeds ibmx" << endl;
            exit (0);
        }

        iBoundNum=j;

        #ifdef VERBOSE
        cout <<" VdMS> After surface elaboration iBoundNum= " << iBoundNum << " and iBoundNumsurf= " << iBoundNumsurf << endl;
        #endif // VERBOSE

    }

    delete [] ibnd;
    if(bndeps != NULL) free_pt4d<delphi_integer>(bndeps,iGrid.nX+1,iGrid.nY+1,iGrid.nZ+1,3);

    //fScale bondary grid point positions relative to acc data
    if (isolv&&(irea||logs||lognl||isen||isch))
    {
        #ifdef VERBOSE
        cout <<" VdMS> Scaling boundary grid points ..." << endl;
        #endif // VERBOSE

        scspos_v.assign(iBoundNum,sgrid_temp_real);
        scspos=&scspos_v[0]-1;

        for(j=1; j<=iBoundNum; j++)
            scspos[j]=optCast<delphi_real,delphi_integer>(ibgrd[j]);

        scsnor_v.assign(iBoundNum,sgrid_temp_real);
        scsnor=&scsnor_v[0]-1;

        atsurf_v.assign(iBoundNum,0);
        atsurf=&atsurf_v[0]-1;

        atndx_v.assign(iBoundNum,0);
        atndx=&atndx_v[0]-1;

        sclbp();

        #ifdef VERBOSE
        cout << " Info> " << iall << " points had to be assigned by global comparison" << endl;
        #endif // VERBOSE

        if (!isite && scsnor_v.size() >0 ) vector<SGrid <delphi_real> >().swap(scsnor_v); //no need to deallocate for vectors
    }

    if (isrf)
    {
        if (bOnlyMol)
        {
            get_pt3d <SGrid <delphi_integer> > (egrid,iGrid.nX+1,iGrid.nY+1,iGrid.nZ+1);

            if(egrid != NULL) free_pt3d(egrid, iGrid.nX + 1, iGrid.nY + 1, iGrid.nZ + 1);
        }
        else
        {
            #ifdef VERBOSE
            cout << " WARNING !!! msrf routine cannot be run" << endl;
            cout << " WARNING !!! because there are also geometric objects" << endl;
            #endif
        }
    }

    if (!isitsf&&!isite&&!(isch&&scrgfrm!=0))
    {
        if(atndx_v.size()>0) vector <delphi_integer>().swap(atndx_v);
        if(atsurf_v.size()>0) vector <delphi_integer>().swap(atsurf_v);

    }// if

    if(iab1 != NULL) free_pt3d<delphi_integer>(iab1,lcb1+1,mcb1+1,ncb1+1);
    if(iab2 != NULL) free_pt3d<delphi_integer>(iab2,lcb1+1,mcb1+1,ncb1+1);

    if(cbn1 != NULL) free_pt3d<delphi_integer>(cbn1,lcb+1,mcb+1,ncb+1);
    if(cbn2 != NULL) free_pt3d<delphi_integer>(cbn2,lcb+1,mcb+1,ncb+1);

    if(r0.size()>0) vector <delphi_real>().swap(r0);
    if(r02.size()>0) vector <delphi_real>().swap(r02);
    if(rs2.size()>0) vector <delphi_real>().swap(rs2);
    if(ast.size()>0) vector <delphi_integer>().swap(ast);

    ibgrd_v.resize(iBoundNum);
    ibgrd=&ibgrd_v[0]-1;

    //At the very end of things
    #ifdef VERBOSE
    cout << " VdMS> MS creation done" << endl;
    #endif
}// void vwtms;
