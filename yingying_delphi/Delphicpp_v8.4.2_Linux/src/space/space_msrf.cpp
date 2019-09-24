#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"

#include "space.h"

using namespace std;

void CDelphiSpace::msrf()
{
    delphi_integer itot, vtot;
    string fixyn;
    string vdat1, line[6], fname;
    bool fix;
    SGrid <delphi_real> v1,v2,vxyz;
    SGrid <delphi_integer> iv123;

    delphi_integer mxvtx,i,ia1,ia2,ia3,iate,it,imxvtx;
    delphi_integer ib,iv1,iv2,iv3,j,k,mxtri,ntot2;
    delphi_real aa,area,cc,bb,rad,rad2,ss;
    delphi_real tar,tne4,vmg;

    iv1=0;
    iv2=0;
    iv3=0;
    vtot=0;
    itot=0;
    ntot2=0;

    mxvtx=iBoundNum*2.2;
    mxtri=2*mxvtx;

    //epsdim=iNatom+iNObject+2; replaced in space module semi global

    for(k=1; k<=iGrid.nZ; k++)
    {
        for(j=1; j<=iGrid.nY; j++)
        {
            for(i=1; i<=iGrid.nX; i++)
            {
                //changed from div mod, it should only serve the means
                egrid[i][j][k]=iepsmp[i][j][k]/epsdim;
            }
        }
    }

    for(i=1; i<=iGrid.nX; i++)
    {
        for(j=1; j<=iGrid.nY; j++)
        {
            egrid[i][1][j].nX = egrid[i][1][j].nY;
            egrid[i][j][1].nX = egrid[i][j][1].nZ;
        }
    }

    for(k=2; k<=iGrid.nZ; k++)
    {
        for(j=2; j<=iGrid.nY; j++)
        {
            for(i=2; i<=iGrid.nX; i++)
            {
                iate = 0;
                if (egrid[i][j][k].nX > 0) iate = iate + 1;
                if (egrid[i][j][k].nY > 0) iate = iate + 1;
                if (egrid[i][j][k].nZ > 0) iate = iate + 1;
                if (egrid[i-1][j][k].nX > 0) iate = iate + 1;
                if (egrid[i][j-1][k].nY > 0) iate = iate + 1;
                if (egrid[i][j][k-1].nZ > 0) iate = iate + 1;

                if (iate<=3)
                    egrid[i][j][k].nX = 0;
                else
                    egrid[i][j][k].nX = 1;
            }
        }
    }

    vindx.assign(mxtri+1, sgrid_temp_int);
    vert.assign(mxvtx+1, sgrid_temp_real);

    #ifdef VERBOSE
    cout << "mxtri,mxvtx= " << " " << mxtri << " " << mxvtx << endl;
    #endif

    vdat1="./";

    if (vtot>mxvtx)
    {
        cout <<"vtot = " << vtot << " > mxvtx = " << mxvtx << endl;
        cout <<"increase mxvtx in msrf.f" << endl;
        exit(0);
    }

    for(ib=1; ib<=vtot; ib++)
    {
        vert[ib]=vert[ib]/2.;
    }

    itot = itot/3;

    //fScale boundary grid pointeger positions relative to acc data
    #ifdef VERBOSE
    cout <<"scaling vertices" << endl;
    #endif

    vnorm.assign(vtot+1, sgrid_temp_real);
    vnorm2.assign(vtot+1, sgrid_temp_real);

    //fix holes and make vertex to triangle arrays allocate hole-fixing arrays next hole-fixing variables
    if (vtot < mxvtx/2)
        imxvtx = vtot*2;
    else
        imxvtx = mxvtx;

    vtlen.assign(imxvtx+1,0);
    vtpnt.assign(imxvtx+1,0);
    vtlst.assign(6*imxvtx+1,0);

    get_pt2d <delphi_integer> (tmlst,9+1,imxvtx+1);

    if (ntot2 > 0)
        fix = true;
    else
        fix = false;

    while(fix)
    {
        if (ntot2 > 0)
            fix = true;
        else
            fix = false;
    }

    if (itot>mxtri)
    {
        cout <<"itot = " << itot << " > mxtri = " << mxtri << endl;
        cout <<"increase mxtri in msrf.f" << endl;
        exit(0);
    }

    if(vtlen.size()>0) vector <delphi_integer>().swap(vtlen);
    if(vtlst.size()>0) vector <delphi_integer>().swap(vtlst);
    if(vtpnt.size()>0) vector <delphi_integer>().swap(vtpnt);
    if(tmlst != NULL) free_pt2d(tmlst,9+1,imxvtx+1);


    #ifdef VERBOSE
    cout << "number of vertices = " << vtot << endl;
    cout << "number of triangles = " << itot << endl;
    #endif

    //calculate area
    area=0.0;

    for(it=1; it<=itot; it++)
    {
        iv123=vindx[it];
        v1=vert[iv2]-vert[iv1];
        v2=vert[iv3]-vert[iv1];

        //2011-05-20 Vector product defined in operators_on_coordinates module
        vxyz=optCross(v1,v2);
        vmg=sqrt(optDot(vxyz,vxyz));
        tar=vmg/2.;
        vxyz=vnorm[iv1]+vnorm[iv2]+vnorm[iv3];
        vmg=sqrt(optDot(vxyz,vxyz));

        vnorm2[iv1]=vnorm2[iv1]+(vxyz/vmg);
        vnorm2[iv2]=vnorm2[iv2]+(vxyz/vmg);
        vnorm2[iv3]=vnorm2[iv3]+(vxyz/vmg);

        //calculate spherical triangle area if appropriate
        ia1=atndx[iv1];
        ia2=atndx[iv2];
        ia3=atndx[iv3];

        if (ia1>0)
        {
            if (ia1==ia2&&ia1==ia3)
            {
                rad=sDelPhiPDB[ia1].radius;
                rad2=rad*rad;
                aa=optSum((vert[iv2]-vert[iv1])*(vert[iv2]-vert[iv1]));
                bb=optSum((vert[iv3]-vert[iv2])*(vert[iv3]-vert[iv2]));
                cc=optSum((vert[iv1]-vert[iv3])*(vert[iv1]-vert[iv3]));

                aa=acos(1.-aa/(2.*rad2));
                bb=acos(1.-bb/(2.*rad2));
                cc=acos(1.-cc/(2.*rad2));
                ss=(aa+bb+cc)*.5;
                tne4=sqrt(tan(ss*.5)*tan((ss-aa)*.5)*tan((ss-bb)*.5)*tan((ss-cc)*.5));
                tar=4.*atan(tne4)*rad2;
            }
        }

        area=area+tar;
    }

    for(i=1; i<=vtot; i++)
    {
        vmg=sqrt( optDot( vnorm2[i],vnorm2[i] ) );
        vnorm2[i]=vnorm2[i]/vmg;
    }

    #ifdef VERBOSE
    cout <<"MS area = " << area << endl;
    #endif

    //2011-05-26 Subroutine wrtsurf is short and called only once,
    // therefore put the body of the void here

    //-------- Begin of wrtsurf body -----------------------------

    if (!ibem)
    {
        ofstream surfile;
        surfile.open ("grasp.srf");

        fname="grasp.srf";

        #ifdef VERBOSE
        cout << "writing GRASP file to " << fname << endl;
        #endif

        line[1]="format=2";
        line[2]="vertices,normals,triangles";

        surfile << "format=2" << endl;
        surfile << "vertices,normals,triangles" << endl;
        surfile << endl;
        surfile <<  setw(6) << vtot << setw(6) << itot << setw(6) << iGrid << setw(12) << setprecision(6) << fScale;
        surfile << setw(6) << setw(12) << setprecision(6) << cOldMid.nX
                << setw(6) << setw(12) << setprecision(6) << cOldMid.nY
                << setw(6) << setw(12) << setprecision(6) << cOldMid.nZ;

        #ifdef VERBOSE
        cout << "writing data for" << vtot << " vertices and" << itot << " triangles" << endl;
        #endif

        for(i=1;i<=mxvtx;i++)
            surfile << vert[i] << endl;

        for(i=1;i<=vtot;i++)
            surfile << vnorm[i] << endl;

        for(i=1;i<=mxtri;i++)
            surfile << vindx[i] << endl;

        surfile.close();

        #ifdef VERBOSE
        cout << "finished writing " << fname << endl;
        #endif
    }
    else
    {
        ofstream surfile;
        surfile.open ("bem.srf");

        surfile << vtot << " " << itot << endl;

        for(i=1; i<=vtot; i++)
            surfile << vert[i] << endl;

        for(i=1; i<=itot; i++)
            surfile << vindx[i] << endl;

        for(i=1; i<=vtot; i++)
            surfile << vnorm[i] << endl;

        surfile.close();
    }
    //-------- End of wrtsurf body --------------------------------

}// void msrf;
