#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"


#include "../space/space.h"

using namespace std;

void CDelphiSpace::cubedata(delphi_real fac, delphi_real cbln)
{
    delphi_real off;
    SGrid <delphi_real> xyzp;
    off=0.1;

    // fac=1.0 was passed as parameter, now removed passing process and set fac=1.0 below:
    xyzo=cMin-((fac*cbln)+off);
    xyzp=cMax+((fac*cbln)+off);

    lmncb=optCast<delphi_integer,delphi_real> ((xyzp-xyzo)/cbln);
    lcb=lmncb.nX;
    mcb=lmncb.nY;
    ncb=lmncb.nZ;

    cbai=1./cbln;

}// void cubedata;

//---------------------------------------------------------------------
void CDelphiSpace::cube()
{
    // rda changed to sDelPhiPDB[].radius
    //here rda equals rad3
    //2011-05-27 Arrays declared in pointers and allocated in calling
    //void
    //delphi_integer cbn1[0:lcb,0:mcb,0:ncb],cbn2[0:lcb,0:mcb,0:ncb];
    //delphi_integer ***cbn1,***cbn2;

    get_pt3d<delphi_integer>(cbn1,lcb+1,mcb+1,ncb+1);
    get_pt3d<delphi_integer>(cbn2,lcb+1,mcb+1,ncb+1);

    delphi_integer newatm,itmp,kind;

    /**
     * creating a set of fictitious atoms occupying little cubes
     */
    SGrid <delphi_integer> icbn[iNatom+(iNObject-numbmol)*(lcb+1)*(mcb+1)*(ncb+1)+1];

    //iatmobj connects the fictitious atoms to the objects
    delphi_integer iatmobj[(iNObject-numbmol)*(lcb+1)*(mcb+1)*(ncb+1)];

    string strtmp;

    SGrid <delphi_real> xq,xyz;
    SGrid <delphi_integer> ixyz;

    delphi_real shift;
    delphi_real cost;

    delphi_real prbrd,cbln;
    delphi_integer i,j,k,ii,ix,iy,iz,jx,jy,jz,icum;

    for(i=0;i<=iNatom+(iNObject-numbmol)*(lcb+1)*(mcb+1)*(ncb+1);i++)
        icbn[i]=sgrid_temp_int;

    prbrd=fRadPrb[1];

    cbln=1./cbai;

    for (i=0;i<=lcb;i++)
    {
        for (j=0;j<=mcb;j++)
        {
            for (k=0;k<=ncb;k++)
            {
                cbn1[i][j][k]=1; //this is in [x][y][z]
                cbn2[i][j][k]=0; //this is in [x][y][z]
            }
        }
    }

    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            xyz=(xn1[i]-xyzo)*cbai; //xyzo variable is declared in pointers module, cbai in qlog
            ixyz=optCast<delphi_integer,delphi_real>(xyz);

            ix=ixyz.nX;
            iy=ixyz.nY;
            iz=ixyz.nZ;

            for(jz=iz-1; jz<=iz+1; jz++)
                for(jy=iy-1; jy<=iy+1; jy++)
                    for(jx=ix-1; jx<=ix+1; jx++)
                        cbn2[jx][jy][jz]=cbn2[jx][jy][jz]+1;

            icbn[i]=ixyz;
        }
    }

    newatm=iNatom;
    cost=cbln*.87+1./fScale;
    shift=cost+prbrd;

    for(ii=1; ii<=iNObject; ii++) // icbn will contain also coord center of fictious atoms
    {
        strtmp=dataobject_v[(ii-1)*2];

        kind = atoi(strtmp.substr(15,3).c_str());
        if (strtmp.substr(0,4) != "is a" && kind != 2)
        {
            itmp=ii+iNatom;

            for(iz=0; iz<=ncb; iz++)
            {
                for(iy=0; iy<=mcb; iy++)
                {
                    for(ix=0; ix<=lcb; ix++)
                    {
                        ixyz=int_coord(ix,iy,iz);
                        xq=((optCast<delphi_real,delphi_integer>(ixyz)+0.5)*cbln)+xyzo;
                    }
                }
            }
        }
    }

    icum=1;
    for(iz=0; iz<=ncb; iz++)
    {
        for(iy=0; iy<=mcb; iy++)
        {
            for(ix=0; ix<=lcb; ix++)
            {
                if (cbn2[ix][iy][iz]>0)
                {
                    cbn1[ix][iy][iz]=icum;
                    icum=icum+cbn2[ix][iy][iz];
                }
            }
        }
    }

    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
        }
    }

    for(i=iNatom+1; i<=newatm; i++) //This will not be excuted
    {
        cout << "Warning: newatm > iNatom" << endl;
        ix=icbn[i].nX;
        iy=icbn[i].nY;
        iz=icbn[i].nZ;
        cbal[cbn1[ix][iy][iz]]=iatmobj[i-iNatom];
        cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
    }

    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //1,0,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //0,-1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-1;
            iy=iy-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
            icbn[i].nY=iy;
        }
    }

    //0,1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }
    }

    //0,0,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy-1;
            iz=iz-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
            icbn[i].nZ=iz;
        }
    }

    //0,0,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iz=iz+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nZ=iz;
        }
    }

    //nn=2
    //1,0,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //-1,0,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //0,1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+1;
            iy=iy+1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
            icbn[i].nY=iy;
        }
    }

    //0,-1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }
    }

    //-1,-1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-1;
            iz=iz-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
            icbn[i].nZ=iz;
        }
    }

    //1,-1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //1,1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }
    }

    //-1,1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //-1,0,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iz=iz-1;
            iy=iy-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
            icbn[i].nZ=iz;
        }
    }

    //1,0,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //0,1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-1;
            iy=iy+1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
            icbn[i].nY=iy;
        }
    }

    //0,-1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }
    }

    //nn=3
    //-1,-1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //1,-1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //1,1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }
    }

    //-1,1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //-1,1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iz=iz+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nZ=iz;
        }
    }

    //1,1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //1,-1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }
    }

    //-1,-1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }
    }

    //reset cbn1
    icum=1;
    for(iz=0; iz<=ncb; iz++)
    {
        for(iy=0; iy<=mcb; iy++)
        {
            for(ix=0; ix<=lcb; ix++)
            {
                if (cbn2[ix][iy][iz]>0)
                {
                    cbn1[ix][iy][iz]=icum;
                    icum=icum+cbn2[ix][iy][iz];
                    cbn2[ix][iy][iz]=icum-1;
                }
            }
        }
    }

    icum=icum-1;
}// void cube;


