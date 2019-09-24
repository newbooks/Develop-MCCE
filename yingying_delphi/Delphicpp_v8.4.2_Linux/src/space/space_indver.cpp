#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"

#include "space.h"

using namespace std;

void CDelphiSpace::indverdata(delphi_real prbrad)
{
    delphi_real jig;
    delphi_real cba;
    delphi_integer idtemp;

    /*
     *cba defines fScale, cba and jig the box limits..
     *note that cba is the minumum box size. it can be bigger if need be. i.e. can refScale to a larger box if need be.
     */
    cba=prbrad;

    jig=2.0*(1.0/fScale);

    //2011-05-24 Using operations on coord and int_coord type
    mnxyz=cMin-(fRMax+prbrad);
    mxxyz=cMax+(fRMax+prbrad);

    /**
     *accessible surface points are prbrad away from the actual surface mid points are at most (1/fScale) away from the actual surface.
     *midpoints must never be in box zero. accessible pointer can be in box zero but not in box -1.
     *
     *e.g. mnx+prbrad==van der Waals surface van der Waals surface - (1/fScale)= min. midpointeger position.
     *therefore mnx+prbrad-(1/fScale) gt 1 and mnx gt 0. the latter is always true since cba is ge. prbrad.
     *therefore we have...
     */
    mnxyz=mnxyz-jig;
    mxxyz=mxxyz+jig;

    //2011-05-24 Cycle is introduced to remove GOTO 100 statement
    while(true)
    {
        mxxyz=mxxyz+(cba-prbrad);
        mnxyz=mnxyz-(cba-prbrad);
        lmncb1=optCast <delphi_integer,delphi_real> (((mxxyz-mnxyz)/cba)+1.);

        /**
         *if points are too widely separated for the fScale and idmax
         *refScale..
         */

        //2011-05-24 idmax is declared as parameter in qlog module
        if (optORGT(lmncb1, delphi_integer(iDMax)))
        {
            idtemp=optMax(lmncb1);
            cba=cba*Int2Float(idtemp+1)/iDMax;
            cout << "initial cube size too small << " << endl;
            cout << "in assigning accessible points to a grid" << endl;
            cout << "therefore rescaling..." << endl;

            //2011-05-24 CYCLE and EXIT statements instead of GOTO
            continue;
        }
        else
        {
            break;
        }
    }

    lcb1=lmncb1.nX;
    mcb1=lmncb1.nY;
    ncb1=lmncb1.nZ;

    //grdi is just the inverse of cba..,used more...
    grdi=1.0/cba;
}// void indverdata;

//---------------------------------------------------------------------

void CDelphiSpace::indver(delphi_integer extot1)
{

    /**
     * program to compile the lists iab1,iab2 and icume for use in
     * nearest vertex work. iexpos are box coordinates of vertices
     * and don't need to be passed if comparisons are done with float
     * angstroms..
     * but its often more convenient to do so since grid points have to
     * be converted anyway...
     */

    //2011-05-24 Arrays are accessible via declaration in pointers
    //module and allocation in calling vwtms void
    get_pt3d<delphi_integer>(iab1,lcb1+1,mcb1+1,ncb1+1);
    get_pt3d<delphi_integer>(iab2,lcb1+1,mcb1+1,ncb1+1);

    //2011-05-27 Declarations added due to IMPLICIT NONE
    delphi_integer i,n,ix,iy,iz,k,j;
    vector < SGrid <delphi_integer> > iexpos; //iexpos now is a local variable;

    //2011-05-24 Array iexpos now is of int_coord type and thus 1D
    iexpos.assign(extot1+1, sgrid_temp_int);

    //initialize grid..
    //2011-05-24 Changed to array operation
    for (i=0; i<=lcb1; i++)
    {
        for (j=0; j<=mcb1; j++)
        {
            for (k=0; k<=ncb1; k++)
            {
                iab1[i][j][k]=1;
                iab2[i][j][k]=0;
            }
        }
    }

    /*
     * make linear arrays for expos
     *find the number of points in each box, put in iab2, make iexpos
     */
    for(i=1; i<=extot1; i++)
    {
        //2011-05-24 Using operations on coord and int_coord type
        //variables defined in module operators_on_coordinates
        iexpos[i]=optCast <delphi_integer,delphi_real> ( (expos[i]-mnxyz)*grdi );
        iab2[iexpos[i].nX][iexpos[i].nY][iexpos[i].nZ]= iab2[iexpos[i].nX][iexpos[i].nY][iexpos[i].nZ]+1;
    }

    //check each box for occupancy, using fill number to mark out space in icume
    n=0;
    for(i=0; i<=lcb1; i++)
    {
        for(j=0; j<=mcb1; j++)
        {
            for(k=0; k<=ncb1; k++)
            {
                //if the box is not empty put start position to n+1 in iab1
                // to n+box occupancy in iab2, overwriting occupancy..
                if (iab2[i][j][k]!=0)
                {
                    iab1[i][j][k]=n+1;
                    n=n+iab2[i][j][k];
                    iab2[i][j][k]=n;
                }
            }
        }
    }

    //fill icume using iab1 and iab2, note that iab1 is used to hold
    //the position in icume, therefore needs to be recalculated..
    for(i=1; i<=extot1; i++)
    {
        ix=iexpos[i].nX;
        iy=iexpos[i].nY;
        iz=iexpos[i].nZ;
        j=iab1[ix][iy][iz];
        icume[j]=i;
        iab1[ix][iy][iz]=iab1[ix][iy][iz]+1;
    }

    //reset iab1 for use in inner loop
    for(i=1; i<=extot1; i++)
    {
        ix=iexpos[i].nX;
        iy=iexpos[i].nY;
        iz=iexpos[i].nZ;
        iab1[ix][iy][iz]=iab1[ix][iy][iz]-1;
    }

    /*
     * icume now contains pointers to each dot inside a particular box and each box has
     * 2 pointers into icume.,a start pointer and a finish pointer
     * use however you want...
     */

    //2011-05-24 Array deallocation is done with DEALLOCATE statement
    if(iexpos.size() > 0 ) vector < SGrid <delphi_integer> > ().swap(iexpos);

}// void indver;

