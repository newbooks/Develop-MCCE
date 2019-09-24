#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"

#include "space.h"

using namespace std;

void CDelphiSpace::epsmak()
{
    SGrid <delphi_real> amin,amax;
    SExtrema<delphi_real> sextrema_temp;
    delphi_real fRMaxTemp;
    delphi_integer i,ix,iy,iz;
    sextrema_temp.nMax=sgrid_temp_real;
    sextrema_temp.nMin=sgrid_temp_real;

    sLimGridUnit.assign(iNObject, {});

    /*
     * here limobject is expressed in grid units
     */
    for(i=0; i<=iNObject-1; i++)
    {
        if (ibctyp != 3)
        {
            // (Parallel)
            sLimGridUnit[i].nMin = (sLimObject[i].nMin - global_cOldMid)*fScale + fRMid - optCast<delphi_real, delphi_integer>(myStart);
            sLimGridUnit[i].nMax = (sLimObject[i].nMax - global_cOldMid)*fScale + fRMid - optCast<delphi_real, delphi_integer>(myStart);
        }
        else
        {
            delphi_integer iGrid1D = pdc->getKey_Val<delphi_integer>("igrid");
            sLimGridUnit[i].nMin = (sLimObject[i].nMin - global_cOldMid)*fScale + (iGrid1D+1)/2.0 - optCast<delphi_real, delphi_integer>(myStart);
            sLimGridUnit[i].nMax = (sLimObject[i].nMax - global_cOldMid)*fScale + (iGrid1D+1)/2.0 - optCast<delphi_real, delphi_integer>(myStart);
        }
    }

    if(bUniformDiel) // ARGO: make it work only for HOMO models
    {
        cout << "not going to calculate boundary elements since" << endl;
        cout << "uniform dielectric" << endl;
        iBoundNum=0;
        return;
    }

    /*
     * lepsx.y.z and uepsx.y.z should be the upper and lower limits of the expanded box. If the molecule
     * is smaller than this then reduce leps and upeps accordingly note leps/ueps not yet defined..
     */

    // 2011-05-10 Converted to SGrid <float> derived type
    amin = sLimGridUnit[0].nMin;
    amax = sLimGridUnit[0].nMax;

    /**
     * find global limits IN GRID UNITS, both, molecule and objects, are considered
     */
    #ifdef VERBOSE
    if(iNObject > 1)
        for(i=1; i<=iNObject-1; i++)
            cout << "to be finished using SGrid <float> operations:" << endl;
    #endif

    fRMaxTemp = rdmx;
    if(fIonStrenth !=0 ) fRMaxTemp=max(fRMaxTemp,fExternRadius);

    fRMaxTemp=fRMaxTemp*fScale;

    // Using operations on SGrid <float> type variables defined in module operators_on_coordinates
    amin=amin-fRMaxTemp;
    amax=amax+fRMaxTemp;

    // determine the limit of limEps, no larger than igrdd and no less than 1.
    LimEps.nMin=optMax(optCast <delphi_integer,delphi_real> (amin) - delphi_integer(2), delphi_integer(1));
    LimEps.nMax=optMin(optCast <delphi_integer,delphi_real> (amax) + delphi_integer(2), iGrid);

    /**
     *point is out of any kind of object (probably in solution):
     *If radprb is less than half of grid spacing, then use old algorithm (sri 29March 93);
     *The new algorithm should be able to handle all scales; still remains to be tested (Sri Apr 95)
     */
    if ( iGaussian==0 )
    {
        setout();

        VdwToMs();
    }
    else
        setGaussian();
}





