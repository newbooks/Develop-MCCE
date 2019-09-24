#ifndef SPACE_H
#define SPACE_H

#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <memory>
#include <valarray>
#include <algorithm>
//#include <deque>

#include <cmath> // std::abs

#include "../interface/interface_abstractmodule.h"
#include "../misc/misc_timer.h"
#include "../misc/misc_grid.h"
#include "../delphi/delphi_constants.h"
#include "../io/io.h"
#include "../interface/interface_datacontainer.h"
#include "space_templates.h"
#include "space_exceptions.h"
#include <stdio.h>
#include <stdlib.h>

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
using namespace std;

class CDelphiSpace:virtual public IAbstractModule
{
private:
    shared_ptr<CTimer> pTimer;

    /*********************************************************************************************
     *                                                                                           *
     *              references to the variables obtained from the data container                 *
     *                                                                                           *
     ********************************************************************************************/

    //++++++++++++++ const references to read-only variables from data container +++++++++++++++//

    delphi_integer iNatom;                           // natom
    const delphi_real&    fScale;                    // scale
    SGrid<delphi_integer> iGrid;                     // igrid (modified in setFocusBndy)
    SGrid<delphi_real> fRMid;
    const delphi_integer& iNObject;                  // nobject
    const delphi_real&    repsout;                   // repsout
    const delphi_real&    repsin;                    // repsin
    const delphi_real& fDebyeLength;
    const delphi_real& fEpsOut;
    const delphi_real& fEpsIn;
    const delphi_real& repsout2;

    //ARGO 14-FB,2016
    //To use for vdwSurf_file to be written as a CUBE file
    const SGrid<delphi_real>& fgBoxCenter;           // oldmid
    const bool& bUniformDiel;
    const bool& bOnlyMol;
    const bool& isolv;
    const bool& irea;
    const bool& logs;
    const bool& lognl;
    const bool& isen;
    const bool& isch;
    const bool& isite;
    const bool& ibem;
    const int& ibctyp;
    const bool& isitsf;
    const bool& bFrcAsPqr;

    const delphi_integer& iNMedia;
    const delphi_integer& numbmol;
    const delphi_integer& scrgfrm;
    SGrid<delphi_real> cOldMid;
    const delphi_real& fIonStrenth;
    const delphi_real& fExternRadius;
    const delphi_real& fRMax;
    vector < delphi_real >& fRadPrb_v;
    const delphi_real* fRadPrb;
    vector <string> & dataobject_v;

    //++++++++++++++++ reference to read-and-write variables from data container +++++++++++++++//
    delphi_integer ndistr;                           // ndistr
    delphi_integer &iBoundNum;                       // ibnum
    delphi_real& rdmx;                               // rdmx
    delphi_integer& nqass;
    delphi_integer& nqgrd;
    bool &iacs;
    bool &isrf;
    SGrid<delphi_real> &cMin;
    SGrid<delphi_real> &cMax;
    delphi_real &qnet;
    delphi_real &qmin;
    delphi_real &qplus;
    SGrid<delphi_real>& cqmin;
    SGrid<delphi_real>& cqplus;
    SGrid<delphi_real>& acenter;                     //for focusing
    vector <delphi_real>& medeps;
    vector < SGrid<delphi_real> >& xn1_v;            //prgfgAtomCoordA
    vector < SGrid<delphi_real> >& xn2_v;            //prgfgAtomCoordG
    vector <char>&  bDebMap_v;                       //idebmap
    vector < SGrid<delphi_integer> > &iEpsMap_v;     //iepsmap
    vector < SGrid<delphi_real> >& fGepsMap_v;       //gepsmap
    vector < SGrid<delphi_real> >& fGepsMap2_v;      //gepsmap2
    vector < delphi_real >& fGDensityMap_v;          //gaussian density map
    vector <delphi_integer>& iAtomMed_v;             //iatmmed
    vector < SExtrema<delphi_real> > sLimObject;     //limobject
    vector <SGrid <delphi_integer> >& ibgrd_v;       //ibgrd
    vector < SGrid <delphi_real> >& scspos_v;
    vector < SGrid <delphi_real> >& chgpos_v;
    vector <delphi_integer>& crgatn_v;
    vector <delphi_integer>& nqgrdtonqass_v;
    vector <delphi_real>& atmeps_v;
    vector< SGridValue<delphi_real> >& atmcrg_v;
    vector< SGridValue<delphi_real> >& chrgv2_v;
    vector < SGrid <delphi_real> >& scsnor_v;
    vector < delphi_integer >& atsurf_v;
    vector < delphi_integer >& atndx_v;
    vector <CAtomPdb> delphipdb;                     //delphipdb
    const delphi_real& cutoff;
    const delphi_real&  sigma;
    int&  inhomo;
    const delphi_real&  srfcut;
    const delphi_real&  dencut;
    const int&  iGaussian;

    //+++++++++++++++++++ NON-refereces: +++++++++++++++++++++++++

    // Argo: added new members to the Struct. This is required for the Born Radius calculations.
    struct delphipdb_struc
    {
        SGrid <delphi_real> xyz;
        delphi_real charge;
        delphi_real radius;
        string atom_resname;        //ARGO
    };
    delphipdb_struc * sDelPhiPDB;

    //################# semi global variables in this class #########################
    delphi_integer iBoundNumsurf,extot,iall,lcb, mcb, ncb, lcb1, mcb1, ncb1;
    delphi_real radpmax,grdi, cbai, sideinter, sidemin;
    SGrid <delphi_real> mnxyz, xyzo, mxxyz;
    SGrid <delphi_integer> lmncb1, lmncb;
    SExtrema <delphi_integer> LimEps;
    delphi_integer extracrg;
    delphi_integer epsdim;

    vector <delphi_integer>  iab1_v, iab2_v, icume, ast, cbn1_v, cbn2_v, cbal, icbn;
    vector <delphi_real> r0,r02,rs2;
    vector < SGrid <delphi_real> > expos;
    vector < SExtrema<delphi_real> > sLimGridUnit;

    vector < SGrid <delphi_real> > vert, vnorm, vnorm2;
    vector < SGrid <delphi_integer> > vindx;
    vector <delphi_integer> vtlen, vtlst, vtpnt;

    delphi_integer ** tmlst;

    SGrid<delphi_real> * xn1;
    SGrid<delphi_real> * xn2;

    SGrid <delphi_integer> * ibgrd;
    SGridValue <delphi_real> * atmcrg;
    SGridValue <delphi_real> * chrgv2;
    SGrid <delphi_real> * scspos;
    SGrid <delphi_real> * chgpos;

    SGrid <delphi_real> * scsnor;
    delphi_integer * atsurf;
    delphi_integer * atndx;

    delphi_integer * crgatn;
    delphi_integer * nqgrdtonqass;
    delphi_real * atmeps;

    SGrid <delphi_integer> *** egrid;
    char *** idebmap;

    //ARGO: defining a new variable to store bools for extended vdw surface for zeta
    char *** zetaSurfMap;
    int& zetaOn;
    delphi_real& zetaDistance;
    vector <char>& zetaSurfMap_v;                    //zetaSurfMap
    const string&  strZetaPhiFile;                   //zphinam
    vector < delphi_real > dMoment;
    vector < delphi_real > qMoment;
    vector <delphi_real>& surf_grid_coords_v;
    vector <delphi_real>& surf_grid_index_v;

    //ARGO Born Radius
    delphi_real obc_cutoff = 6;
    delphi_real obc_a = 1, obc_b = 0.8, obc_g = 4.85;
    delphi_real Roffset = 0.09;
    bool bUseBornRadius = true;

    SGrid <delphi_integer> *** iepsmp;
    SGrid <delphi_real> *** gepsmp;
    SGrid <delphi_real> *** gepsmp2;

    /**
     * Gaussian Density Map array
     * Stores the density of atoms on each grid point
     *
     */
    delphi_real *** gDensityMapOnGridPoint;

    delphi_integer *** cbn1, *** cbn2, *** iab1, *** iab2;
    delphi_integer * iAtomMed;

    SGrid <delphi_real> sgrid_temp_real;
    SGrid <delphi_integer> sgrid_temp_int;
    SGrid <delphi_real> sgrid_rho_real;

    //############### Functions in other files: ######################
    void epsmak();
    void setout();
    void setGaussian();
    void VdwToMs();
    void VdwToMs_piece(bool&, const delphi_integer&, const delphi_integer&, const delphi_integer&,
            const delphi_integer&, const delphi_real&, const SGrid <delphi_real>&, delphi_integer*,
            delphi_integer*, delphi_integer&, delphi_real&);
    void sas();
    void cube();
    void cubedata(delphi_real, delphi_real);
    void indverdata(delphi_real);
    void indver(delphi_integer);
    void sclbp();
    void msrf();
    void crgarr();

    SGrid <delphi_integer> int_coord( const delphi_integer& a, const delphi_integer& b,  const delphi_integer& c);
    SGrid <delphi_real> coord( const delphi_real& a,  const delphi_real& b,  const delphi_real& c);

    SGrid <int> Float2Int( const SGrid <float>& a );
    int Float2Int(const float& a);

    SGrid <float> Int2Float( const SGrid <int>& a );
    float Int2Float(const int& a);

    shared_ptr<IDataContainer> test_pdc;

    //############### For parallel Space module: ######################
    size_t my_id;
    SGrid <delphi_integer> global_iGrid;
    SGrid<delphi_real> global_cOldMid;
    delphi_integer global_iNatom;
    vector<delphi_integer> globalAtomIndex;
    SGrid<delphi_integer> cpuDim;
    SGrid<delphi_integer> myCPULocation;

    //local grid size and location relative to the global grid
    SGrid <delphi_integer> myStart;
    SGrid <delphi_integer> myEnd;
    SGrid <delphi_integer> mySize;
    SGrid <delphi_integer> mySizeCoord;

    SGrid <delphi_real> myStartCoor;
    SGrid <delphi_real> myEndCoor;

    SGrid <delphi_real> halfSize;

    //local net grid size and location relative to the global grid
    SGrid <delphi_integer> myStartNoBuffer;
    SGrid <delphi_integer> myEndNoBuffer;
    SGrid <delphi_integer> mySizeNoBuffer;

    SGrid <delphi_integer> myInternalStartNoBuffer;
    SGrid <delphi_integer> myInternalEndNoBuffer;

    SGrid <delphi_real> myStartCoorNoBuffer;
    SGrid <delphi_real> myEndCoorNoBuffer;

    SGrid <delphi_integer> iGrid_temp;
    SGrid<delphi_integer> *** localGridDim;

    delphi_real max_radius;
    delphi_integer buffer;

public:
    SGrid <delphi_integer> *** iEpsMap;
    CDelphiSpace(shared_ptr<IDataContainer> pdc,shared_ptr<CTimer> pt):
        /*********************************************************************************************
         *                                                                                           *
         *              references to the variables obtained from the data container                 *
         *                                                                                           *
         ********************************************************************************************/

        //++++++++++++++ const references to read-only variables from data container +++++++++++++++//
        IAbstractModule(pdc),
        pTimer(pt),
        iNatom (pdc->getKey_Val<delphi_integer>("natom")),
        fScale (pdc->getKey_constRef<delphi_real>("scale")),
        iNObject (pdc->getKey_constRef<delphi_integer>("nobject")),
        repsout (pdc->getKey_constRef<delphi_real>("repsout")),
        repsout2 (pdc->getKey_constRef<delphi_real>("repsout2")),
        repsin (pdc->getKey_constRef<delphi_real>("repsin")),
        bUniformDiel (pdc->getKey_constRef<bool>("uniformdiel")),
        bOnlyMol (pdc->getKey_constRef<bool>("ionlymol")),
        isolv (pdc->getKey_constRef<bool>("isolv")),
        irea (pdc->getKey_constRef<bool>("irea")),
        
#ifdef PARALLEL_MPI
        logs(true),
#else
        logs (pdc->getKey_constRef<bool>("logs")),
#endif        
        lognl (pdc->getKey_constRef<bool>("lognl")),
        isen (pdc->getKey_constRef<bool>("isen")),
        isch (pdc->getKey_constRef<bool>("isch")),
        fEpsOut(pdc->getKey_constRef<delphi_real>("epsout")),
        fDebyeLength(pdc->getKey_constRef<delphi_real>("deblen")),
        fEpsIn(pdc->getKey_constRef<delphi_real>("epsin")),

        //ARGO-Putting in the reference for fgBoxCenter
        fgBoxCenter(pdc->getKey_constRef< SGrid<delphi_real> >("oldmid")),

        isite (pdc->getKey_constRef<bool>("isite")),
        ibem (pdc->getKey_constRef<bool>("ibem")),
        ibctyp (pdc->getKey_constRef<int>("ibctyp")),
        isitsf (pdc->getKey_constRef<bool>("isitsf")),
        bFrcAsPqr (pdc->getKey_constRef<bool>("frcpqr")),

        // Lin Li : Gaussian:
        cutoff (pdc->getKey_constRef<delphi_real>("cutoff")),
        sigma (pdc->getKey_constRef<delphi_real>("sigma")),
        inhomo (pdc->getKey_Ref<int>("inhomo")),
        srfcut (pdc->getKey_constRef<delphi_real>("srfcut")),
        dencut (pdc->getKey_constRef<delphi_real>("dencut")),
        iGaussian (pdc->getKey_constRef<int>("gaussian")),

        //iTestGloble (pdc->getKey_constRef<int>(" ")),
        iNMedia (pdc->getKey_constRef<delphi_integer>("nmedia")),
        numbmol (pdc->getKey_constRef<delphi_integer>("numbmol")),
        scrgfrm (pdc->getKey_constRef<delphi_integer>("scrgfrm")),

        global_cOldMid(pdc->getKey_Val< SGrid<delphi_real> >("oldmid")),
        fIonStrenth (pdc->getKey_constRef<delphi_real>("rionst")),
        fExternRadius (pdc->getKey_constRef<delphi_real>("exrad")),
        fRMax (pdc->getKey_constRef<delphi_real>("rdmx")),
        fRadPrb_v(pdc->getKey_Ref<vector <delphi_real> >("radprb")),

        dataobject_v(pdc->getKey_Ref< vector<string> >("dataobject")),

        //++++++++++++++++ reference to read-and-write variables from data container +++++++++++++++//
        //Read only
        delphipdb(pdc->getKey_Ref<vector <CAtomPdb> >("delphipdb")),
        rdmx(pdc->getKey_Ref<delphi_real>("rdmx")),
        iacs(pdc->getKey_Ref< bool >("iacs")),
        isrf(pdc->getKey_Ref< bool >("isrf")),
        cMin(pdc->getKey_Ref< SGrid<delphi_real> >("cmin")),
        cMax(pdc->getKey_Ref< SGrid<delphi_real> >("cmax")),

        //global=sum of local
        iBoundNum (pdc->getKey_Ref<delphi_integer>("ibnum")),
        nqass(pdc->getKey_Ref<delphi_integer>("nqass")),
        nqgrd(pdc->getKey_Ref<delphi_integer>("nqgrd")),

        //charges
        qnet(pdc->getKey_Ref< delphi_real >("qnet")),
        qmin(pdc->getKey_Ref< delphi_real >("qmin")),
        qplus(pdc->getKey_Ref< delphi_real >("qplus")),
        cqmin(pdc->getKey_Ref< SGrid<delphi_real> >("cqmin")),
        cqplus(pdc->getKey_Ref< SGrid<delphi_real> >("cqplus")),

        //Read only
        acenter(pdc->getKey_Ref< SGrid<delphi_real> >("acent")),
        //Read only
        medeps(pdc->getKey_Ref < vector < delphi_real > >("medeps")),

        xn1_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("xn1")), //Not used in other modules
        xn2_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("xn2")), //Not used in other modules

        //ARGO: doing the 'pdc' thing just like its done with idebmap
        zetaSurfMap_v( pdc->getKey_Ref< vector< char > >("zetaSurfMap")),
        zetaOn( pdc->getKey_Ref<int>("zetaOn")),
        zetaDistance( pdc->getKey_Ref<delphi_real>("zetaDistance")),
        strZetaPhiFile(pdc->getKey_constRef<string>("zphinam")),
        dMoment( pdc->getKey_Ref< vector<delphi_real> >("dMoment")),
        qMoment( pdc->getKey_Ref< vector<delphi_real> >("qMoment")),
        surf_grid_coords_v( pdc->getKey_Ref< vector<delphi_real> >("surf_grid_coords") ),
        surf_grid_index_v( pdc->getKey_Ref< vector<delphi_real> >("surf_grid_index") ),

        //all grids
        bDebMap_v(pdc->getKey_Ref< vector< char > >("idebmap")),
        iEpsMap_v( pdc->getKey_Ref< vector< SGrid<delphi_integer> > > ("iepsmp")),
        fGepsMap_v( pdc->getKey_Ref< vector< SGrid<delphi_real> > > ("gepsmp")),
        fGepsMap2_v( pdc->getKey_Ref< vector< SGrid<delphi_real> > > ("gepsmp2")),
        fGDensityMap_v(pdc->getKey_Ref< vector< delphi_real> >("gdensity")),

        //Read only
        iAtomMed_v(pdc->getKey_Ref< vector<delphi_integer> >("iatmmed")),
        sLimObject(pdc->getKey_Ref< vector < SExtrema<delphi_real> > >("limobject")),

        //******iBoundNum related********
        //Initial size = 0; Final size = iBoundNum
        ibgrd_v(pdc->getKey_Ref< vector< SGrid<delphi_integer> > >("ibgrd")),
        //Initial size = 0; Final size = iBoundNum, same as ibgrd_v, but it's real, used in energy, but why we need this?
        scspos_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("scspos")),
        //Size = iBoundNum, Atom number at surface, used in energy
        atsurf_v(pdc->getKey_Ref< vector< delphi_integer > >("atsurf")),
        //Not used in other modules, size = iBoundNum
        scsnor_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("scsnor")),
        //Not used in other modules, ignore. Size = iBoundNum, atom number at boundary,
        atndx_v(pdc->getKey_Ref< vector< delphi_integer > >("atndx")),

        //******nqass related*********
        //Size = nqass, stores coordinates of charged atoms
        chgpos_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("chgpos")),
        //Size = nqass
        crgatn_v(pdc->getKey_Ref< vector< delphi_integer > >("crgatn")),
        //Size = nqass, atom
        atmeps_v(pdc->getKey_Ref< vector< delphi_real > >("atmeps")),
        //Size = nqass
        atmcrg_v(pdc->getKey_Ref< vector< SGridValue<delphi_real> > >("atmcrg")),

        //############### For parallel Space module: ######################
#ifdef PARALLEL_MPI
        my_id(pdc->myid()), //if not MPI, define my_id=0
#else 
        my_id(0),           //if not MPI, define my_id=0
#endif // PARALLEL_MPI

        //******** nqgrd related******
        nqgrdtonqass_v(pdc->getKey_Ref< vector< delphi_integer > >("nqgrdtonqass")),
        //Size = nqgrd
        chrgv2_v(pdc->getKey_Ref< vector< SGridValue<delphi_real> > >("chrgv2"))
    {
        delphi_integer iGrid1D = pdc->getKey_Val<delphi_integer>("igrid");
        global_iGrid = { iGrid1D, iGrid1D, iGrid1D };
        global_iNatom= pdc->getKey_Val<delphi_integer>("natom");
        iGrid = { iGrid1D, iGrid1D, iGrid1D };
        cOldMid = global_cOldMid;

        epsdim = global_iNatom + iNObject + 2;

        iAtomMed=&iAtomMed_v[0]-1;

        ndistr = 0;

        test_pdc=pdc;

        // initialize all the pointers to be NULL:
        tmlst=NULL;
        egrid=NULL;
        idebmap=NULL;

        //ARGO
        zetaSurfMap=NULL;

        iepsmp=NULL;
        gepsmp=NULL;
        gepsmp2=NULL;
        gDensityMapOnGridPoint=NULL;

        iab1=NULL;
        iab2=NULL;
        cbn1=NULL;
        cbn2=NULL;
    };


    ~CDelphiSpace()
    {
        //delete [] sDelPhiPDB ;

#ifdef PARALLEL_MPI
        //if (fRadPrb_v.size()      > 0) fRadPrb_v.clear();
        if (dataobject_v.size()   > 0) vector<string>().swap(dataobject_v);                 //dataobject_v.clear();
        if (medeps.size()         > 0) vector<delphi_real>().swap(medeps);                  //medeps.clear();
        if (xn1_v.size()          > 0) vector< SGrid<delphi_real> >().swap(xn1_v);          //xn1_v.clear();
        if (xn2_v.size()          > 0) vector< SGrid<delphi_real> >().swap(xn2_v);          //xn2_v.clear();
        if (bDebMap_v.size()      > 0) vector<char>().swap(bDebMap_v);                      //bDebMap_v.clear();
        if (iEpsMap_v.size()      > 0) vector< SGrid<delphi_integer> >().swap(iEpsMap_v);   //iEpsMap_v.clear();
        if (fGepsMap_v.size()     > 0) vector< SGrid<delphi_real> >().swap(fGepsMap_v);     //fGepsMap_v.clear();
        if (fGepsMap2_v.size()    > 0) vector< SGrid<delphi_real> >().swap(fGepsMap2_v);    //fGepsMap2_v.clear();
        if (fGDensityMap_v.size() > 0) vector<delphi_real>().swap(fGDensityMap_v);          //fGDensityMap_v.clear();
        if (iAtomMed_v.size()     > 0) vector<delphi_integer>().swap(iAtomMed_v);           //iAtomMed_v.clear();
        if (sLimObject.size()     > 0) vector< SExtrema<delphi_real> >().swap(sLimObject);  //sLimObject.clear();  // pdc->getKey_Val
        if (ibgrd_v.size()        > 0) vector< SGrid<delphi_integer> >().swap(ibgrd_v);     //ibgrd_v.clear();
        if (scspos_v.size()       > 0) vector< SGrid<delphi_real> >().swap(scspos_v);       //scspos_v.clear();
        if (chgpos_v.size()       > 0) vector< SGrid <delphi_real> >().swap(chgpos_v);      //chgpos_v.clear();
        if (crgatn_v.size()       > 0) vector<delphi_integer>().swap(crgatn_v);             //crgatn_v.clear();
        if (nqgrdtonqass_v.size() > 0) vector<delphi_integer>().swap(nqgrdtonqass_v);       //nqgrdtonqass_v.clear();
        if (atmeps_v.size()       > 0) vector<delphi_real>().swap(atmeps_v);                //atmeps_v.clear();
        if (atmcrg_v.size()       > 0) vector< SGridValue<delphi_real> >().swap(atmcrg_v);  //atmcrg_v.clear();
        if (chrgv2_v.size()       > 0) vector< SGridValue<delphi_real> >().swap(chrgv2_v);  //chrgv2_v.clear();
        if (scsnor_v.size()       > 0) vector< SGrid<delphi_real> >().swap(scsnor_v);       //scsnor_v.clear();
        if (atsurf_v.size()       > 0) vector<delphi_integer>().swap(atsurf_v);             //atsurf_v.clear();
        if (atndx_v.size()        > 0) vector<delphi_integer>().swap(atndx_v);              //atndx_v.clear();
        if (delphipdb.size()      > 0) vector<CAtomPdb>().swap(delphipdb);                  //delphipdb.clear();  // pdc->getKey_Val
        if (zetaSurfMap_v.size()  > 0) vector<char>().swap(zetaSurfMap_v);                  //zetaSurfMap_v.clear();

        /*
         * semi global variables in this class
         */
        if (iab1_v.size()         > 0) vector<delphi_integer>().swap(iab1_v);
        if (iab2_v.size()         > 0) vector<delphi_integer>().swap(iab2_v);
        if (icume.size()          > 0) vector<delphi_integer>().swap(icume);
        if (ast.size()            > 0) vector<delphi_integer>().swap(ast);
        if (cbn1_v.size()         > 0) vector<delphi_integer>().swap(cbn1_v);
        if (cbn2_v.size()         > 0) vector<delphi_integer>().swap(cbn2_v);
        if (cbal.size()           > 0) vector<delphi_integer>().swap(cbal);
        if (icbn.size()           > 0) vector<delphi_integer>().swap(icbn);
        if (r0.size()             > 0) vector<delphi_real>().swap(r0);
        if (r02.size()            > 0) vector<delphi_real>().swap(r02);
        if (rs2.size()            > 0) vector<delphi_real>().swap(rs2);
        if (expos.size()          > 0) vector<SGrid <delphi_real>>().swap(expos);
        if (sLimGridUnit.size()   > 0) vector<SExtrema<delphi_real>>().swap(sLimGridUnit);
        if (vert.size()           > 0) vector<SGrid<delphi_real>>().swap(vert);
        if (vnorm.size()          > 0) vector<SGrid<delphi_real>>().swap(vnorm);
        if (vnorm2.size()         > 0) vector<SGrid<delphi_real>>().swap(vnorm2);
        if (vindx.size()          > 0) vector<SGrid<delphi_integer>>().swap(vindx);
        if (vtlen.size()          > 0) vector<delphi_integer>().swap(vtlen);
        if (vtlst.size()          > 0) vector<delphi_integer>().swap(vtlst);
        if (globalAtomIndex.size()> 0) vector<delphi_integer>().swap(globalAtomIndex);

        /*
         * For parallel Space module
         */
        if (vtpnt.size()          > 0) vector<delphi_integer>().swap(vtpnt);

#endif
    };

    virtual void run();
    virtual void validateInput();

    inline bool checkNotBufferGrid(SGrid<delphi_integer> localGrid)
    {
        return ( optANDGE<delphi_integer>(localGrid, myInternalStartNoBuffer) && optANDLE<delphi_integer>(localGrid, myInternalEndNoBuffer) );
    }

    inline bool checkNotBufferCoor(SGrid<delphi_real> localGrid)
    {
        return ( optANDGE<delphi_real>(localGrid, myStartCoorNoBuffer) && optANDLT<delphi_real>(localGrid, myEndCoorNoBuffer) );
    }

    static inline bool compareCoorX(SGrid<delphi_real> a, SGrid<delphi_real> b)
    {
        return a.nX < b.nX;
    }

    static inline bool compareCoorY(SGrid<delphi_real> a, SGrid<delphi_real> b)
    {
        return a.nY < b.nY;
    }

    static inline bool compareCoorZ(SGrid<delphi_real> a, SGrid<delphi_real> b)
    {
        return a.nZ < b.nZ;
    }


    inline SGrid<delphi_integer> convLocalGridToGlobal(SGrid<delphi_integer> localGrid)
            {
        return localGrid + myStart;
            }

    inline SGrid<delphi_real> convLocalGridToGlobal(SGrid<delphi_real> localGrid)
            {
        return localGrid + optCast<delphi_real, delphi_integer>(myStart);
            }

    inline SGrid<delphi_real> convLocalCoorToGlobal(SGrid<delphi_real> localCoor)
            {
        return localCoor + myStartCoor ;
            }

    inline delphi_integer convLocalAtomIndexToGlobal(delphi_integer localIndex)
    {
        return (localIndex < 1) ? 0 : globalAtomIndex[localIndex];
    }

#ifdef PARALLEL_MPI
    template <typename T> bool checkValue(T a, T b) { return (a == b); };

    template <typename T> bool checkValue(const string strKey)
    {
        cout << "Node "<< pdc->myid()<<"  comparing value " << strKey <<"...";

        T a = pdc->readGlobalVar<T>(strKey);

        T b = pdc->getKey_Val<T>(strKey);

        bool success = (a == b);
        if (success)
        {
            cout << "success Global=" << a << "  local=" << b << endl;
        }
        else
        {
            cout << "failed Global=" <<a<<"  local="<<b<< endl;
        }

        return success;
    };

    template <typename T> bool checkVector(vector<T> a, vector<T> b)
    {
        if (a.size() != b.size()) { return false; };
        typename vector<T>::iterator ia = a.begin();
        typename vector<T>::iterator ib=b.begin();
        for (; ia != a.end(); ia++, ib++)
        {
            if (*ia != *ib) { return false; }
        }
        return true;
    };

    template <typename T> bool checkVector(const string strKey)
    {
        cout << "Node " << pdc->myid() << "comparing vector " << strKey << " size: ";

        vector<T> &b = pdc->getKey_Ref<vector<T>>(strKey);

        vector<T> a(b.size());
        pdc->readGlobalVector1D<T>(strKey, 0, b.size(),a);

        bool success = true;

        for (typename vector<T>::iterator ia = a.begin(), ib = b.begin(); ia != a.end(); ia++, ib++)
        {
            if (*ia != *ib) { success = false; }
        }

        if (success)
        {
            cout << "success   ";
        }
        else
        {
            cout << "failed   ";
        }

        cout << "global=" << pdc->sizeofGlobal1D<T>(strKey)<< "   local=" << b.size() << endl;

        return success;
    }

    bool checkString(const string strKey)
    {
        string local = pdc->getKey_constRef<string>(strKey);
        string global = pdc->readGlobalChar(strKey);
        cout << "comparing string " << strKey;

        if(local== global)cout << "success  " ;
        else cout << "failed   " ;

        cout<<" global=" << global << "  local=" << local <<endl;
    }

    template <typename T> void sum_global_size(ddm::Shared<size_t> &global,
            dplt::global1D<delphi_integer> &size, vector<T> &vec)
    {
        size[pdc->myid()] = vec.size();
        if (pdc->myid() == 0)
        {
            delphi_integer total = 0;
            for (int i = 0; i < pdc->num_procs(); i++)
            {
                total += size[i];
            }
            global.set(total);
        }
    }

    template <typename T> void sum_global_value (ddm::Shared<size_t> &global,
            dplt::global1D<delphi_integer> &value, T &local)
    {
        value[pdc->myid()] = local;
        if (pdc->myid() == 0)
        {
            delphi_integer total = 0;
            for (int i = 0; i < pdc->num_procs(); i++)
            {
                total += value[i];
            }
            global.set(total);
        }
    }

    void write_to_global();
    void reallocate_global(delphi_integer usingCPU);
    void merge_global() {};

#else // NOT PARALLEL_MPI
    void write_to_global() {};
    void reallocate_global(delphi_integer usingCPU) {};
    void merge_global() {};
#endif //PARALLEL_MPI

    SGrid<delphi_integer> calcCPUDim();

    SGrid<delphi_integer> calcMyCPULocation();

    SGrid<delphi_integer> calcGridDim();

    SGrid<delphi_integer> calcLocalGridStartEnd();

    SGrid<delphi_integer> calcLocalGridStartEnd_fast();

    bool calcLocalGridStartEnd_naive();

    bool assignBuffer();

    void verify();

    void split();

    void mpi_run(){};

};

#endif // SPACE_H
