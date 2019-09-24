#ifndef SOLVER_FASTSOR_H
#define SOLVER_FASTSOR_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <memory>
#include <cmath>      // std::abs
//#include <deque>

#include "../interface/interface_abstractmodule.h"
#include "../misc/misc_timer.h"
#include "../delphi/delphi_constants.h"
#include "../io/io.h"
#include "solver_exceptions.h"

#ifdef PARALLEL_MPI
#include <mpi.h>
#endif

using namespace std;

class CDelphiFastSOR: virtual public IAbstractModule {
    
private:
    shared_ptr<CTimer> pTimer;

    /*
     * +++++++++++++++ references to the variables obtained from the data container +++++++++++++++
     *
     * const references to read-only variables from data container
     */
    
    //----- uniform parameters
    const string& strBioModel;                                          // biomodel
    const string& strNumSolver;                                         // solver
    
    //----- set by Statements
    const delphi_integer& iAtomNum;                                     // natom
    const bool& bCrgInterplateType;                                     // isph true-spherical 
                                                                        //      false-cubic
    const bool& bSpectralRadius;                                        // iuspec
    const delphi_real& fSpectralRadius;                                 // uspec
    vector<char>& rgbPeriodicBndy;                                      // iper
    const int& iNonIterateNum;                                          // nnit
    const delphi_real& fIonRadius;                                      // exrad
    const bool& bFixedRelaxParam;                                       // icheb
    const bool& bAutoConverge;                                          // iautocon
    const delphi_real& fGridConverge;                                   // gten
    const delphi_real& fRmsc;                                           // res1
    const delphi_real& fMaxc;                                           // res2
    const bool& bLogGraph;                                              // igraph
    const bool& bLogPotential;                                          // ipoten
    const bool& bManualRelaxParam;                                      // imanual
    const int& phiintype;                                               // for focusing
    const delphi_real& repsout;                                         // repsout
    const delphi_real& repsin;                                          // repsin
    
    //----- io file names
    const string& strEpsFile;                                           // epsnam
    const string& strPhiiFile;                                          // phiinam
    
    //----- set by functions
    const bool& bGridCrgOut;                                            // iwgcrg
    const bool& bEpsOut;                                                // epswrt
    const bool& bIonsEng;                                               // logions
    const int& iGaussian;                                               // Gaussian
    const int& inhomo;                                                  // inhomo

    //----- set by DelPhi
    const delphi_real& fEpsOut;                                         // epsout
    const delphi_real& fDebyeLength;                                    // deblen
    const delphi_real& fScale;                                          // scale
    const delphi_real& fEpsIn;                                          // epsin
    delphi_real fIonStrength;                                           // rionst
    const int& iDirectEpsMap;                                           // idirectalg
    const delphi_real& fEPKT;                                           // epkt
    const SGrid<delphi_real>& fgBoxCenter;                              // oldmid
    const SGrid<delphi_real>& fgAcenter;                                // acent
    const SGrid<delphi_real>& fgPotentialDrop;                          // vdrop
    const delphi_real& fTaylorCoeff2;                                   // chi2
    const delphi_real& fTaylorCoeff3;                                   // chi3
    const delphi_real& fTaylorCoeff4;                                   // chi4
    const delphi_real& fTaylorCoeff5;                                   // chi5
    const bool& uniformdiel;                                            //uniformdiel for Gaussian
    const delphi_real& fsrfcut;                                         //gaussian surface cutoff
    
    //----- set by IO class
    const delphi_integer& iMediaNum;                                    // nmedia
    const delphi_integer& iObjectNum;                                   // nobject
    vector<delphi_real>& prgfMediaEps;                                  // medeps(0:nmediamax)
    
    //----- set by Surface class
    const delphi_integer& iCrgGridNum;                                  // nqass
    const delphi_real& fMinusCrg;                                       // qmin
    const delphi_real& fPlusCrg;                                        // qplus
    const SGrid<delphi_real>& fgPlusCrgCenter;                          // cqplus
    const SGrid<delphi_real>& fgMinusCrgCenter;                         // cqmin
    const delphi_integer& iBndyGridNum;                                 // ibnum
    vector<SGrid<delphi_integer> >& prgigBndyGrid;                      // ibgrd(ibnum)
    const delphi_integer& iCrg2GridNum;                                 // nqgrd
    vector<SGridValue<delphi_real> >& prggvCrg2Grid;                    // chrgv2(nqgrd)
    vector<delphi_integer>& prgiCrg2GridMap;                            // nqgrdtonqass(nqgrd)
    vector<delphi_real>& prgfAtomEps;                                   // atmeps(nqass)
    vector<SGridValue<delphi_real> >& prggvCrgedAtom;                   // atmcrg(nqass)
    //const vector< SGrid<delphi_real> >& prgfgSurfCrgA;                // scspos(ibnum)

    vector<SGrid<delphi_integer> >& prgigEpsMap;                        // iepsmp(igrid,igrid,igrid)
    vector<SGrid<delphi_real> >& gepsmp2;                               // gepsmp2(igrid,igrid,igrid)
                                                                        // LinLi: Gaussian
    vector<char> prgbDielecMap;                                         // idebmap(igrid,igrid,igrid)
                                                                        // Zhe: modified for testing
    vector<delphi_real>& gaussianDensityMap;                            // gaussianDensityMap(igrid,igrid,igrid)

    /*
     * reference to read-and-write variables from data container
     */
    delphi_integer& iGrid;                                              // igrid (modified in setFocusBndy)
    int& iBndyType;                                                     // ibctyp (modified in setBndy)
    delphi_real& fNetCrg;                                               // qnet 
                                                                        // (modified in setCoulombBndy and setDipolarBndy)
    int& iLinIterateNum;                                                // nlit
    int& iIterateInterval;                                              // icon1
    int& iConvergeFract;                                                // icon2
    bool& bDbOut;                                                       // idbwrt
    delphi_real& fRelaxParam;                                           // relpar
    
    //----- out into data container
    delphi_integer& iDielecBndyOdd;                                     // icount2b
    delphi_integer& iCrgedGridSum;                                      // icount1b
    vector<delphi_real>& prgfGridCrg;                                   // gchrg(icount1b)
    vector<SGrid<delphi_integer> >& prgigGridCrgPose;                   // gchrgp(icount1b)
    delphi_integer& iCrgBndyGridNum;                                    // ibc
    vector<SDoubleGridValue>& prgdgvCrgBndyGrid;                        // cgbp(ibc)
    vector<delphi_real>& prgfPhiMap;                                    // phimap(igrid,igrid,igrid)
    vector<delphi_real>& phimap_pre_v;                                  // phimap_pre
    
    /*
     * +++++++++++++++ variables defined in this class +++++++++++++++
     *
     * const class variables
     */
    const delphi_integer iTotalGridNum;                                 // ngp=igrid*igrid*igrid+1
    const delphi_integer iHalfGridNum;                                  // nhgp=ngp/2
    const delphi_integer iBndyDielecNum;                                // nsp=2*(ibnum+1)
    const delphi_real fDebFct;                                          // debfct
    const delphi_real fEpsDiff;                                         // difeps
    const delphi_real fSixEps;                                          // sixeps
    const delphi_integer iEpsDim;                                       // epsdim
    const int nxran;                                                    // nxran
    const int nyran;                                                    // nyran
    bool debug_solver;
    delphi_real fsrfdens;                                               // Gaussian *sollution* density corresponding 
                                                                        // to the srfcut dielectricsonstant 
                                                                        // (denssolution = 1- densprotein)
    
    /*
     * local variables in this class
     *
     * using either std:vector or std:deque we can construct idpos, db etc. w/o realignment
     */
    delphi_integer iDielecBndyEven;                                     // icount2a
    vector<delphi_integer> prgiBndyDielecIndex;                         // idpos(nsp)
    vector<vector<delphi_real> > prgfBndyDielec;                        // db(6,nsp) 
                                                                        // <--> prgfBndyDielec[iBndyDielecNum][6]
    vector<delphi_real> prgfSaltMap1;                                   // sf1(nhgp)
    vector<delphi_real> prgfSaltMap2;                                   // sf2(nhgp)

    delphi_integer iCrgedGridEven;                                      // icount1a
    vector<delphi_integer> prgiCrgPose;                                 // iqpos(icount1b)
    vector<delphi_real> prgfCrgValA;                                    // qval(icount1b)
    vector<delphi_real> prgfCrgValG;                                    // gval(icount1b)

    //======Used in the Gaussian Salt nonliear Iterator
    vector<delphi_real> gaussianBoundaryDensity;                        // gaussianBoundaryDielec[nBound]
    vector<vector<delphi_real> > gaussianBoundaryDielec;                // gaussianBoundaryDielec[nBound][6]
    vector<delphi_real> gaussianChargeDensity;                          // gaussianBoundaryDielec[nCharge]
    vector<vector<delphi_real> > gaussianChargeDielec;                  // gaussianBoundaryDielec[nCharge][6]
    vector<delphi_real> gaussianBoundaryNonlinear;                      // gaussianBoundaryNonlinear1(nBound)
    vector<delphi_real> gaussianChargeNonlinear;                        // gaussianChargeNonlinear1(nCharge)

    delphi_real fSpec;                                                  // spec
    vector<delphi_integer> ibndx, ibndy, ibndz;
    delphi_integer idif1x, idif2x, inc1xa, inc1xb, inc2xa, inc2xb;
    delphi_integer idif1y, idif2y, inc1ya, inc1yb, inc2ya, inc2yb;
    delphi_integer idif1z, idif2z, inc1za, inc1zb, inc2za, inc2zb;
    vector<delphi_integer> sta1, sta2, fi1, fi2;
    delphi_integer lat1, lat2, long1, long2;
    vector<delphi_real> phimap1, phimap2;
    vector<delphi_real> bndx1, bndx2, bndx3, bndx4;
    vector<delphi_real> qmap1, qmap2;
    vector<delphi_real> debmap1, debmap2;
    delphi_real om1, om2, om3, om4, sixth;

    const char * infoString = " Info> ";
    const char * timeString = " Time> ";
    const char * bndcString = " Bndc> ";
    size_t MAXWIDTH;

#ifdef MCCE
    SMCCE* pmcce;
#endif

#ifdef PRIME
    SPrime* pPrime;
#endif

    void setBndy();                             // subroutine setbc()

    //+++++ choices of boundary conditions
    bool isDipolarBndy(delphi_real *** phimap); // ibctyp = 2
    bool isFocusBndy(delphi_real *** phimap);   // ibctyp = 3
    bool isCoulombBndy(delphi_real *** phimap); // ibctyp = 4
 
    delphi_real calculateRelaxFactor();         // subroutine relfac()

    void conplt(const vector<delphi_real>& array, const string& title, const int& iclr, const int& iscl, 
                const int& imk, const int& iplt, const char symb, const int& ixmin, const int& ixmax, 
                vector<string>& iplot, delphi_real& ymin, delphi_real& ymax);
    void postItr(const vector<delphi_real>& rmaxl, const vector<delphi_real>& rmsl);

    //shared_ptr<IDataContainer> solver_pdc;

    delphi_real calcExpSolvE(delphi_real);
    delphi_real calcSinh(delphi_real);
    delphi_real calcPhiMinusSinh(delphi_real);
    delphi_real calcExp(delphi_real);

    void setDielecBndySaltMap();                // subroutine mkdbsf()
    void setCrg();                              // subroutine setcrg()

    void itit();                                // subroutine itit()
    void nitit(const delphi_real& qfact);       // subroutine nitit(qfact)
    void initOddEvenItr(const int& forWhom);
    void itrOddPoints(const int& forWhom, const int&);
    void itrEvenPoints(const int& forWhom, const int&);

#ifdef PARALLEL_MPI
    /*
     * +++++++++++++++ local variables for MPI usage +++++++++++++++
     */

    /*
     * mpi variables for general information
     */
    int mpi_num_procs;                              // total # of processes
    int mpi_num_workers;                            // total # of slave processes = total # of processes - 1
    int mpi_rank;                                   // process rank or id

    /*
     * mpi variables for MPI_Send and MPI_Recv
     */
    int mpi_sender;                                 // sender id
    int mpi_sendtag;                                // tag associated with one sent message
    int mpi_recver;                                 // receiver id
    int mpi_recvtag;                                // tag associated with one received message
    MPI_Status mpi_stat;                            // mpi status

    delphi_real *mpi_phimap1 = NULL, *mpi_phimap2 = NULL; // pointers to phimap1.data() and phimap2.data()
    delphi_real *mpi_sf1     = NULL, *mpi_sf2     = NULL; // pointers to prgfSaltMap1(sf1) and prgfSaltMap1(sf2)
    delphi_real *mpi_qmap1   = NULL, *mpi_qmap2   = NULL;
    delphi_real *mpi_debmap1 = NULL, *mpi_debmap2 = NULL;

    /*
     * offsets in phimap1-2 to the beginning of the vector
     */ 
    delphi_integer mpi_wrphimap1_start, mpi_wrphimap1_end; // not used: mpi_wrphimap1_offset
    delphi_integer mpi_wrphimap2_start, mpi_wrphimap2_end; // not used: mpi_wrphimap2_offset
    delphi_integer mpi_wrqmap1_start, mpi_wrqmap2_start;
    delphi_integer mpi_wrdebmap1_start, mpi_wrdebmap2_start;
    
    /*
     * mpi variables for constructing phimap1 and phimap2
     */
    delphi_integer *mpi_wrstar1 = NULL, *mpi_wrstar2 = NULL;      // mpi_wrstar1[mpi_num_procs], sf1(WRstar1:WRfinl1)
                                                                  // mpi_wrstar2[mpi_num_procs], sf2(WRstar2:WRfinl2)
    delphi_integer *mpi_wrfinl1 = NULL, *mpi_wrfinl2 = NULL;      // mpi_wrfinl1[mpi_num_procs], mpi_wrfinl2[mpi_num_procs]
    delphi_integer *mpi_wrlen1 = NULL, *mpi_wrlen2 = NULL;        // mpi_wrlen1[mpi_num_procs], mpi_wrlen2[mpi_num_procs]
    delphi_integer *mpi_wrnstafi1 = NULL, *mpi_wrnstafi2 = NULL;  // mpi_wrnstafi1[mpi_num_procs], size of local sta1 and fi1

    delphi_integer mpi_mrlen1, mpi_mrlen2, mpi_wrsplit1, mpi_wrsplit2;

    /*
     * mpi variables for local charged positions
     */
    delphi_integer mpi_wricount1a, mpi_wricount1b, mpi_wricount2a, mpi_wricount2b;

    /*
     * used for MPI_Scatterv and MPI_Gatherv
     */
    delphi_integer *mpi_recvcounts1  = NULL, *mpi_recvdispls1  = NULL;
    delphi_integer *mpi_sendcounts1l = NULL, *mpi_senddispls1l = NULL;
    delphi_integer *mpi_sendcounts1r = NULL, *mpi_senddispls1r = NULL;
    delphi_integer *mpi_recvcounts2  = NULL, *mpi_recvdispls2  = NULL; 
    delphi_integer *mpi_sendcounts2l = NULL, *mpi_senddispls2l = NULL;
    delphi_integer *mpi_sendcounts2r = NULL, *mpi_senddispls2r = NULL;

    /*
     * used for MPI_Sendrecv
     */
    delphi_integer mpi_ltostart1, mpi_ltosize1, mpi_lfromstart1, mpi_lfromsize1;
    delphi_integer mpi_rtostart1, mpi_rtosize1, mpi_rfromstart1, mpi_rfromsize1;
    delphi_integer mpi_ltostart2, mpi_ltosize2, mpi_lfromstart2, mpi_lfromsize2;
    delphi_integer mpi_rtostart2, mpi_rtosize2, mpi_rfromstart2, mpi_rfromsize2;

    /*
     * used for MPI_Win_Create
     */
    MPI_Aint mpi_zerodispl, mpi_todispl1, mpi_fromdispl1, mpi_todispl2, mpi_fromdispl2;
    MPI_Group mpi_wholegroup, mpi_startgroup, mpi_postgroup;
    int mpi_startrank, mpi_postrank, mpi_sizeofdouble;

    /*
     * mpi variables for one-sided memory access
     */
    MPI_Win mpi_towin1, mpi_towin2; // not used: mpi_fromwin1, mpi_fromwin2;
    delphi_integer mpi_towin1start, mpi_towin2start;
    MPI_Aint mpi_towin1size, mpi_towin2size;
    
    /*
     * mpi version of private functions: itit, nitit, initOddEvenItr, itrEvenPoints, itrOddPoints
     */
    void mpi_setDielecBndySaltMap();
    void mpi_setCrg();

    int  mpi_msgtag(const int& mpi_sender, int& mpi_recver) { return mpi_sender * 100 + mpi_recver;}
    void mpi_itit();
    void mpi_nitit(const delphi_real& qfact);
    void mpi_initOddEvenItr(const int& forWhom);
    void mpi_itrEvenPoints(const int& forWhom, const int& flag);
    void mpi_itrOddPoints(const int& forWhom, const int& flag);
#endif    
    
public:

    /**
     *  constructor
     */
    CDelphiFastSOR(shared_ptr<IDataContainer> pdc, shared_ptr<CTimer> pt) :
    /*
     * ++++++++++++++ references to the variables obtained from the data container +++++++++++++++
     *
     * const references to read-only variables from data container
     */
    IAbstractModule(pdc), 
    pTimer(pt),
    
    //----- uniform parameters
    strBioModel(pdc->getKey_constRef < string > ("biomodel")), 
    strNumSolver(pdc->getKey_constRef < string > ("solver")),
    
    //----- set by Statements
    iAtomNum(pdc->getKey_constRef < delphi_integer > ("natom")), 
    bCrgInterplateType(pdc->getKey_constRef<bool>("isph")), 
    bSpectralRadius(pdc->getKey_constRef<bool>("iuspec")), 
    fSpectralRadius(pdc->getKey_constRef < delphi_real > ("uspec")), 
    rgbPeriodicBndy(pdc->getKey_Ref < vector<char> > ("iper")), 
    iNonIterateNum(pdc->getKey_constRef<int>("nnit")), 
    fIonRadius(pdc->getKey_constRef < delphi_real > ("exrad")), 
    bFixedRelaxParam(pdc->getKey_constRef<bool>("icheb")), 
    bAutoConverge(pdc->getKey_constRef<bool>("iautocon")), 
    fGridConverge(pdc->getKey_constRef < delphi_real > ("gten")), 
    fRmsc(pdc->getKey_constRef < delphi_real > ("res1")), 
    fMaxc(pdc->getKey_constRef < delphi_real > ("res2")), 
    bLogGraph(pdc->getKey_constRef<bool>("igraph")), 
    bLogPotential(pdc->getKey_constRef<bool>("ipoten")), 
    bManualRelaxParam(pdc->getKey_constRef<bool>("imanual")), 
    phiintype(pdc->getKey_constRef<int>("phiintype")), 
    repsout(pdc->getKey_constRef < delphi_real > ("repsout")), 
    repsin(pdc->getKey_constRef < delphi_real > ("repsin")),
    
    //----- io file names
    strEpsFile(pdc->getKey_constRef < string > ("epsnam")), 
    strPhiiFile(pdc->getKey_constRef < string > ("phiinam")),
    
    //----- set by functions
    bDbOut(pdc->getKey_Ref<bool>("idbwrt")), 
    bGridCrgOut(pdc->getKey_constRef<bool>("iwgcrg")), 
    bEpsOut(pdc->getKey_constRef<bool>("epswrt")), 
    bIonsEng(pdc->getKey_constRef<bool>("logions")), 
    iGaussian(pdc->getKey_constRef<int>("gaussian")), 
    inhomo(pdc->getKey_constRef<int>("inhomo")), 
    fsrfcut(pdc->getKey_constRef < delphi_real > ("srfcut")),

    //---------- set for Gaussian Salt by space class
    gaussianDensityMap(pdc->getKey_Ref < vector < delphi_real >> ("gdensity")),

    //----- set by DelPhi
    fEpsOut(pdc->getKey_constRef < delphi_real > ("epsout")), 
    fDebyeLength(pdc->getKey_constRef < delphi_real > ("deblen")), 
    fScale(pdc->getKey_constRef < delphi_real > ("scale")), 
    fEpsIn(pdc->getKey_constRef < delphi_real > ("epsin")), 
    fIonStrength(pdc->getKey_Val < delphi_real > ("rionst")), 
    iDirectEpsMap(pdc->getKey_constRef<int>("idirectalg")), 
    fEPKT(pdc->getKey_constRef < delphi_real > ("epkt")), 
    fgBoxCenter(pdc->getKey_constRef < SGrid<delphi_real> > ("oldmid")), 
    fgAcenter(pdc->getKey_constRef < SGrid<delphi_real> > ("acent")), 
    fgPotentialDrop(pdc->getKey_constRef < SGrid<delphi_real> > ("vdrop")), 
    fTaylorCoeff2(pdc->getKey_constRef < delphi_real > ("chi2")), 
    fTaylorCoeff3(pdc->getKey_constRef < delphi_real > ("chi3")), 
    fTaylorCoeff4(pdc->getKey_constRef < delphi_real > ("chi4")), 
    fTaylorCoeff5(pdc->getKey_constRef < delphi_real > ("chi5")), 
    uniformdiel(pdc->getKey_constRef<bool>("uniformdiel")),
    
    //----- set by IO class
    iMediaNum(pdc->getKey_constRef < delphi_integer > ("nmedia")), 
    iObjectNum(pdc->getKey_constRef < delphi_integer > ("nobject")), 
    prgfMediaEps(pdc->getKey_Ref < vector<delphi_real> > ("medeps")),
    
    //----- set by Surface class
    iCrgGridNum(pdc->getKey_constRef < delphi_integer > ("nqass")), 
    fMinusCrg(pdc->getKey_constRef < delphi_real > ("qmin")), 
    fPlusCrg(pdc->getKey_constRef < delphi_real > ("qplus")), 
    fgPlusCrgCenter(pdc->getKey_constRef < SGrid<delphi_real> > ("cqplus")), 
    fgMinusCrgCenter(pdc->getKey_constRef < SGrid<delphi_real> > ("cqmin")), 
    iBndyGridNum(pdc->getKey_constRef < delphi_integer > ("ibnum")), 
    prgigBndyGrid(pdc->getKey_Ref < vector<SGrid<delphi_integer> > > ("ibgrd")), 
    gepsmp2(pdc->getKey_Ref < vector<SGrid<delphi_real> >> ("gepsmp2")),
    prgbDielecMap(pdc->getKey_Ref< vector<char> >("idebmap")),
    iCrg2GridNum(pdc->getKey_constRef < delphi_integer > ("nqgrd")), 
    prggvCrg2Grid(pdc->getKey_Ref < vector<SGridValue<delphi_real> > > ("chrgv2")), 
    prgiCrg2GridMap(pdc->getKey_Ref < vector<delphi_integer> > ("nqgrdtonqass")), 
    prgfAtomEps(pdc->getKey_Ref < vector<delphi_real> > ("atmeps")), 
    prggvCrgedAtom(pdc->getKey_Ref < vector<SGridValue<delphi_real> > > ("atmcrg")),
    //prgfgSurfCrgA(pdc->getKey_constRef< vector< SGrid<delphi_real> > >("scspos")),

    prgigEpsMap(pdc->getKey_Ref < vector<SGrid<delphi_integer> > > ("iepsmp")),

    /*
     * reference to read-and-write variables from data container
     */
    iGrid(pdc->getKey_Ref < delphi_integer > ("igrid")),                             // modified in setFocusBndy
    iBndyType(pdc->getKey_Ref<int>("ibctyp")),                                       // modified in setBndy
    fNetCrg(pdc->getKey_Ref < delphi_real > ("qnet")),                               // modified in setCoulombBndy 
                                                                                     // and setDipolarBndy
    iLinIterateNum(pdc->getKey_Ref<int>("nlit")), 
    iIterateInterval(pdc->getKey_Ref<int>("icon1")), 
    iConvergeFract(pdc->getKey_Ref<int>("icon2")), 
    fRelaxParam(pdc->getKey_Ref < delphi_real > ("relpar")),
    
    //----- out into data container
    iDielecBndyOdd(pdc->getKey_Ref < delphi_integer > ("icount2b")), 
    iCrgedGridSum(pdc->getKey_Ref < delphi_integer > ("icount1b")),                  // modified in setcrg
    iCrgBndyGridNum(pdc->getKey_Ref < delphi_integer > ("ibc")),                     // modified in setcrg
    
    prgfGridCrg(pdc->getKey_Ref < vector<delphi_real> > ("gchrg")),                  // modified in setcrg
    prgigGridCrgPose(pdc->getKey_Ref < vector<SGrid<delphi_integer> > > ("gchrgp")), // modified in setcrg
    prgdgvCrgBndyGrid(pdc->getKey_Ref < vector<SDoubleGridValue> > ("cgbp")),        // modified in setcrg
    prgfPhiMap(pdc->getKey_Ref < vector<delphi_real> > ("phimap")), 
    phimap_pre_v(pdc->getKey_Ref < vector<delphi_real> > ("phimap_pre")),

    /*
     * +++++++++++++++ variables defined in this class +++++++++++++++
     *
     * const class variables
     */
    iTotalGridNum(iGrid * iGrid * iGrid + 1), 
    iHalfGridNum(iTotalGridNum / 2), 
    iBndyDielecNum(2 * (iBndyGridNum + 1)), 
    fDebFct(fEpsOut / (fDebyeLength * fScale * fDebyeLength * fScale)), 
    fEpsDiff(fEpsIn - fEpsOut), 
    fSixEps(fEpsOut * 6.0), 
    iEpsDim(iAtomNum + iObjectNum + 2), 
    nxran(60), 
    nyran(20) 
    {        
        if (0 != strBioModel.compare("PBE") && 0 != strNumSolver.compare("DELPHI")) 
            throw CUnknownBioModelSolver(strBioModel, strNumSolver);

        //++++++++++++++++++++++++++++ local variables in this class ++++++++++++++++++++++++++++//
        iDielecBndyEven = 0;
        iCrgedGridEven  = 0;
        fSpec           = 0.0;
        fsrfdens        = (fsrfcut - repsin) / (repsout - repsin);
        //solver_pdc    = pdc;
        
#ifdef PARALLEL_MPI             
        mpi_num_procs    = pdc->num_procs();
        mpi_num_workers  = mpi_num_procs - 1;
        mpi_rank         = pdc->myid(); 
        
        mpi_wrstar1      = new delphi_integer[mpi_num_procs];
        mpi_wrstar2      = new delphi_integer[mpi_num_procs];
        mpi_wrfinl1      = new delphi_integer[mpi_num_procs];
        mpi_wrfinl2      = new delphi_integer[mpi_num_procs];
        mpi_wrlen1       = new delphi_integer[mpi_num_procs];
        mpi_wrlen2       = new delphi_integer[mpi_num_procs];
        mpi_wrnstafi1    = new delphi_integer[mpi_num_procs]; mpi_wrnstafi1[0] = iGrid - 1;
        mpi_wrnstafi2    = new delphi_integer[mpi_num_procs]; mpi_wrnstafi2[0] = iGrid - 1;
        
        /*
         * for MPI_GATHERV use
         */
        mpi_recvcounts1 = new delphi_integer[mpi_num_procs]; mpi_recvcounts1[0] = 0;
        mpi_recvdispls1 = new delphi_integer[mpi_num_procs]; mpi_recvdispls1[0] = 0;
        mpi_recvcounts2 = new delphi_integer[mpi_num_procs]; mpi_recvcounts2[0] = 0;
        mpi_recvdispls2 = new delphi_integer[mpi_num_procs]; mpi_recvdispls2[0] = 0;

        /*
         * for MPI_SCATTERV use
         */
        mpi_sendcounts1l = new delphi_integer[mpi_num_procs]; mpi_sendcounts1l[0] = 0;
        mpi_senddispls1l = new delphi_integer[mpi_num_procs]; mpi_senddispls1l[0] = 0;
        mpi_sendcounts1r = new delphi_integer[mpi_num_procs]; mpi_sendcounts1r[0] = 0;
        mpi_senddispls1r = new delphi_integer[mpi_num_procs]; mpi_senddispls1r[0] = 0;

        mpi_sendcounts2l = new delphi_integer[mpi_num_procs]; mpi_sendcounts2l[0] = 0;
        mpi_senddispls2l = new delphi_integer[mpi_num_procs]; mpi_senddispls2l[0] = 0;
        mpi_sendcounts2r = new delphi_integer[mpi_num_procs]; mpi_sendcounts2r[0] = 0;
        mpi_senddispls2r = new delphi_integer[mpi_num_procs]; mpi_senddispls2r[0] = 0;

        mpi_wricount1a = 0; mpi_wricount1b = 0;
        mpi_wricount2a = 0; mpi_wricount2b = 0;

        mpi_ltostart2 = 0; mpi_ltosize2 = 0; mpi_rtostart2 = 0; mpi_rtosize2 = 0;
        mpi_ltostart1 = 0; mpi_ltosize1 = 0; mpi_rtostart1 = 0; mpi_rtosize1 = 0;

        mpi_lfromstart2 = 0; mpi_lfromsize2 = 0; mpi_rfromstart2 = 0; mpi_rfromsize2 = 0;
        mpi_lfromstart1 = 0; mpi_lfromsize1 = 0; mpi_rfromstart1 = 0; mpi_rfromsize1 = 0;

        mpi_zerodispl = 0; mpi_todispl1 = 0; mpi_fromdispl1 = 0; mpi_todispl2 = 0; mpi_fromdispl2 = 0;                
#endif         
    };

    ~CDelphiFastSOR() 
    {
        #ifdef DEBUG_OBJECT
        cout << endl;
        cout << "****************************************************************\n";
        cout << "*                CDelphiFastSOR is destroyed                   *\n";
        cout << "****************************************************************\n";
        #endif

#ifdef PARALLEL_MPI
        delete[] mpi_wrstar1;     mpi_wrstar1       = NULL;
        delete[] mpi_wrstar2;     mpi_wrstar2       = NULL;
        delete[] mpi_wrfinl1;     mpi_wrfinl1       = NULL;
        delete[] mpi_wrfinl2;     mpi_wrfinl2       = NULL;
        delete[] mpi_wrlen1;      mpi_wrlen1        = NULL;
        delete[] mpi_wrlen2;      mpi_wrlen2        = NULL;
        delete[] mpi_wrnstafi1;   mpi_wrnstafi1     = NULL;
        delete[] mpi_wrnstafi2;   mpi_wrnstafi2     = NULL;

        delete[] mpi_recvcounts1; mpi_recvcounts1   = NULL;
        delete[] mpi_recvdispls1; mpi_recvdispls1   = NULL;
        delete[] mpi_recvcounts2; mpi_recvcounts2   = NULL;
        delete[] mpi_recvdispls2; mpi_recvdispls2   = NULL;

        delete[] mpi_sendcounts1l; mpi_sendcounts1l = NULL;
        delete[] mpi_senddispls1l; mpi_senddispls1l = NULL;
        delete[] mpi_sendcounts1r; mpi_sendcounts1r = NULL;
        delete[] mpi_senddispls1r; mpi_senddispls1r = NULL;

        delete[] mpi_sendcounts2l; mpi_sendcounts2l = NULL;
        delete[] mpi_senddispls2l; mpi_senddispls2l = NULL;
        delete[] mpi_sendcounts2r; mpi_sendcounts2r = NULL;
        delete[] mpi_senddispls2r; mpi_senddispls2r = NULL;
        
        mpi_phimap1   = NULL;
        mpi_phimap2   = NULL;
        mpi_sf1       = NULL; 
        mpi_sf2       = NULL;
        mpi_qmap1     = NULL; 
        mpi_qmap2     = NULL;
        mpi_debmap1   = NULL; 
        mpi_debmap2   = NULL;

        /*
         * clear allocated memory in local data containers
         */
        if (prgfMediaEps.size()       > 0)
            vector<delphi_real>().swap(prgfMediaEps);               //prgfMediaEps.clear();       // not sync'd
        if (prgigBndyGrid.size()      > 0)
            vector<SGrid<delphi_integer>>().swap(prgigBndyGrid);    //prgigBndyGrid.clear();      // not sync'd
        if (prggvCrg2Grid.size()      > 0)
            vector<SGridValue<delphi_real>>().swap(prggvCrg2Grid);  //prggvCrg2Grid.clear();      // not sync'd
        if (prgiCrg2GridMap.size()    > 0)
            vector<delphi_integer>().swap(prgiCrg2GridMap);         //prgiCrg2GridMap.clear();    // not sync'd
        if (prgfAtomEps.size()        > 0)
            vector<delphi_real>().swap(prgfAtomEps);                //prgfAtomEps.clear();        // not sync'd
        if (prggvCrgedAtom.size()     > 0)
            vector<SGridValue<delphi_real>>().swap(prggvCrgedAtom); //prggvCrgedAtom.clear();     // master::mpi_run()
        if (prgigEpsMap.size()        > 0)
            vector<SGrid<delphi_integer>>().swap(prgigEpsMap);      //prgigEpsMap.clear();        // not sync'd
        if (gepsmp2.size()            > 0)
            vector<SGrid<delphi_real>>().swap(gepsmp2);             //gepsmp2.clear();            // not sync'd
        if (prgbDielecMap.size()      > 0)
            vector<char>().swap(prgbDielecMap);                     //prgbDielecMap.clear();      // master::mpi_run()
        if (gaussianDensityMap.size() > 0)
            vector<delphi_real>().swap(gaussianDensityMap);         //gaussianDensityMap.clear(); // not sync'd
        if (prgfGridCrg.size()        > 0)
            vector<delphi_real>().swap(prgfGridCrg);                //prgfGridCrg.clear();        // master::mpi_setcrg()
        if (prgigGridCrgPose.size()   > 0)
            vector<SGrid<delphi_integer>>().swap(prgigGridCrgPose); //prgigGridCrgPose.clear();   // master::mpi_setcrg()
        if (prgdgvCrgBndyGrid.size()  > 0)
            vector<SDoubleGridValue>().swap(prgdgvCrgBndyGrid);     //prgdgvCrgBndyGrid.clear();  // master::mpi_setcrg()
        if (prgfPhiMap.size()         > 0)
            vector<delphi_real>().swap(prgfPhiMap);                 //prgfPhiMap.clear();         // master
        if (phimap_pre_v.size()       > 0)
            vector<delphi_real>().swap(phimap_pre_v);               //phimap_pre_v.clear();       // master
#endif    

        /*
         * free local vectors allocated in this class
         */
        if (prgiBndyDielecIndex.size()       > 0) vector<delphi_integer>().swap(prgiBndyDielecIndex);
        if (prgfBndyDielec.size()            > 0) vector<vector<delphi_real>>().swap(prgfBndyDielec);
        if (prgfSaltMap1.size()              > 0) vector<delphi_real>().swap(prgfSaltMap1);
        if (prgfSaltMap2.size()              > 0) vector<delphi_real>().swap(prgfSaltMap2);
        if (prgiCrgPose.size()               > 0) vector<delphi_integer>().swap(prgiCrgPose);
        if (prgfCrgValA.size()               > 0) vector<delphi_real>().swap(prgfCrgValA);
        if (prgfCrgValG.size()               > 0) vector<delphi_real>().swap(prgfCrgValG);
        if (gaussianBoundaryDensity.size()   > 0) vector<delphi_real>().swap(gaussianBoundaryDensity);
        if (gaussianBoundaryDielec.size()    > 0) vector<vector<delphi_real>>().swap(gaussianBoundaryDielec);
        if (gaussianChargeDensity.size()     > 0) vector<delphi_real>().swap(gaussianChargeDensity);
        if (gaussianChargeDielec.size()      > 0) vector<vector<delphi_real>>().swap(gaussianChargeDielec);
        if (gaussianBoundaryNonlinear.size() > 0) vector<delphi_real>().swap(gaussianBoundaryNonlinear);
        if (gaussianChargeNonlinear.size()   > 0) vector<delphi_real>().swap(gaussianChargeNonlinear);
        if (ibndx.size()                     > 0) vector<delphi_integer>().swap(ibndx);
        if (ibndy.size()                     > 0) vector<delphi_integer>().swap(ibndy);
        if (ibndz.size()                     > 0) vector<delphi_integer>().swap(ibndz);
        if (sta1.size()                      > 0) vector<delphi_integer>().swap(sta1);
        if (sta2.size()                      > 0) vector<delphi_integer>().swap(sta2);
        if (fi1.size()                       > 0) vector<delphi_integer>().swap(fi1);
        if (fi2.size()                       > 0) vector<delphi_integer>().swap(fi2);
        if (phimap1.size()                   > 0) vector<delphi_real>().swap(phimap1);
        if (phimap2.size()                   > 0) vector<delphi_real>().swap(phimap2);
        if (bndx1.size()                     > 0) vector<delphi_real>().swap(bndx1);
        if (bndx2.size()                     > 0) vector<delphi_real>().swap(bndx2);
        if (bndx3.size()                     > 0) vector<delphi_real>().swap(bndx3);
        if (bndx4.size()                     > 0) vector<delphi_real>().swap(bndx4);
        if (qmap1.size()                     > 0) vector<delphi_real>().swap(qmap1);
        if (qmap2.size()                     > 0) vector<delphi_real>().swap(qmap2);
        if (debmap1.size()                   > 0) vector<delphi_real>().swap(debmap1);
        if (debmap2.size()                   > 0) vector<delphi_real>().swap(debmap2);
    };

    virtual void validateInput();

#ifdef MCCE
    void getMCCE(SMCCE* mcce_data) { pmcce = mcce_data; }
#endif

#ifdef PRIME
    void getPRIME(shared_ptr<SPrime> param) { pPrime = &*param; }
#endif

    virtual void run();

#ifdef PARALLEL_MPI
    /*
     * implementation of pure abstract function mpirun() in class IAbstractModule
     */
    void reallocate_global(delphi_integer usingCPU)
    {

        ddm::Shared<size_t> globalSize;
        dplt::globalVec1D<delphi_integer> distSize(pdc->num_procs());
        ddm::barrier();

        /*
         * prgfGridCrg(pdc->getKey_Ref< vector<delphi_real> >("gchrg")) - modified in setcrg
         */
        pdc->deallocateGlobal1D<delphi_real>("gchrg"); 
        if (pdc->myid() == usingCPU) globalSize.set(prgfGridCrg.size());
        ddm::barrier();
        pdc->allocateGlobal1D<delphi_real>("gchrg", globalSize.get());
        ddm::barrier();

        /*
         * prgigGridCrgPose(pdc->getKey_Ref< vector< SGrid<delphi_integer> > >("gchrgp")) - modified in setcrg
         */
        pdc->deallocateGlobal1D<SGrid<delphi_integer>>("gchrgp");
        if (pdc->myid() == usingCPU) globalSize.set(prgigGridCrgPose.size());
        ddm::barrier();
        pdc->allocateGlobal1D<SGrid<delphi_integer>>("gchrgp", globalSize.get());
        ddm::barrier();

        /*
         * prgdgvCrgBndyGrid(pdc->getKey_Ref< vector<SDoubleGridValue> >("cgbp")) - modified in setcrg 
         */
        pdc->deallocateGlobal1D<SDoubleGridValue>("cgbp");
        if (pdc->myid() == usingCPU) globalSize.set(prgdgvCrgBndyGrid.size());
        ddm::barrier();
        pdc->allocateGlobal1D<SDoubleGridValue>("cgbp", globalSize.get());
        ddm::barrier();

        /*
         * prgfPhiMap(pdc->getKey_Ref< vector<delphi_real> >("phimap")) 
         */
        pdc->deallocateGlobal1D<delphi_real>("phimap");
        if (pdc->myid() == usingCPU) globalSize.set(prgfPhiMap.size());
        ddm::barrier();
        pdc->allocateGlobal1D<delphi_real>("phimap", globalSize.get());
        ddm::barrier();

        /*
         * phimap_pre_v(pdc->getKey_Ref< vector<delphi_real> >("phimap_pre")) 
         */
        pdc->deallocateGlobal1D<delphi_real>("phimap_pre");
        if (pdc->myid() == usingCPU) globalSize.set(phimap_pre_v.size());
        ddm::barrier();
        pdc->allocateGlobal1D<delphi_real>("phimap_pre", globalSize.get());
        ddm::barrier();
    };

    void write_to_global() 
    {
        // iDielecBndyOdd(pdc->getKey_Ref<delphi_integer>("icount2b"))
        pdc->writeGlobalVar<delphi_integer>("icount2b", iDielecBndyOdd);
        
        // iCrgedGridSum(pdc->getKey_Ref<delphi_integer>("icount1b")) - modified in setcrg
        pdc->writeGlobalVar<delphi_integer>("icount1b", iCrgedGridSum);                   
        
        // iCrgBndyGridNum(pdc->getKey_Ref<delphi_integer>("ibc")) - modified in setcrg
        pdc->writeGlobalVar<delphi_integer>("ibc", iCrgBndyGridNum);                    
        
        // prgfGridCrg(pdc->getKey_Ref< vector<delphi_real> >("gchrg")) - modified in setcrg
        pdc->writeGlobalVector1D<delphi_real>("gchrg", 0, prgfGridCrg.size(), prgfGridCrg); 
        
        // prgigGridCrgPose(pdc->getKey_Ref< vector< SGrid<delphi_integer> > >("gchrgp")) - modified in setcrg
        pdc->writeGlobalVector1D<SGrid<delphi_integer>>("gchrgp", 0, prgigGridCrgPose.size(), prgigGridCrgPose); 
        
        // prgdgvCrgBndyGrid(pdc->getKey_Ref< vector<SDoubleGridValue> >("cgbp")) - modified in setcrg 
        pdc->writeGlobalVector1D<SDoubleGridValue>("cgbp", 0, prgdgvCrgBndyGrid.size(), prgdgvCrgBndyGrid);
        
        // prgfPhiMap(pdc->getKey_Ref< vector<delphi_real> >("phimap"))
        pdc->writeGlobalVector1D<delphi_real>("phimap", 0, prgfPhiMap.size(), prgfPhiMap);                  
        
        // phimap_pre_v(pdc->getKey_Ref< vector<delphi_real> >("phimap_pre")) 
        pdc->writeGlobalVector1D<delphi_real>("phimap_pre", 0, phimap_pre_v.size(), phimap_pre_v);               
    };    
    
    void merge_global() {};
    
    virtual void mpi_run();     
#endif  

};

#endif // SOLVER_FASTSOR_H
