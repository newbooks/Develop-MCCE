/*
 * energy.h
 *
 *  Created on: Feb 23, 2014
 *      Author: Lin Wang, lwang3@g.clemson.edu
 *
 *  This file declares the public and private variables in CDelphiEnergy Class.
 *
 */

#ifndef ENERGY_H_
#define ENERGY_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include "../interface/interface_abstractmodule.h"
#include "../misc/misc_timer.h"
#include "../misc/misc_grid.h"
#include "energy_exceptions.h"

#ifdef PARALLEL_MPI
#include <mpi.h>
#endif

#define NUMPRECISION 2

using namespace std;

class CDelphiEnergy: virtual public IAbstractModule {

private:

    shared_ptr<CTimer> pTimer;

    const int& iBndyType;
    const string& strEnergyFile;
    const string& strScrgFile;
    const bool& bEngOut;
    const bool& bSurfEngOut;
    const bool& bSurfCrgOut;
    const int& iSurfCrgFormatOut;
    const bool& bGridEng;
    const bool& bSolvEng;
    const bool& bAnalySurfEng;
    const bool& bAnalyEng;
    const bool& bIonsEng;
    const bool& bCoulombEng;
    const bool& bPotentiallnSite;
    const bool& bReactFieldlnFRC;
    const bool& bBuffz;
    const delphi_integer& iCrgGridNum;
    const delphi_real& fEpsOut;
    const bool& bNonlinearEng;
    const delphi_integer& iMediaNum;
    const int& iGaussian; //for Gaussian
    const int& inhomo;

    const delphi_real& prgrMediaEps;
    const delphi_integer& prgiAtomMediaNum;
    const delphi_integer& iTotalBdyGridNum;
    const delphi_integer& iCrgBdyGrid;
    vector<SGrid<delphi_integer> >& prgigBndyGrid;
    vector<SGridValue<delphi_real> >& prggvAtomicCrg;
    const delphi_integer& iGrid;
    const delphi_integer& iCrgedGridB;
    const delphi_real& fEpsin;
    //const bool& DEVELOPER;    // Developer is not from datacontainer now. It changed to the flag.
    vector<delphi_real>& prgfPhimap;
    const delphi_real& fScale;
    const SGrid<delphi_real>& fgPotentialDrop;
    vector<SGrid<delphi_real> >& prgfgCrgPoseA;
    vector<SGrid<delphi_integer> >& prgigGridCrgPose;
    vector<delphi_real>& prgfGridCrg;
    delphi_real fIonStrength;
    const SGrid<delphi_real>& fgBoxCenter;
    vector<char>& prgbDielecMap;
    vector<delphi_real>& prgfAtomEps;
    const SExtrema<delphi_integer>& ieBuffz;

    const delphi_integer& iObjectNum;
    const delphi_integer& iAtomNum;
    const delphi_integer& iDielecBndyOdd;
    const delphi_real& fEPKT;
    vector<CAtomPdb>& prgapAtomPdb;
    vector<SDoubleGridValue>& prgdgvCrgBndyGrid;
    vector<SGrid<delphi_real> >& prgfgSurfCrgA;
    vector<delphi_integer>& prgiCrgAt;
    vector<delphi_integer>& atsurf;

    const delphi_real& fTaylorCoeff1;
    const delphi_real& fTaylorCoeff2;
    const delphi_real& fTaylorCoeff3;
    const delphi_real& fTaylorCoeff4;
    const delphi_real& fTaylorCoeff5;

    const char* infoString;
    const char* timeString;
    const char* enerString;
    size_t MAXWIDTH;
    size_t NUMWIDTH;
    vector<SGridValue<delphi_real> > sout;
    SGrid<delphi_integer> lim_min, lim_max;

    vector<delphi_real>& schrg;
    delphi_real& ergg;
    delphi_real& ergc;
    delphi_real& ergs;
    delphi_real& ergr;
    delphi_real& ergions;
    delphi_real& ergsgaussian;

    void energy_clb(delphi_real& ergc);
    void energy_clbmedia(delphi_real& ergc);
    void energy_clbnonl(delphi_real& ergc, delphi_real& ergest, int& igridout);
    void energy_clbtotal(delphi_real& ergest, delphi_real& ergc);
    void energy_react(delphi_real& ergs, delphi_real& ergas, int& iisitpot);
    void energy_nonl(delphi_real& ergnl, int& igridout);

#ifdef PARALLEL_MPI    
    /*
     * mpi version of local functions
     */
    void mpi_energy_clb(delphi_real& ergc);
    void mpi_energy_clbmedia(delphi_real& ergc);
    void mpi_energy_clbnonl(delphi_real& ergc, delphi_real& ergest, int& igridout);
    void mpi_energy_clbtotal(delphi_real& ergest, delphi_real& ergc);
    void mpi_energy_react(delphi_real& ergs, delphi_real& ergas, int& iisitpot);
#endif
    
public:
    CDelphiEnergy(shared_ptr<IDataContainer> pdc, shared_ptr<CTimer> pt) :

    IAbstractModule(pdc), 
    pTimer(pt),

    iBndyType(pdc->getKey_constRef<int>("ibctyp")), 
    strEnergyFile(pdc->getKey_constRef < string > ("nrgnam")), 
    strScrgFile(pdc->getKey_constRef < string > ("scrgnam")), 
    bEngOut(pdc->getKey_constRef<bool>("inrgwrt")), 
    bSurfEngOut(pdc->getKey_constRef<bool>("isen")), 
    bSurfCrgOut(pdc->getKey_constRef<bool>("isch")), 
    iSurfCrgFormatOut(pdc->getKey_constRef<int>("scrgfrm")), 
    bGridEng(pdc->getKey_constRef<bool>("logg")), 
    bSolvEng(pdc->getKey_constRef<bool>("logs")), 
    bAnalySurfEng(pdc->getKey_constRef<bool>("logas")), 
    bAnalyEng(pdc->getKey_constRef<bool>("loga")), 
    bIonsEng(pdc->getKey_constRef<bool>("logions")), 
    bCoulombEng(pdc->getKey_constRef<bool>("logc")), 
    bPotentiallnSite(pdc->getKey_constRef<bool>("isitpot")), 
    bReactFieldlnFRC(pdc->getKey_constRef<bool>("irea")), 
    bBuffz(pdc->getKey_constRef<bool>("ibufz")), 
    iCrgGridNum(pdc->getKey_constRef < delphi_integer > ("nqass")), 
    fEpsOut(pdc->getKey_constRef < delphi_real > ("epsout")), 
    bNonlinearEng(pdc->getKey_constRef<bool>("lognl")), 
    iMediaNum(pdc->getKey_constRef < delphi_integer > ("nmedia")), 
    iGaussian(pdc->getKey_constRef<int>("gaussian")), 
    inhomo(pdc->getKey_constRef<int>("inhomo")), 
    prgrMediaEps(pdc->getKey_constRef < delphi_real > ("medeps")), 
    prgiAtomMediaNum(pdc->getKey_constRef < delphi_integer > ("iatmmed")), 
    iTotalBdyGridNum(pdc->getKey_constRef < delphi_integer > ("ibnum")), 
    iCrgBdyGrid(pdc->getKey_constRef < delphi_integer > ("ibc")), 
    prgigBndyGrid(pdc->getKey_Ref < vector<SGrid<delphi_integer> > > ("ibgrd")),
    iGrid(pdc->getKey_constRef < delphi_integer > ("igrid")), 
    iCrgedGridB(pdc->getKey_constRef < delphi_integer > ("icount1b")), 
    fEpsin(pdc->getKey_constRef < delphi_real > ("epsin")),
    //DEVELOPER(pdc->getKey_constRef<bool>("ideveloper")),
    prgfPhimap(pdc->getKey_Ref < vector<delphi_real> > ("phimap")), // phimap read function to convert to 3d array. //
    fScale(pdc->getKey_constRef < delphi_real > ("scale")), 
    fgPotentialDrop(pdc->getKey_constRef < SGrid<delphi_real> > ("vdrop")), 
    prgfgCrgPoseA(pdc->getKey_Ref < vector<SGrid<delphi_real> > > ("chgpos")),
    prgfGridCrg(pdc->getKey_Ref < vector<delphi_real> > ("gchrg")), 
    prgigGridCrgPose(pdc->getKey_Ref < vector<SGrid<delphi_integer> > > ("gchrgp")),
    prggvAtomicCrg(pdc->getKey_Ref < vector<SGridValue<delphi_real> > > ("atmcrg")),
    fIonStrength(pdc->getKey_Val < delphi_real > ("rionst")), 
    fgBoxCenter(pdc->getKey_constRef < SGrid<delphi_real> > ("oldmid")), 
    prgbDielecMap(pdc->getKey_Ref < vector<char> > ("idebmap")), 
    prgfAtomEps(pdc->getKey_Ref < vector<delphi_real> > ("atmeps")),
    ieBuffz(pdc->getKey_constRef < SExtrema<delphi_integer> > ("buffz")),
    iObjectNum(pdc->getKey_constRef<delphi_integer>("nobject")),
    iAtomNum(pdc->getKey_constRef<delphi_integer>("natom")),
    iDielecBndyOdd(pdc->getKey_constRef<delphi_integer>("icount2b")),
    fEPKT(pdc->getKey_constRef < delphi_real > ("epkt")), 
    prgdgvCrgBndyGrid(pdc->getKey_Ref < vector<SDoubleGridValue> > ("cgbp")), 
    prgapAtomPdb(pdc->getKey_Ref < vector<CAtomPdb> > ("delphipdb")), 
    prgfgSurfCrgA(pdc->getKey_Ref < vector<SGrid<delphi_real> > > ("scspos")),
    prgiCrgAt(pdc->getKey_Ref < vector<delphi_integer> > ("crgatn")),
    atsurf(pdc->getKey_Ref < vector<delphi_integer> > ("atsurf")),
    fTaylorCoeff1(pdc->getKey_constRef < delphi_real > ("chi1")), 
    fTaylorCoeff2(pdc->getKey_constRef < delphi_real > ("chi2")), 
    fTaylorCoeff3(pdc->getKey_constRef < delphi_real > ("chi3")), 
    fTaylorCoeff4(pdc->getKey_constRef < delphi_real > ("chi4")), 
    fTaylorCoeff5(pdc->getKey_constRef < delphi_real > ("chi5")), 
    
    schrg(pdc->getKey_Ref < vector<delphi_real> > ("schrg")), 
    ergg(pdc->getKey_Ref < delphi_real > ("ergg")), 
    ergc(pdc->getKey_Ref < delphi_real > ("ergc")), 
    ergs(pdc->getKey_Ref < delphi_real > ("ergs")), 
    ergr(pdc->getKey_Ref < delphi_real > ("ergr")), 
    ergions(pdc->getKey_Ref < delphi_real > ("ergions")), 
    ergsgaussian(pdc->getKey_Ref < delphi_real > ("ergsgaussian"))  
    {
        #ifdef DEBUG_OBJECT
        cout << endl;
        cout << "****************************************************************\n";
        cout << "*                CDelphiEnergy is constructed                  *\n";
        cout << "****************************************************************\n";
        #endif

#ifdef PARALLEL_MPI     
        mpi_num_procs   = pdc->num_procs();
        mpi_num_workers = mpi_num_procs - 1;
        mpi_rank        = pdc->myid();  
#endif  
    };

    ~CDelphiEnergy() 
    {
        #ifdef DEBUG_OBJECT
        cout << endl;
        cout << "****************************************************************\n";
        cout << "*                 CDelphiEnergy is destroyed                   *\n";
        cout << "****************************************************************\n";
        #endif

#ifdef PARALLEL_MPI
        /*
         * clear allocated memory in local data containers
         */
        if (prgigBndyGrid.size()     > 0) vector< SGrid<delphi_integer> >().swap(prgigBndyGrid);    //prgigBndyGrid.clear();
        if (prggvAtomicCrg.size()    > 0) vector< SGridValue<delphi_real> >().swap(prggvAtomicCrg); //prggvAtomicCrg.clear();
        if (prgfPhimap.size()        > 0) vector<delphi_real>().swap(prgfPhimap);                   //prgfPhimap.clear();
        if (prgfgCrgPoseA.size()     > 0) vector< SGrid<delphi_real> >().swap(prgfgCrgPoseA);       //prgfgCrgPoseA.clear();
        if (prgigGridCrgPose.size()  > 0) vector< SGrid<delphi_integer> >().swap(prgigGridCrgPose); //prgigGridCrgPose.clear();
        if (prgfGridCrg.size()       > 0) vector<delphi_real>().swap(prgfGridCrg);                  //prgfGridCrg.clear();
        if (prgbDielecMap.size()     > 0) vector<char>().swap(prgbDielecMap);                       //prgbDielecMap.clear();
        if (prgfAtomEps.size()       > 0) vector<delphi_real>().swap(prgfAtomEps);                  //prgfAtomEps.clear();
        if (prgapAtomPdb.size()      > 0) vector<CAtomPdb>().swap(prgapAtomPdb);                    //prgapAtomPdb.clear();
        if (prgdgvCrgBndyGrid.size() > 0) vector<SDoubleGridValue>().swap(prgdgvCrgBndyGrid);       //prgdgvCrgBndyGrid.clear();
        if (prgfgSurfCrgA.size()     > 0) vector< SGrid<delphi_real> >().swap(prgfgSurfCrgA);       //prgfgSurfCrgA.clear();
        if (prgiCrgAt.size()         > 0) vector<delphi_integer>().swap(prgiCrgAt);                 //prgiCrgAt.clear();
        if (atsurf.size()            > 0) vector<delphi_integer>().swap(atsurf);                    //atsurf.clear();
        if (schrg.size()             > 0) vector<delphi_real>().swap(schrg);                        //schrg.clear();
#endif
    };

    virtual void validateInput() { };
    virtual void run();

#ifdef PARALLEL_MPI  
    /*
     * mpi variables for general information
     */
    int mpi_num_procs;      // total # of processes
    int mpi_num_workers;    // total # of slave processes = total # of processes - 1
    int mpi_rank;           // process rank or id
    MPI_Status mpi_stat;    // mpi status
    
    void reallocate_global() 
    {
        ddm::Shared<size_t> globalSize;
        
        pdc->deallocateGlobal1D<delphi_real>("schrg");
        if (pdc->myid() == 0) globalSize.set(schrg.size());
        ddm::barrier();
        pdc->allocateGlobal1D<delphi_real>("schrg", globalSize.get());
    };
    
    void write_to_global()
    {
        pdc->writeGlobalVector1D<delphi_real>("schrg",0, schrg.size(), schrg);
        pdc->writeGlobalVar<delphi_real>("ergg",         ergg);
        pdc->writeGlobalVar<delphi_real>("ergc",         ergc);
        pdc->writeGlobalVar<delphi_real>("ergs",         ergs);
        pdc->writeGlobalVar<delphi_real>("ergr",         ergr);
        pdc->writeGlobalVar<delphi_real>("ergions",      ergions);
        pdc->writeGlobalVar<delphi_real>("ergsgaussian", ergsgaussian);
    };

    void merge_global() {};
    
    void mpi_run(); 
#endif     

};

#endif /* ENERGY_H_ */
