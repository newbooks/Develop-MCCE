#include "../dplt/dplt.h"
#include "delphi_datamarshal.h"
#include "delphi_data.h"


void CDelphiData::setMap_dplt()
{
    //--------------------- uniform parameters ---------------------//
    //myData["biomodel"] = pddm->strBioModel;
    createGlobalChar("biomodel", pddm->strBioModel.size());
    //myData["solver"] = pddm->strNumSolver;
    createGlobalChar("solver", pddm->strNumSolver.size());
    //--------------------- set by Statements ----------------------//
    //myData["iautocon"] = pddm->bAutoConverge;
    createGlobalObj<ddm::Shared<bool>>("iautocon");
    //myData["ibctyp"] = pddm->iBndyType;
    createGlobalObj<ddm::Shared<int>>("ibctyp");
    //myData["perfil"] = pddm->fPercentageFill;
    createGlobalObj<ddm::Shared<delphi_real>>("perfil");
    //myData["icheb"] = pddm->bFixedRelaxParam;
    createGlobalObj<ddm::Shared<bool>>("icheb");
    //myData["isrf"] = pddm->bOutGraspSurf;
    createGlobalObj<ddm::Shared<bool>>("isrf");
    //myData["icon2"] = pddm->iConvergeFract;
    createGlobalObj<ddm::Shared<int>>("icon2");
    //myData["icon1"] = pddm->iIterateInterval;
    createGlobalObj<ddm::Shared<int>>("icon1");
    //myData["iexun"] = pddm->bExitUniformDielect;
    createGlobalObj<ddm::Shared<bool>>("iexun");
    //myData["repsout"] = pddm->fExDielec; 
    createGlobalObj<ddm::Shared<delphi_real>>("repsout");
    //myData["isph"] = pddm->bCrgInterplateType;
    createGlobalObj<ddm::Shared<bool>>("isph");
    //myData["gten"] = pddm->fGridConverge;
    createGlobalObj<ddm::Shared<delphi_real>>("gten");
    //myData["igrid"] = pddm->iGrid;
    createGlobalObj<ddm::Shared<delphi_integer>>("igrid");
    //myData["repsin"] = pddm->fInDielec;
    createGlobalObj<ddm::Shared<delphi_real>>("repsin");
    //myData["conc"] = pddm->vctfSalt; // std::vector
    createGlobalObj<dplt::global1D<delphi_real>>("conc");
    //myData["exrad"] = pddm->fIonRadius;
    createGlobalObj<ddm::Shared<delphi_real>>("exrad");
    //myData["nlit"] = pddm->iLinIterateNum;
    createGlobalObj<ddm::Shared<int>>("nlit");
    //myData["igraph"] = pddm->bLogGraph;
    createGlobalObj<ddm::Shared<bool>>("igraph");
    //myData["ipoten"] = pddm->bLogPotential;
    createGlobalObj<ddm::Shared<bool>>("ipoten");
    //myData["res2"] = pddm->fMaxc;
    createGlobalObj<ddm::Shared<delphi_real>>("res2");
    //myData["nnit"] = pddm->iNonIterateNum;
    createGlobalObj<ddm::Shared<int>>("nnit");
    //myData["iper"] = pddm->vctbPeriodicBndy; // std::vector
    createGlobalObj<dplt::global1D<char>>("iper");
    //myData["iconc"] = pddm->bOutCrgDensity;
    createGlobalObj<ddm::Shared<bool>>("iconc");
    //myData["radprb"] = pddm->vctfProbeRadius; // std::vector  
    createGlobalObj<dplt::global1D<delphi_real>>("radprb");
    //myData["uspec"] = pddm->fSpectralRadius;
    createGlobalObj<ddm::Shared<delphi_real>>("uspec");
    //myData["relpar"] = pddm->fRelaxParam;
    createGlobalObj<ddm::Shared<delphi_real>>("relpar");
    //myData["res1"] = pddm->fRmsc;
    createGlobalObj<ddm::Shared<delphi_real>>("res1");
    //myData["scale"] = pddm->fScale;
    createGlobalObj<ddm::Shared<delphi_real>>("scale");
    //myData["isolv"] = pddm->bSolvePB;
    createGlobalObj<ddm::Shared<bool>>("isolv");
    //myData["ival"] = pddm->vctiValence1; // std::vector
    createGlobalObj<dplt::global1D<int>>("ival");
    //myData["ival2"] = pddm->vctiValence2; // std::vector
    createGlobalObj<dplt::global1D<int>>("ival2");
    //myData["atompotdist"] = pddm->fPotentialUpperBond;
    createGlobalObj<ddm::Shared<delphi_real>>("atompotdist");
    //myData["temperature"] = pddm->fTemper;
    createGlobalObj<ddm::Shared<delphi_real>>("temperature");
    //myData["vdrop"] = pddm->gfPotentialDrop; // SGrid<delphi_real>
    createGlobalObj<ddm::Shared<SGrid<delphi_real>>>("vdrop");
    //myData["iuspec"] = pddm->bSpectralRadius;
    createGlobalObj<ddm::Shared<bool>>("iuspec");
    //myData["imanual"] = pddm->bManualRelaxParam;
    createGlobalObj<ddm::Shared<bool>>("imanual");
    //myData["phiintype"] = pddm->iPhiInType; //for focusing
    createGlobalObj<ddm::Shared<int>>("phiintype");

    //---------------------Gaussian & MEMPOT ----------------------//
    //Lin Li: for Gaussian & MEMPOT options
    //myData["cutoff"] = pddm->fCutoff;
    createGlobalObj<ddm::Shared<delphi_real>>("cutoff");
    //myData["sigma"] = pddm->fSigma;
    createGlobalObj<ddm::Shared<delphi_real>>("sigma");
    //myData["inhomo"] = pddm->iInhomo;
    createGlobalObj<ddm::Shared<int>>("inhomo");
    //myData["srfcut"] = pddm->fSrfcut;
    createGlobalObj<ddm::Shared<delphi_real>>("srfcut");
    //myData["gaussian"] = pddm->iGaussian;
    createGlobalObj<ddm::Shared<int>>("gaussian");
    //myData["gepsmp"] = pddm->vctgfGepsMap;
    createGlobalObj<dplt::global1D<SGrid<delphi_real>>>("gepsmp"); //Vector
    //myData["gepsmp2"] = pddm->vctgfGepsMap2;
    createGlobalObj<dplt::global1D<SGrid<delphi_real>>>("gepsmp2");//Vector
    //myData["ergsgaussian"] = pddm->fErgsgaussian;
    createGlobalObj<ddm::Shared<delphi_real>>("ergsgaussian");
    //myData["radipz"] = pddm->fRadipz;
    createGlobalObj<ddm::Shared<delphi_real>>("radipz");

    //-------------------Gaussian Salt----------------------------//
    //myData["gdensity"] = pddm->vctfgaussianSoluteDensity;
    createGlobalObj<dplt::global1D<delphi_real>>("gdensity"); //Vector
    //myData["gdtype"] = pddm->iGaussianDensity;
    createGlobalObj<ddm::Shared<delphi_integer>>("gdtype");
    //myData["getype"] = pddm->iGaussianEnergy;
    createGlobalObj<ddm::Shared<delphi_integer>>("getype");

    //-------------------------- io file names ------------------------//
    //  myData["prmnam"]    = pddm->strParamFile[0]; // (not to be mapped)
    //  myData["siznam"] = pddm->strSizeFile;
    createGlobalChar("siznam", pddm->strSizeFile.size());
    //  myData["crgnam"] = pddm->strCrgFile;
    createGlobalChar("crgnam", pddm->strCrgFile.size());
    //  myData["pdbnam"] = pddm->strPdbFile;
    createGlobalChar("pdbnam", pddm->strPdbFile.size());
    //  myData["phinam"] = pddm->strPhiFile;
    createGlobalChar("phinam", pddm->strPhiFile.size());
    //  myData["frcinam"] = pddm->strFrciFile;
    createGlobalChar("frcinam", pddm->strFrciFile.size());
    //  myData["frcnam"] = pddm->strFrcFile;
    createGlobalChar("frcnam", pddm->strFrcFile.size());
    //  myData["epsnam"] = pddm->strEpsFile;
    createGlobalChar("epsnam", pddm->strEpsFile.size());
    //  myData["phiinam"] = pddm->strPhiiFile;
    createGlobalChar("phiinam", pddm->strPhiiFile.size());
    //  myData["mpdbnam"] = pddm->strModifiedPdbFile;
    createGlobalChar("mpdbnam", pddm->strModifiedPdbFile.size());
    //  myData["updbnam"] = pddm->strUnformatPdbFile;
    createGlobalChar("updbnam", pddm->strUnformatPdbFile.size());
    //  myData["ufrcnam"] = pddm->strUnformatFrcFile;
    createGlobalChar("ufrcnam", pddm->strUnformatFrcFile.size());
    //  myData["srfnam"] = pddm->strGraspFile;
    createGlobalChar("srfnam", pddm->strGraspFile.size());
    //  myData["nrgnam"] = pddm->strEnergyFile;
    createGlobalChar("nrgnam", pddm->strEnergyFile.size());
    //  myData["scrgnam"] = pddm->strScrgFile;
    createGlobalChar("scrgnam", pddm->strScrgFile.size());
    //  myData["zphinam"] = pddm->strZetaPhiFile;   //ARGO for ZPHI output file :: ZETA
    createGlobalChar("zphinam", pddm->strZetaPhiFile.size());
    //

    //  myData["offset"] = pddm->gfOffCenter; // SGrid<delphi_real>
    createGlobalObj<ddm::Shared<SGrid<delphi_real>>>("offset");
    //

    // * set by ACENTER or ACENT function

    //  myData["acent"] = pddm->gfAcent; // SGrid<delphi_real>
    createGlobalObj<ddm::Shared<SGrid<delphi_real>>>("acent");
    //  myData["iacent"] = pddm->bIsAcent;
    createGlobalObj<ddm::Shared<bool>>("iacent");
    //

    //  * set by READ or IN function

    //  myData["pdbfrm"] = pddm->iPdbFormatIn;
    createGlobalObj<ddm::Shared<int>>("pdbfrm");
    //  myData["ipdbrd"] = pddm->bPdbUnformatIn;
    createGlobalObj<ddm::Shared<bool>>("ipdbrd");
    //

    //   set by WRITE or OUT function

    //  myData["phiwrt"] = pddm->bPhimapOut;
    createGlobalObj<ddm::Shared<bool>>("phiwrt");
    //  myData["phifrm"] = pddm->iPhiFormatOut;
    createGlobalObj<ddm::Shared<int>>("phifrm");
    //  myData["ibios"] = pddm->bBiosystemOut;
    createGlobalObj<ddm::Shared<bool>>("ibios");
    //  myData["ibem"] = pddm->bBemSrfOut;
    createGlobalObj<ddm::Shared<bool>>("ibem");
    //  myData["isite"] = pddm->bSiteOut;
    createGlobalObj<ddm::Shared<bool>>("isite");
    //  myData["frcfrm"] = pddm->iFrcFormatOut;
    createGlobalObj<ddm::Shared<int>>("frcfrm");
    //  myData["epswrt"] = pddm->bEpsOut;
    createGlobalObj<ddm::Shared<bool>>("epswrt");
    //  myData["iatout"] = pddm->bModPdbOut;
    createGlobalObj<ddm::Shared<bool>>("iatout");
    //  myData["mpdbfrm"] = pddm->iModPdbFormatOut;
    createGlobalObj<ddm::Shared<int>>("mpdbfrm");
    //  myData["ipdbwrt"] = pddm->bUnformatPdbOut;
    createGlobalObj<ddm::Shared<bool>>("ipdbwrt");
    //  myData["ifrcwrt"] = pddm->bUnformatFrcOut;
    createGlobalObj<ddm::Shared<bool>>("ifrcwrt");
    //  myData["inrgwrt"] = pddm->bEngOut;
    createGlobalObj<ddm::Shared<bool>>("inrgwrt");
    //  myData["iwgcrg"] = pddm->bGridCrgOut;
    createGlobalObj<ddm::Shared<bool>>("iwgcrg");
    //  myData["iacs"] = pddm->bHsurf2DatOut;
    createGlobalObj<ddm::Shared<bool>>("iacs");
    //  myData["idbwrt"] = pddm->bDbOut;
    createGlobalObj<ddm::Shared<bool>>("idbwrt");
    //  myData["isen"] = pddm->bSurfEngOut;
    createGlobalObj<ddm::Shared<bool>>("isen");
    //  myData["isch"] = pddm->bSurfCrgOut;
    createGlobalObj<ddm::Shared<bool>>("isch");
    //  myData["scrgfrm"] = pddm->iSurfCrgFormatOut;
    createGlobalObj<ddm::Shared<int>>("scrgfrm");
    //  myData["zphi_out"] = pddm->bZetaPhiOut; //zeta
    createGlobalObj<ddm::Shared<bool>>("zphi_out");
    //

    // * set by ENERGY function

    //  myData["logg"] = pddm->bGridEng;
    createGlobalObj<ddm::Shared<bool>>("logg");
    //  myData["logs"] = pddm->bSolvEng;
    createGlobalObj<ddm::Shared<bool>>("logs");
    //  myData["logas"] = pddm->bAnalySurfEng;
    createGlobalObj<ddm::Shared<bool>>("logas");
    //  myData["loga"] = pddm->bAnalyEng;
    createGlobalObj<ddm::Shared<bool>>("loga");
    //  myData["logions"] = pddm->bIonsEng;
    createGlobalObj<ddm::Shared<bool>>("logions");
    //  myData["logc"] = pddm->bCoulombEng;
    createGlobalObj<ddm::Shared<bool>>("logc");
    //

    //  * set by SITE function: all MUST be initialized to to false

    //  myData["isita"] = pddm->bAtomInSite;
    createGlobalObj<ddm::Shared<bool>>("isita");
    //  myData["isitq"] = pddm->bCrgInSite;
    createGlobalObj<ddm::Shared<bool>>("isitq");
    //  myData["isitp"] = pddm->bGridPotentialInSite;
    createGlobalObj<ddm::Shared<bool>>("isitp");
    //  myData["isitap"] = pddm->bAtomPotentialInSite;
    createGlobalObj<ddm::Shared<bool>>("isitap");
    //  myData["isitdeb"] = pddm->bDebyeFractionInSite;
    createGlobalObj<ddm::Shared<bool>>("isitdeb");
    //  myData["isitf"] = pddm->bFieldInSite;
    createGlobalObj<ddm::Shared<bool>>("isitf");
    //  myData["isitr"] = pddm->bReactPotentialInSite;
    createGlobalObj<ddm::Shared<bool>>("isitr");
    //  myData["isitc"] = pddm->bCoulombPotentialInSite;
    createGlobalObj<ddm::Shared<bool>>("isitc");
    //  myData["isitx"] = pddm->bAtomCoordInSite;
    createGlobalObj<ddm::Shared<bool>>("isitx");
    //  myData["isiti"] = pddm->bSaltInSite;
    createGlobalObj<ddm::Shared<bool>>("isiti");
    //  myData["isitt"] = pddm->bTotalPotentialInSite;
    createGlobalObj<ddm::Shared<bool>>("isitt");
    //  myData["isitrf"] = pddm->bReactForceInSite;
    createGlobalObj<ddm::Shared<bool>>("isitrf");
    //  myData["isitcf"] = pddm->bCoulombForceInSite;
    createGlobalObj<ddm::Shared<bool>>("isitcf");
    //  myData["isitmd"] = pddm->bMDInSite;
    createGlobalObj<ddm::Shared<bool>>("isitmd");
    //  myData["isitsf"] = pddm->bSurfCrgInSite;
    createGlobalObj<ddm::Shared<bool>>("isitsf");
    //  myData["isittf"] = pddm->bTotalForceInSite;
    createGlobalObj<ddm::Shared<bool>>("isittf");
    //  myData["isitpot"] = pddm->bPotentialInSite;
    createGlobalObj<ddm::Shared<bool>>("isitpot");
    //  myData["irea"] = pddm->bReactFieldInFRC;
    createGlobalObj<ddm::Shared<bool>>("irea");
    //  myData["iself"] = pddm->bPDB2FRCInSite;
    createGlobalObj<ddm::Shared<bool>>("iself");
    //

    //  * set by BUFFZ function

    //myData["buffz"] = pddm->eiBuffz;
    createGlobalObj<ddm::Shared<SExtrema<delphi_integer>>>("buffz");
    //  myData["ibufz"] = pddm->bIsBuffz;
    createGlobalObj<ddm::Shared<bool>>("ibufz");
    //

    //  * set by SURFACE function

    //  myData["isurftype"] = pddm->iTypeSurf;
    createGlobalObj<ddm::Shared<int>>("isurftype");
    //  //----------------------- set by DelPhi ------------------------//
    //  myData["deblen"] = pddm->fDebyeLength;
    createGlobalObj<ddm::Shared<delphi_real>>("deblen");
    //  myData["epsout"] = pddm->fEpsOut;
    createGlobalObj<ddm::Shared<delphi_real>>("epsout");
    //  myData["cran"] = pddm->gfCoordinateRange;  //SGrid<delphi_real>
    createGlobalObj<ddm::Shared<SGrid<delphi_real>>>("cran");
    //  myData["pmid"] = pddm->gfGeometricCenter;
    createGlobalObj<ddm::Shared<SGrid<delphi_real>>>("pmid");
    //  myData["oldmid"] = pddm->gfBoxCenter;
    createGlobalObj<ddm::Shared<SGrid<delphi_real>>>("oldmid");
    //  myData["rionst"] = pddm->fIonStrength;
    createGlobalObj<ddm::Shared<delphi_real>>("rionst");
    //  myData["chi1"] = pddm->fTaylorCoeff1;
    createGlobalObj<ddm::Shared<delphi_real>>("chi1");
    //  myData["chi2"] = pddm->fTaylorCoeff2;
    createGlobalObj<ddm::Shared<delphi_real>>("chi2");
    //  myData["chi3"] = pddm->fTaylorCoeff3;
    createGlobalObj<ddm::Shared<delphi_real>>("chi3");
    //  myData["chi4"] = pddm->fTaylorCoeff4;
    createGlobalObj<ddm::Shared<delphi_real>>("chi4");
    //  myData["chi5"] = pddm->fTaylorCoeff5;
    createGlobalObj<ddm::Shared<delphi_real>>("chi5");
    //  myData["lognl"] = pddm->bNonlinearEng;
    createGlobalObj<ddm::Shared<bool>>("lognl");
    //  myData["epkt"] = pddm->fEPKT;
    createGlobalObj<ddm::Shared<delphi_real>>("epkt");
    //  myData["epsin"] = pddm->fEpsIn;
    createGlobalObj<ddm::Shared<delphi_real>>("epsin");
    //  myData["ifrcrd"] = pddm->bFrcUnformatIn;
    createGlobalObj<ddm::Shared<bool>>("ifrcrd");
    //  myData["idirectalg"] = pddm->iDirectEpsMap;
    createGlobalObj<ddm::Shared<int>>("idirectalg");
    //  myData["numbmol"] = pddm->iMoleculeNum;
    createGlobalObj<ddm::Shared<delphi_integer>>("numbmol");
    //  myData["rdmx"] = pddm->fMaxRadius;
    createGlobalObj<ddm::Shared<delphi_real>>("rdmx");
    //  myData["uniformdiel"] = pddm->bUniformDielec;
    createGlobalObj<ddm::Shared<bool>>("uniformdiel");
    //  myData["limobject"] = pddm->vctefExtrema;     // std::vector< SExtrema<delphi_real> >
    createGlobalObj<dplt::global1D<SExtrema<delphi_real>>>("limobject");
    //  myData["xn1"] = pddm->vctgfAtomCoordA;  // std::vector< SGrid<delphi_real> >
    createGlobalObj<dplt::global1D<SGrid<delphi_real>>>("xn1");
    //  myData["xn2"] = pddm->vctgfAtomCoordG;  // std::vector< SGrid<delphi_real> >
    createGlobalObj<dplt::global1D<SGrid<delphi_real>>>("xn2");
    //----------------------- set by IO class ------------------------//
    //  myData["resnummax"] = pddm->iResidueNum;
    createGlobalObj<ddm::Shared<delphi_integer>>("resnummax");
    //  myData["nmedia"] = pddm->iMediaNum;
    createGlobalObj<ddm::Shared<delphi_integer>>("nmedia");
    //  myData["medeps"] = pddm->vctfMediaEps;     // std::vector<delphi_real>
    createGlobalObj<dplt::global1D<delphi_real>>("medeps");
    //  myData["nobject"] = pddm->iObjectNum;
    createGlobalObj<ddm::Shared<delphi_integer>>("nobject");
    //  myData["dataobject"] = pddm->vctstrObject;     // std::vector<string> no need to broadcast, used in space only
    createGlobal1DString("dataobject", pddm->vctstrObject);
    //  myData["natom"] = pddm->iAtomNum;
    createGlobalObj<ddm::Shared<delphi_integer>>("natom");
    //  myData["delphipdb"] = pddm->vctapAtomPdb;     // std::vector<CAtomPdb>
    createGlobalObj<dplt::global1D<CAtomPdb>>("delphipdb");
    //  myData["iatmmed"] = pddm->vctiAtomMediaNum; // std::vector<delphi_integer>
    createGlobalObj<dplt::global1D<delphi_integer>>("iatmmed");
    //  myData["ionlymol"] = pddm->bOnlyMolecule;
    createGlobalObj<ddm::Shared<bool>>("ionlymol");

    //------------------- set by Surface class ---------------------//
    //  myData["nqass"] = pddm->iCrgGridNum;
    createGlobalObj<ddm::Shared<delphi_integer>>("nqass");
    //  myData["qnet"] = pddm->fNetCrg;
    createGlobalObj<ddm::Shared<delphi_real>>("qnet");
    //  myData["qmin"] = pddm->fMinusCrg;
    createGlobalObj<ddm::Shared<delphi_real>>("qmin");
    //  myData["qplus"] = pddm->fPlusCrg;
    createGlobalObj<ddm::Shared<delphi_real>>("qplus");
    //  myData["cqplus"] = pddm->gfPlusCrgCenter;
    createGlobalObj<ddm::Shared<SGrid<delphi_real>>>("cqplus");
    //  myData["cqmin"] = pddm->gfMinusCrgCenter;
    createGlobalObj<ddm::Shared<SGrid<delphi_real>>>("cqmin");
    //  myData["cmin"] = pddm->gfMinCoordinate;
    createGlobalObj<ddm::Shared<SGrid<delphi_real>>>("cmin");
    //  myData["cmax"] = pddm->gfMaxCoordinate;
    createGlobalObj<ddm::Shared<SGrid<delphi_real>>>("cmax");
    //  myData["ibnum"] = pddm->iBndyGridNum;
    createGlobalObj<ddm::Shared<delphi_integer>>("ibnum");
    //  myData["iepsmp"] = pddm->vctgiEpsMap;      // std::vector< SGrid<delphi_integer> >
    createGlobalObj<dplt::global1D<SGrid<delphi_integer> >>("iepsmp");
    //  myData["idebmap"] = pddm->vctbDielecMap;    // std::vector<bool>  here uses char
    createGlobalObj<dplt::global1D<char>>("idebmap");
    //
    //                                              //ARGO: For zeta surface map
    //  myData["zetaSurfMap"] = pddm->vctZetaSurfMap;    // std::vector<bool> here uses char
    createGlobalObj<dplt::global1D<char>>("zetaSurfMap");
    //  myData["zetaOn"] = pddm->zetaOn;  //int
    createGlobalObj<ddm::Shared<int>>("zetaOn");
    //  myData["zetaDistance"] = pddm->zetaDistance;  //delphi_real
    createGlobalObj<ddm::Shared<delphi_real>>("zetaDistance");

    //  myData["ibgrd"] = pddm->vctgiBndyGrid;    // std::vector< SGrid<delphi_integer> >
    createGlobalObj<dplt::global1D<SGrid<delphi_integer>>>("ibgrd");
    //  myData["nqgrd"] = pddm->iCrg2GridNum;
    createGlobalObj<ddm::Shared<delphi_integer>>("nqgrd");
    //  myData["chrgv2"] = pddm->vctgvfCrg2Grid;   // std::vector< SGridValue<delphi_real> >
    createGlobalObj<dplt::global1D<SGridValue<delphi_real>>>("chrgv2");
    //  myData["nqgrdtonqass"] = pddm->vctiCrg2GridMap;  // std::vector<delphi_integer>
    createGlobalObj<dplt::global1D<delphi_integer>>("nqgrdtonqass");
    //  myData["atmcrg"] = pddm->vctgvfAtomCrg;    // std::vector< SGridValue<delphi_real> >
    createGlobalObj<dplt::global1D<SGridValue<delphi_real>>>("atmcrg");
    //  myData["chgpos"] = pddm->vctgfCrgPoseA;    // std::vector< SGrid<delphi_real> >
    createGlobalObj<dplt::global1D<SGrid<delphi_real>>>("chgpos");
    //  myData["scspos"] = pddm->vctgfSurfCrgA;    // std::vector< SGrid<delphi_real> >
    createGlobalObj<dplt::global1D<SGrid<delphi_real>>>("scspos");
    //  myData["crgatn"] = pddm->vctiCrgAt;        // std::vector<delphi_integer>
    createGlobalObj<dplt::global1D<delphi_integer>>("crgatn");
    //  myData["atsurf"] = pddm->vctiAtSurf;       // std::vector<delphi_integer>
    createGlobalObj<dplt::global1D<delphi_integer>>("atsurf");
    //  myData["atndx"] = pddm->vctiAtNdx;        // std::vector<delphi_integer>
    createGlobalObj<dplt::global1D<delphi_integer>>("atndx");
    //  myData["scsnor"] = pddm->vctgfSurfCrgE;    // std::vector< SGrid<delphi_real> >
    createGlobalObj<dplt::global1D<SGrid<delphi_real>>>("scsnor");
    //  myData["atmeps"] = pddm->vctfAtomEps;      // std::vector<delphi_real>
    createGlobalObj<dplt::global1D<delphi_real>>("atmeps");
    //------------------- set by Solver class ---------------------//
    //  myData["icount2b"] = pddm->iDielecBndySum;
    createGlobalObj<ddm::Shared<delphi_integer>>("icount2b");
    //  myData["icount1b"] = pddm->iCrgedGridSum;
    createGlobalObj<ddm::Shared<delphi_integer>>("icount1b");
    //  myData["gchrg"] = pddm->vctfGridCrg;      // std::vector<delphi_real>
    createGlobalObj<dplt::global1D<delphi_real>>("gchrg");
    //  myData["gchrgp"] = pddm->vctgiGridCrgPose; // std::vector< SGrid<delphi_integer> >
    createGlobalObj<dplt::global1D<SGrid<delphi_integer>>>("gchrgp");
    //  myData["ibc"] = pddm->iCrgBdyGrid;
    createGlobalObj<ddm::Shared<delphi_integer>>("ibc");
    //  myData["cgbp"] = pddm->vctdgvCrgBndyGrid;// std::vector<SDoubleGridValue>
    createGlobalObj<dplt::global1D<SDoubleGridValue>>("cgbp");
    //  myData["phimap"] = pddm->vctfPhiMap;       // std::vector<delphi_real>
    createGlobalObj<dplt::global1D<delphi_real>>("phimap");
    //  myData["phimap_pre"] = pddm->vctfPhiMap_Pre;       // previous phimap for focusing
    createGlobalObj<dplt::global1D<delphi_real>>("phimap_pre");
    //
    //------------------- set by Energy class ---------------------//
    //  myData["schrg"] = pddm->vctfSurfCrgE;     // std::vector<delphi_real>
    createGlobalObj<dplt::global1D<delphi_real>>("schrg");
    //  myData["ergg"] = pddm->fEngGrid;
    createGlobalObj<ddm::Shared<delphi_real>>("ergg");
    //  myData["ergc"] = pddm->fEngCoul;
    createGlobalObj<ddm::Shared<delphi_real>>("ergc");
    //  myData["ergs"] = pddm->fEngCorrect;
    createGlobalObj<ddm::Shared<delphi_real>>("ergs");
    //  myData["ergr"] = pddm->fEngReact;
    createGlobalObj<ddm::Shared<delphi_real>>("ergr");
    //  myData["ergions"] = pddm->fEngIons;
    createGlobalObj<ddm::Shared<delphi_real>>("ergions");
    //  //------------------- New var by Interface  ---------------------//
#ifdef PRIME
    //  myData["vecfrcin"] = pddm->strCommFRCIn;    // update: Dec 19, 2014 by Lin Wang
#endif
    //  myData["isitcomm"] = pddm->bCommFRCIn;      // update: Dec 19, 2014 by Lin Wang
    createGlobalObj<ddm::Shared<bool>>("isitcomm");
    ddm::barrier();
}

void CDelphiData::allocate_dplt()
{
    ddm::Shared<size_t> globalSize;

    //Update variables set by statements
    if (myid() == 0) globalSize.set(pddm->vctfSalt.size());
    ddm::barrier();
    allocateGlobal1D<delphi_real>("conc", globalSize.get());
    ddm::barrier();

    if (myid() == 0) globalSize.set(pddm->vctbPeriodicBndy.size());
    ddm::barrier();
    allocateGlobal1D<char>("iper", globalSize.get());
    ddm::barrier();

    if (myid() == 0) globalSize.set(pddm->vctfProbeRadius.size());
    ddm::barrier();
    allocateGlobal1D<delphi_real>("radprb", globalSize.get());
    ddm::barrier();

    if (myid() == 0) globalSize.set(pddm->vctiValence1.size());
    ddm::barrier();
    allocateGlobal1D<int>("ival", globalSize.get());
    ddm::barrier();

    if (myid() == 0) globalSize.set(pddm->vctiValence2.size());
    ddm::barrier();
    allocateGlobal1D<int>("ival2", globalSize.get());
    ddm::barrier();

    //From Delphi Class
    if (myid() == 0) globalSize.set(pddm->vctefExtrema.size());
    ddm::barrier();
    allocateGlobal1D<SExtrema<delphi_real>>("limobject", globalSize.get());
    ddm::barrier();

    if (myid() == 0) globalSize.set(pddm->vctgfAtomCoordA.size());
    ddm::barrier();
    allocateGlobal1D<SGrid<delphi_real>>("xn1", globalSize.get());
    ddm::barrier();

    if (myid() == 0) globalSize.set(pddm->vctgfAtomCoordG.size());
    ddm::barrier();
    allocateGlobal1D<SGrid<delphi_real>>("xn2", globalSize.get());
    ddm::barrier();
    //allocateGlobal1D<SExtrema<delphi_real>>("limobject", 1);
    //allocateGlobal1D<SGrid<delphi_real>>("xn1", 1);
    //allocateGlobal1D<SGrid<delphi_real>>("xn2", 1);

    // From IO Class
    if (myid() == 0) globalSize.set(pddm->vctfMediaEps.size());
    ddm::barrier();
    allocateGlobal1D<delphi_real>("medeps", globalSize.get());
    ddm::barrier();

    if (myid() == 0) globalSize.set(pddm->vctapAtomPdb.size());
    ddm::barrier();
    allocateGlobal1D<CAtomPdb>("delphipdb", globalSize.get());
    ddm::barrier();

    if (myid() == 0) globalSize.set(pddm->vctiAtomMediaNum.size());
    ddm::barrier();
    allocateGlobal1D<delphi_integer>("iatmmed", globalSize.get());
    ddm::barrier();
}

void CDelphiData::updateParameters_dplt()
{
    //Update variables set by statements
    // iautocon     bAutoConverge
    writeGlobalVar<bool>("iautocon", pddm->bAutoConverge);
    // ibctyp     iBndyType
    writeGlobalVar<int>("ibctyp", pddm->iBndyType);
    // perfil     fPercentageFill
    writeGlobalVar<delphi_real>("perfil", pddm->fPercentageFill);
    // icheb     bFixedRelaxParam
    writeGlobalVar<bool>("icheb", pddm->bFixedRelaxParam);
    // isrf     bOutGraspSurf
    writeGlobalVar<bool>("isrf", pddm->bOutGraspSurf);
    // icon2     iConvergeFract
    writeGlobalVar<int>("icon2", pddm->iConvergeFract);
    // icon1     iIterateInterval
    writeGlobalVar<int>("icon1", pddm->iIterateInterval);
    // iexun     bExitUniformDielect
    writeGlobalVar<bool>("iexun", pddm->bExitUniformDielect);
    // repsout     fExDielec
    writeGlobalVar<delphi_real>("repsout", pddm->fExDielec);
    // isph     bCrgInterplateType
    writeGlobalVar<bool>("isph", pddm->bCrgInterplateType);
    // gten     fGridConverge
    writeGlobalVar<delphi_real>("gten", pddm->fGridConverge);
    // repsin     fInDielec
    writeGlobalVar<delphi_real>("repsin", pddm->fInDielec);
    // igrid     iGrid
    writeGlobalVar<delphi_integer>("igrid", pddm->iGrid);

    // exrad     fIonRadius
    writeGlobalVar<delphi_real>("exrad", pddm->fIonRadius);

    // conc     rgfSaltConcent[2]
    writeGlobalVector1D<delphi_real>("conc", 0, pddm->vctfSalt.size(), pddm->vctfSalt);

    // nlit     iLinIterateNum
    writeGlobalVar<int>("nlit", pddm->iLinIterateNum);
    // igraph     bLogGraph
    writeGlobalVar<bool>("igraph", pddm->bLogGraph);
    // ipoten     bLogPotential
    writeGlobalVar<bool>("ipoten", pddm->bLogPotential);
    // res2     fMaxc
    writeGlobalVar<delphi_real>("res2", pddm->fMaxc);

    // nnit     iNonIterateNum
    writeGlobalVar<int>("nnit", pddm->iNonIterateNum);
    // iper     rgbPeriodicBndy[6]  
    writeGlobalVector1D<char>("iper", 0, pddm->vctbPeriodicBndy.size(), pddm->vctbPeriodicBndy);
    // iconc     bOutCrgDensity
    writeGlobalVar<bool>("iconc", pddm->bOutCrgDensity);
    // radprb     rgfProbeRadius[2]
    writeGlobalVector1D<delphi_real>("radprb", 0, pddm->vctfProbeRadius.size(), pddm->vctfProbeRadius);
    // uspec     fSpectralRadius
    writeGlobalVar<delphi_real>("uspec", pddm->fSpectralRadius);
    // relpar     fRelaxParam
    writeGlobalVar<delphi_real>("relpar", pddm->fRelaxParam);
    // res1     fRmsc
    writeGlobalVar<delphi_real>("res1", pddm->fRmsc);
    // scale     fScale
    writeGlobalVar<delphi_real>("scale", pddm->fScale);
    // isolv     bSolvePB
    writeGlobalVar<bool>("isolv", pddm->bSolvePB);
    // ival     rgiValence1[2]
    writeGlobalVector1D<int>("ival", 0, pddm->vctiValence1.size(), pddm->vctiValence1);
    // ival2     rgiValence2[2]
    writeGlobalVector1D<int>("ival2", 0, pddm->vctiValence2.size(), pddm->vctiValence2);
    // atompotdist     fPotentialUpperBond
    writeGlobalVar<delphi_real>("atompotdist", pddm->fPotentialUpperBond);
    // temperature     fTemper
    writeGlobalVar<delphi_real>("temperature", pddm->fTemper);
    // vdrop     fgPotentialDrop
    writeGlobalVar<SGrid<delphi_real>>("vdrop", pddm->gfPotentialDrop);
    // iuspec     bSpectralRadius
    writeGlobalVar<bool>("iuspec", pddm->bSpectralRadius);
    // imanual     bManualRelaxParam
    writeGlobalVar<bool>("imanual", pddm->bManualRelaxParam);
    // sigma     fSigma
    writeGlobalVar<delphi_real>("sigma", pddm->fSigma);
    // cutoff     fCutoff
    writeGlobalVar<delphi_real>("cutoff", pddm->fCutoff);
    // srfcut     fSrfcut
    writeGlobalVar<delphi_real>("srfcut", pddm->fSrfcut);

    if(debug_data) cout << "*                 Delphi data broadcast statements done        *\n";

    //File names

    // siznam     strSizeFile
    writeGlobalChar("siznam", pddm->strSizeFile);
    // crgnam     strCrgFile
    writeGlobalChar("crgnam", pddm->strCrgFile);
    // pdbnam     strPdbFile
    writeGlobalChar("pdbnam", pddm->strPdbFile);
    // phinam     strPhiFile
    writeGlobalChar("phinam", pddm->strPhiFile);
    // frcinam     strFrciFile
    writeGlobalChar("frcinam", pddm->strFrciFile);
    // frcnam     strFrcFile
    writeGlobalChar("frcnam", pddm->strFrcFile);
    // epsnam     strEpsFile
    writeGlobalChar("epsnam", pddm->strEpsFile);
    // phiinam     strPhiiFile
    writeGlobalChar("phiinam", pddm->strPhiiFile);
    // mpdbnam     strModifiedPdbFile
    writeGlobalChar("mpdbnam", pddm->strModifiedPdbFile);
    // updbnam     strUnformatPdbFile
    writeGlobalChar("updbnam", pddm->strUnformatPdbFile);
    // ufrcnam     strUnformatFrcFile
    writeGlobalChar("ufrcnam", pddm->strUnformatFrcFile);
    // srfnam     strGraspFile
    writeGlobalChar("srfnam", pddm->strGraspFile);
    // nrgnam     strEnergyFile
    writeGlobalChar("nrgnam", pddm->strEnergyFile);
    // scrgnam     strScrgFile
    writeGlobalChar("scrgnam", pddm->strScrgFile);

    if (debug_data)cout << "*                 Delphi data broadcast file names done        *\n";

    //      Set by Functions
    //myData["biomodel"] = pddm->strBioModel;
    writeGlobalChar("biomodel", pddm->strBioModel);
    //myData["solver"] = pddm->strNumSolver;
    writeGlobalChar("solver", pddm->strNumSolver);
    // offset     fgOffCenter
    writeGlobalVar<SGrid<delphi_real>>("offset", pddm->gfOffCenter);
    // acent     fgAcent
    writeGlobalVar<SGrid<delphi_real>>("acent", pddm->gfAcent);
    // iacent     bIsAcent
    writeGlobalVar<bool>("iacent", pddm->bIsAcent);
    // pdbfrm     iPdbFormatIn
    writeGlobalVar<int>("pdbfrm", pddm->iPdbFormatIn);
    // ipdbrd     bPdbUnformatIn
    writeGlobalVar<bool>("ipdbrd", pddm->bPdbUnformatIn);
    // phiwrt     bPhimapOut
    writeGlobalVar<bool>("phiwrt", pddm->bPhimapOut);
    // phifrm     iPhiFormatOut
    writeGlobalVar<int>("phifrm", pddm->iPhiFormatOut);
    // ibios     bBiosystemOut
    writeGlobalVar<bool>("ibios", pddm->bBiosystemOut);
    // ibem     bBemSrfOut
    writeGlobalVar<bool>("ibem", pddm->bBemSrfOut);
    // isite     bSiteOut
    writeGlobalVar<bool>("isite", pddm->bSiteOut);
    // frcfrm     iFrcFormatOut
    writeGlobalVar<int>("frcfrm", pddm->iFrcFormatOut);
    // epswrt     bEpsOut
    writeGlobalVar<bool>("epswrt", pddm->bEpsOut);
    // iatout     bModPdbOut
    writeGlobalVar<bool>("iatout", pddm->bModPdbOut);
    // mpdbfrm     iModPdbFormatOut
    writeGlobalVar<int>("mpdbfrm", pddm->iModPdbFormatOut);
    // ipdbwrt     bUnformatPdbOut
    writeGlobalVar<bool>("ipdbwrt", pddm->bUnformatPdbOut);
    // ifrcwrt     bUnformatFrcOut
    writeGlobalVar<bool>("ifrcwrt", pddm->bUnformatFrcOut);
    // inrgwrt     bEngOut
    writeGlobalVar<bool>("inrgwrt", pddm->bEngOut);
    // iwgcrg     bGridCrgOut
    writeGlobalVar<bool>("iwgcrg", pddm->bGridCrgOut);
    // iacs     bHsurf2DatOut
    writeGlobalVar<bool>("iacs", pddm->bHsurf2DatOut);
    // idbwrt     bDbOut
    writeGlobalVar<bool>("idbwrt", pddm->bDbOut);
    // isen     bSurfEngOut
    writeGlobalVar<bool>("isen", pddm->bSurfEngOut);
    // isch     bSurfCrgOut
    writeGlobalVar<bool>("isch", pddm->bSurfCrgOut);
    // scrgfrm     iSurfCrgFormatOut
    writeGlobalVar<int>("scrgfrm", pddm->iSurfCrgFormatOut);
    // logg     bGridEng
    writeGlobalVar<bool>("logg", pddm->bGridEng);
    // logs     bSolvEng
    writeGlobalVar<bool>("logs", pddm->bSolvEng);
    // logas     bAnalySurfEng
    writeGlobalVar<bool>("logas", pddm->bAnalySurfEng);
    // loga     bAnalyEng
    writeGlobalVar<bool>("loga", pddm->bAnalyEng);
    // logions     bIonsEng
    writeGlobalVar<bool>("logions", pddm->bIonsEng);
    // logc     bCoulombEng
    writeGlobalVar<bool>("logc", pddm->bCoulombEng);
    // isita     bAtomInSite
    writeGlobalVar<bool>("isita", pddm->bAtomInSite);
    // isitq     bCrgInSite
    writeGlobalVar<bool>("isitq", pddm->bCrgInSite);
    // isitp     bGridPotentialInSite
    writeGlobalVar<bool>("isitp", pddm->bGridPotentialInSite);
    // isitap     bAtomPotentialInSite
    writeGlobalVar<bool>("isitap", pddm->bAtomPotentialInSite);
    // isitdeb     bDebyeFractionInSite
    writeGlobalVar<bool>("isitdeb", pddm->bDebyeFractionInSite);
    // isitf     bFieldInSite
    writeGlobalVar<bool>("isitf", pddm->bFieldInSite);
    // isitr     bReactPotentialInSite
    writeGlobalVar<bool>("isitr", pddm->bReactPotentialInSite);
    // isitc     bCoulombPotentialInSite
    writeGlobalVar<bool>("isitc", pddm->bCoulombPotentialInSite);
    // isitx     bAtomCoordInSite
    writeGlobalVar<bool>("isitx", pddm->bAtomCoordInSite);
    // isiti     bSaltInSite
    writeGlobalVar<bool>("isiti", pddm->bSaltInSite);
    // isitt     bTotalPotentialInSite
    writeGlobalVar<bool>("isitt", pddm->bTotalPotentialInSite);
    // isitrf     bReactForceInSite
    writeGlobalVar<bool>("isitrf", pddm->bReactForceInSite);
    // isitcf     bCoulombForceInSite
    writeGlobalVar<bool>("isitcf", pddm->bCoulombForceInSite);
    // isitmd     bMDInSite
    writeGlobalVar<bool>("isitmd", pddm->bMDInSite);
    // isitsf     bSurfCrgInSite
    writeGlobalVar<bool>("isitsf", pddm->bSurfCrgInSite);
    // isittf     bTotalForceInSite
    writeGlobalVar<bool>("isittf", pddm->bTotalForceInSite);
    // isitpot     bPotentialInSite
    writeGlobalVar<bool>("isitpot", pddm->bPotentialInSite);
    // irea     bReactFieldInFRC
    writeGlobalVar<bool>("irea", pddm->bReactFieldInFRC);
    // iself     bPDB2FRCInSite
    writeGlobalVar<bool>("iself", pddm->bPDB2FRCInSite);
    // bufz     ieBuffz
    writeGlobalVar<SExtrema<delphi_integer>>("buffz", pddm->eiBuffz);
    // ibufz     bBuffz
    writeGlobalVar<bool>("ibufz", pddm->bIsBuffz);
    // isurftype     iTypeSurf
    writeGlobalVar<int>("isurftype", pddm->iTypeSurf);

    if (debug_data)cout << "*                 Delphi data broadcast functions done        *\n";


    //From Delphi Class
    // deblen     fDebyeLength
    writeGlobalVar<delphi_real>("deblen", pddm->fDebyeLength);
    // epsout     fEpsOut
    writeGlobalVar<delphi_real>("epsout", pddm->fEpsOut);
    // cran     fgCoordinateRange
    writeGlobalVar<SGrid<delphi_real>>("cran", pddm->gfCoordinateRange);
    // pmid     fgGeometricCenter
    writeGlobalVar<SGrid<delphi_real>>("pmid", pddm->gfGeometricCenter);
    // oldmid     fgBoxCenter
    writeGlobalVar<SGrid<delphi_real>>("oldmid", pddm->gfBoxCenter);
    // rionst     fIonStrength
    writeGlobalVar<delphi_real>("rionst", pddm->fIonStrength);
    // chi1     fTaylorCoeff1
    writeGlobalVar<delphi_real>("chi1", pddm->fTaylorCoeff1);
    // chi2     fTaylorCoeff2
    writeGlobalVar<delphi_real>("chi2", pddm->fTaylorCoeff2);
    // chi3     fTaylorCoeff3
    writeGlobalVar<delphi_real>("chi3", pddm->fTaylorCoeff3);
    // chi4     fTaylorCoeff4
    writeGlobalVar<delphi_real>("chi4", pddm->fTaylorCoeff4);
    // chi5     fTaylorCoeff5
    writeGlobalVar<delphi_real>("chi5", pddm->fTaylorCoeff5);
    // lognl     bNonlinearEng
    writeGlobalVar<bool>("lognl", pddm->bNonlinearEng);
    // epkt     fEPKT
    writeGlobalVar<delphi_real>("epkt", pddm->fEPKT);
    // epsin     fEpsIn
    writeGlobalVar<delphi_real>("epsin", pddm->fEpsIn);
    // ifrcrd     bFrcUnformatIn
    writeGlobalVar<bool>("ifrcrd", pddm->bFrcUnformatIn);
    // idirectalg     iDirectEpsMap
    writeGlobalVar<int>("idirectalg", pddm->iDirectEpsMap);
    // numbmol     iMoleculeNum
    writeGlobalVar<delphi_integer>("numbmol", pddm->iMoleculeNum);
    // rdmx     fMaxRadius
    writeGlobalVar<delphi_real>("rdmx", pddm->fMaxRadius);
    // uniformdiel     bUniformDielec
    writeGlobalVar<bool>("uniformdiel", pddm->bUniformDielec);
    // limobject(:)     prgfeExtrema
    writeGlobalVector1D<SExtrema<delphi_real>>("limobject", 0, pddm->vctefExtrema.size(), pddm->vctefExtrema);
    // xn1(:)     prgfgAtomCoordA
    writeGlobalVector1D<SGrid<delphi_real>>("xn1", 0, pddm->vctgfAtomCoordA.size(), pddm->vctgfAtomCoordA);
    // xn2(:)     prgfgAtomCoordG
    writeGlobalVector1D<SGrid<delphi_real>>("xn2", 0, pddm->vctgfAtomCoordG.size(), pddm->vctgfAtomCoordG);
    //"dataobject pddm->vctstrObject    
    writeGlobal1DString("dataobject", pddm->vctstrObject);

    if (debug_data)cout << "*                 Delphi data broadcast Delphi class done        *\n";

    // From IO Class
    //  myData["resnummax"] = pddm->iResidueNum;
    writeGlobalVar<delphi_integer>("resnummax", pddm->iResidueNum);
    //  myData["nmedia"] = pddm->iMediaNum;
    writeGlobalVar<delphi_integer>("nmedia", pddm->iMediaNum);
    //  myData["medeps"] = pddm->vctfMediaEps;     // std::vector<delphi_real>
    writeGlobalVector1D<delphi_real>("medeps", 0, pddm->vctfMediaEps.size(), pddm->vctfMediaEps);
    //  myData["nobject"] = pddm->iObjectNum;
    writeGlobalVar<delphi_integer>("nobject", pddm->iObjectNum);
    //  myData["dataobject"] = pddm->vctstrObject;     // std::vector<string> //no need to broadcast, used in space only
    //  myData["natom"] = pddm->iAtomNum;
    writeGlobalVar<delphi_integer>("natom", pddm->iAtomNum);
    //  myData["delphipdb"] = pddm->vctapAtomPdb;     // std::vector<CAtomPdb>
    writeGlobalVector1D<CAtomPdb>("delphipdb", 0, pddm->vctapAtomPdb.size(), pddm->vctapAtomPdb);
    //  myData["iatmmed"] = pddm->vctiAtomMediaNum; // std::vector<delphi_integer>
    writeGlobalVector1D<delphi_integer>("iatmmed", 0, pddm->vctiAtomMediaNum.size(), pddm->vctiAtomMediaNum);
    //  myData["ionlymol"] = pddm->bOnlyMolecule;
    writeGlobalVar<bool>("ionlymol", pddm->bOnlyMolecule);

    if (debug_data)cout << "*                 Delphi data broadcast IO class done        *\n";
}

void CDelphiData::sync_pre_space()
{
    if (debug_data)cout << "*      Sync before Space Module        *\n";
    //iGrid(pdc->getKey_Ref<delphi_integer>("igrid"))
    getKey_Ref<delphi_integer>("igrid") = readGlobalVar<delphi_integer>("igrid");
    //checkValue<bool>("iautocon");
    getKey_Ref<bool>("iautocon") = readGlobalVar<bool>("iautocon");
    //iNatom(getKey_Val<delphi_integer>("natom")),
    getKey_Ref<delphi_integer>("natom") = readGlobalVar<delphi_integer>("natom");
    //fScale(getKey_constRef<delphi_real>("scale")),
    getKey_Ref<delphi_real>("scale") = readGlobalVar<delphi_real>("scale");
    //iNObject(getKey_constRef<delphi_integer>("nobject")),
    getKey_Ref<delphi_integer>("nobject") = readGlobalVar<delphi_integer>("nobject");
    //repsout(getKey_constRef<delphi_real>("repsout")),
    getKey_Ref<delphi_real>("repsout") = readGlobalVar<delphi_real>("repsout");
    //repsin(getKey_constRef<delphi_real>("repsin")),
    getKey_Ref<delphi_real>("repsin") = readGlobalVar<delphi_real>("repsin");
    ////ndistr (getKey_constRef<delphi_integer>("ndistr")),
    //getKey_Ref<delphi_integer>("ndistr") = readGlobalVar<delphi_integer>("ndistr");
    //bUniformDiel(getKey_constRef<bool>("uniformdiel")),
    getKey_Ref<bool>("uniformdiel") = readGlobalVar<bool>("uniformdiel");
    //bOnlyMol(getKey_constRef<bool>("ionlymol")),
    getKey_Ref<bool>("ionlymol") = readGlobalVar<bool>("ionlymol");
    //isolv(getKey_constRef<bool>("isolv")),
    getKey_Ref<bool>("isolv") = readGlobalVar<bool>("isolv");
    //irea(getKey_constRef<bool>("irea")),
    getKey_Ref<bool>("irea") = readGlobalVar<bool>("irea");
    //logs(getKey_constRef<bool>("logs")),
    getKey_Ref<bool>("logs") = readGlobalVar<bool>("logs");
    //lognl(getKey_constRef<bool>("lognl")),
    getKey_Ref<bool>("lognl") = readGlobalVar<bool>("lognl");
    //isen(getKey_constRef<bool>("isen")),
    getKey_Ref<bool>("isen") = readGlobalVar<bool>("isen");
    //isch(getKey_constRef<bool>("isch")),
    getKey_Ref<bool>("isch") = readGlobalVar<bool>("isch");
    //fEpsOut(getKey_constRef<delphi_real>("epsout")),
    getKey_Ref<delphi_real>("epsout") = readGlobalVar<delphi_real>("epsout");
    //fDebyeLength(getKey_constRef<delphi_real>("deblen")),
    getKey_Ref<delphi_real>("deblen") = readGlobalVar<delphi_real>("deblen");
    //fEpsIn(getKey_constRef<delphi_real>("epsin")),
    getKey_Ref<delphi_real>("epsin") = readGlobalVar<delphi_real>("epsin");

    ////Argo: for EPSMAP in conjunction with CONVOLUTION
    //bEpsOut(getKey_constRef<bool>("epswrt")),
    getKey_Ref<bool>("epswrt") = readGlobalVar<bool>("epswrt");
    //strEpsFile(getKey_constRef<string>("epsnam")),
    getKey_Ref<string>("epsnam") = readGlobalChar("epsnam");
    ////ARGO-Putting in the reference for fgBoxCenter
    //fgBoxCenter(getKey_constRef< SGrid<delphi_real> >("oldmid")),
    getKey_Ref< SGrid<delphi_real>>("oldmid") = readGlobalVar< SGrid<delphi_real>>("oldmid");
    ////
    //isite(getKey_constRef<bool>("isite")),
    getKey_Ref<bool>("isite") = readGlobalVar<bool>("isite");
    //ibem(getKey_constRef<bool>("ibem")),
    getKey_Ref<bool>("ibem") = readGlobalVar<bool>("ibem");
    //ibctyp(getKey_constRef<int>("ibctyp")),
    getKey_Ref<delphi_integer>("ibctyp") = readGlobalVar<delphi_integer>("ibctyp");
    //isitsf(getKey_constRef<bool>("isitsf")),
    getKey_Ref<bool>("isitsf") = readGlobalVar<bool>("isitsf");

    //// Lin Li : Gaussian:
    //cutoff(getKey_constRef<float>("cutoff")),
    getKey_Ref<delphi_real>("cutoff") = readGlobalVar<delphi_real>("cutoff");
    //sigma(getKey_constRef<float>("sigma")),
    getKey_Ref<delphi_real>("sigma") = readGlobalVar<delphi_real>("sigma");
    //inhomo(getKey_Ref<int>("inhomo")),
    getKey_Ref<delphi_integer>("inhomo") = readGlobalVar<delphi_integer>("inhomo");
    //srfcut(getKey_constRef<float>("srfcut")),
    getKey_Ref<delphi_real>("srfcut") = readGlobalVar<delphi_real>("srfcut");
    //iGaussian(getKey_constRef<int>("gaussian")),
    getKey_Ref<delphi_integer>("gaussian") = readGlobalVar<delphi_integer>("gaussian");

    ////bDebug (getKey_constRef<bool>("debug")),

    ////iTestGloble (getKey_constRef<int>(" ")),
    //iNMedia(getKey_constRef<delphi_integer>("nmedia")),
    getKey_Ref<delphi_integer>("nmedia") = readGlobalVar<delphi_integer>("nmedia");
    //numbmol(getKey_constRef<delphi_integer>("numbmol")),
    getKey_Ref<delphi_integer>("numbmol") = readGlobalVar<delphi_integer>("numbmol");
    //scrgfrm(getKey_constRef<delphi_integer>("scrgfrm")),
    getKey_Ref<delphi_integer>("scrgfrm") = readGlobalVar<delphi_integer>("scrgfrm");

    ////ndistr (getKey_constRef<delphi_integer>("ndistr")),

    ////fRMid (getKey_constRef<delphi_real>("rmid")),
    //global_cOldMid(getKey_Val< SGrid<delphi_real> >("oldmid")),

    //fIonStrenth(getKey_constRef<delphi_real>("rionst")),
    getKey_Ref<delphi_real>("rionst") = readGlobalVar<delphi_real>("rionst");
    //fExternRadius(getKey_constRef<delphi_real>("exrad")),
    getKey_Ref<delphi_real>("exrad") = readGlobalVar<delphi_real>("exrad");
    //fRMax(getKey_constRef<delphi_real>("rdmx")),
    getKey_Ref<delphi_real>("rdmx") = readGlobalVar<delphi_real>("rdmx");
    //fRadPrb_v(getKey_Ref<vector <delphi_real> >("radprb")),
    readGlobalVector1D<delphi_real>("radprb", 0, sizeofGlobal1D<delphi_real>("radprb"), getKey_Ref<vector <delphi_real> >("radprb"));
    //dataobject_v(getKey_Ref< vector<string> >("dataobject")),

    if (debug_data)
    {
        for (int i = 0; i < getKey_Ref< vector<string> >("dataobject").size(); i++)
            cout << "*  before sync     *\n" << getKey_Ref< vector<string> >("dataobject")[i];
    }

    readGlobal1DString("dataobject", getKey_Ref< vector<string> >("dataobject"));

    if (debug_data)
    {
        for (int i = 0; i < getKey_Ref< vector<string> >("dataobject").size(); i++)
            cout << "*  after sync     *\n" << getKey_Ref< vector<string> >("dataobject")[i];
    }
    ////++++++++++++++++ reference to read-and-write variables from data container +++++++++++++++//

    //delphipdb(getKey_Ref<vector <CAtomPdb> >("delphipdb")),
    //cout << " on node " << myid() << " reading delphipdb size = " << sizeofGlobal1D<CAtomPdb>("delphipdb") << endl; 
    readGlobalVector1D<CAtomPdb>("delphipdb", 0, sizeofGlobal1D<CAtomPdb>("delphipdb"), getKey_Ref<vector <CAtomPdb> >("delphipdb"));
    //iBoundNum(getKey_Ref<delphi_integer>("ibnum")),
    getKey_Ref<delphi_integer>("ibnum") = readGlobalVar<delphi_integer>("ibnum");
    //rdmx(getKey_Ref<delphi_real>("rdmx")),
    getKey_Ref<delphi_real>("rdmx") = readGlobalVar<delphi_real>("rdmx");
    //nqass(getKey_Ref<delphi_integer>("nqass")),
    getKey_Ref<delphi_integer>("nqass") = readGlobalVar<delphi_integer>("nqass");
    //nqgrd(getKey_Ref<delphi_integer>("nqgrd")),
    getKey_Ref<delphi_integer>("nqgrd") = readGlobalVar<delphi_integer>("nqgrd");
    //iacs(getKey_Ref< bool >("iacs")),
    getKey_Ref<bool>("iacs") = readGlobalVar<bool>("iacs");
    //isrf(getKey_Ref< bool >("isrf")),
    getKey_Ref<bool>("isrf") = readGlobalVar<bool>("isrf");
    //cMin(getKey_Ref< SGrid<delphi_real> >("cmin")),
    getKey_Ref<SGrid<delphi_real>>("cmin") = readGlobalVar<SGrid<delphi_real>>("cmin");
    //cMax(getKey_Ref< SGrid<delphi_real> >("cmax")),
    getKey_Ref<SGrid<delphi_real>>("cmax") = readGlobalVar<SGrid<delphi_real>>("cmax");
    //qnet(getKey_Ref< delphi_real >("qnet")),
    getKey_Ref<delphi_real>("qnet") = readGlobalVar<delphi_real>("qnet");
    //qmin(getKey_Ref< delphi_real >("qmin")),
    getKey_Ref<delphi_real>("qmin") = readGlobalVar<delphi_real>("qmin");
    //qplus(getKey_Ref< delphi_real >("qplus")),
    getKey_Ref<delphi_real>("qplus") = readGlobalVar<delphi_real>("qplus");
    //cqmin(getKey_Ref< SGrid<delphi_real> >("cqmin")),
    getKey_Ref<SGrid<delphi_real>>("cqmin") = readGlobalVar<SGrid<delphi_real>>("cqmin");
    //cqplus(getKey_Ref< SGrid<delphi_real> >("cqplus")),
    getKey_Ref<SGrid<delphi_real>>("cqplus") = readGlobalVar<SGrid<delphi_real>>("cqplus");
    //acenter(getKey_Ref< SGrid<delphi_real> >("acent")),
    getKey_Ref<SGrid<delphi_real>>("acent") = readGlobalVar<SGrid<delphi_real>>("acent");
    //medeps(getKey_Ref < vector < delphi_real > >("medeps")),
    readGlobalVector1D<delphi_real>("medeps", 0, sizeofGlobal1D<delphi_real>("medeps"), getKey_Ref < vector < delphi_real > >("medeps"));
    //xn1_v(getKey_Ref< vector< SGrid<delphi_real> > >("xn1")),
    readGlobalVector1D<SGrid<delphi_real> >("xn1", 0, sizeofGlobal1D<SGrid<delphi_real> >("xn1"), getKey_Ref < vector < SGrid<delphi_real>  > >("xn1"));
    //xn2_v(getKey_Ref< vector< SGrid<delphi_real> > >("xn2")),
    readGlobalVector1D<SGrid<delphi_real> >("xn2", 0, sizeofGlobal1D<SGrid<delphi_real> >("xn2"), getKey_Ref < vector < SGrid<delphi_real>  > >("xn2"));
    //bDebMap_v(getKey_Ref< vector< char > >("idebmap")),
    readGlobalVector1D<char>("idebmap", 0, sizeofGlobal1D<char>("idebmap"), getKey_Ref< vector< char > >("idebmap"));

    ////ARGO: doing the 'pdc' thing just like its done with idebmap
    //zetaSurfMap_v(getKey_Ref< vector< char > >("zetaSurfMap")),
    readGlobalVector1D<char>("zetaSurfMap", 0, sizeofGlobal1D<char>("zetaSurfMap"), getKey_Ref< vector< char > >("zetaSurfMap"));
    //zetaOn(getKey_Ref<int>("zetaOn")),
    getKey_Ref<int>("zetaOn") = readGlobalVar<int>("zetaOn");
    //zetaDistance(getKey_Ref<delphi_real>("zetaDistance")),
    getKey_Ref<delphi_real>("zetaDistance") = readGlobalVar<delphi_real>("zetaDistance");
    //strZetaPhiFile(getKey_constRef<string>("zphinam")),
    getKey_Ref<string>("zphinam") = readGlobalChar("zphinam");

    //iEpsMap_v(getKey_Ref< vector< SGrid<delphi_integer> > >("iepsmp")),
    readGlobalVector1D<SGrid<delphi_integer>>("iepsmp", 0, sizeofGlobal1D<SGrid<delphi_integer>>("iepsmp"), getKey_Ref< vector< SGrid<delphi_integer> > >("iepsmp"));
    //fGepsMap_v(getKey_Ref< vector< SGrid<delphi_real> > >("gepsmp")),
    readGlobalVector1D<SGrid<delphi_real>>("gepsmp", 0, sizeofGlobal1D<SGrid<delphi_real>>("gepsmp"), getKey_Ref< vector< SGrid<delphi_real> > >("gepsmp"));
    //fGepsMap2_v(getKey_Ref< vector< SGrid<delphi_real> > >("gepsmp2")),
    readGlobalVector1D<SGrid<delphi_real>>("gepsmp2", 0, sizeofGlobal1D<SGrid<delphi_real>>("gepsmp2"), getKey_Ref< vector< SGrid<delphi_real> > >("gepsmp2"));
    //fGDensityMap_v(getKey_Ref< vector< delphi_real> >("gdensity")),
    readGlobalVector1D<delphi_real>("gdensity", 0, sizeofGlobal1D<delphi_real>("gdensity"), getKey_Ref< vector< delphi_real> >("gdensity"));
    //iAtomMed_v(getKey_Ref< vector<delphi_integer> >("iatmmed")),
    readGlobalVector1D<delphi_integer>("iatmmed", 0, sizeofGlobal1D<delphi_integer>("iatmmed"), getKey_Ref< vector<delphi_integer> >("iatmmed"));
    //sLimObject(getKey_Val< vector < SExtrema<delphi_real> > >("limobject")),
    readGlobalVector1D<SExtrema<delphi_real>>("limobject", 0, sizeofGlobal1D<SExtrema<delphi_real>>("limobject"), getKey_Ref< vector < SExtrema<delphi_real> > >("limobject"));
    //ibgrd_v(getKey_Ref< vector< SGrid<delphi_integer> > >("ibgrd")),
    readGlobalVector1D<SGrid<delphi_integer>>("ibgrd", 0, sizeofGlobal1D<SGrid<delphi_integer>>("ibgrd"), getKey_Ref< vector< SGrid<delphi_integer> > >("ibgrd"));
    //scspos_v(getKey_Ref< vector< SGrid<delphi_real> > >("scspos")),
    readGlobalVector1D<SGrid<delphi_real>>("scspos", 0, sizeofGlobal1D<SGrid<delphi_real>>("scspos"), getKey_Ref< vector< SGrid<delphi_real> > >("scspos"));
    //chgpos_v(getKey_Ref< vector< SGrid<delphi_real> > >("chgpos")),
    readGlobalVector1D<SGrid<delphi_real>>("chgpos", 0, sizeofGlobal1D<SGrid<delphi_real>>("chgpos"), getKey_Ref< vector< SGrid<delphi_real> > >("chgpos"));
    //crgatn_v(getKey_Ref< vector< delphi_integer > >("crgatn")),
    readGlobalVector1D<delphi_integer>("crgatn", 0, sizeofGlobal1D<delphi_integer>("crgatn"), getKey_Ref< vector< delphi_integer > >("crgatn"));
    //nqgrdtonqass_v(getKey_Ref< vector< delphi_integer > >("nqgrdtonqass")),
    readGlobalVector1D<delphi_integer>("nqgrdtonqass", 0, sizeofGlobal1D<delphi_integer>("nqgrdtonqass"), getKey_Ref< vector< delphi_integer > >("nqgrdtonqass"));
    //atmeps_v(getKey_Ref< vector< delphi_real > >("atmeps")),
    readGlobalVector1D<delphi_real>("atmeps", 0, sizeofGlobal1D<delphi_real>("atmeps"), getKey_Ref< vector< delphi_real > >("atmeps"));
    //atmcrg_v(getKey_Ref< vector< SGridValue<delphi_real> > >("atmcrg")),
    readGlobalVector1D<SGridValue<delphi_real>>("atmcrg", 0, sizeofGlobal1D<SGridValue<delphi_real>>("atmcrg"), getKey_Ref< vector< SGridValue<delphi_real> > >("atmcrg"));
    //chrgv2_v(getKey_Ref< vector< SGridValue<delphi_real> > >("chrgv2")),
    readGlobalVector1D<SGridValue<delphi_real>>("chrgv2", 0, sizeofGlobal1D<SGridValue<delphi_real>>("chrgv2"), getKey_Ref< vector< SGridValue<delphi_real> > >("chrgv2"));

    //scsnor_v(getKey_Ref< vector< SGrid<delphi_real> > >("scsnor")),
    readGlobalVector1D<SGrid<delphi_real> >("scsnor", 0, sizeofGlobal1D<SGrid<delphi_real> >("scsnor"), getKey_Ref< vector< SGrid<delphi_real> > >("scsnor"));
    //atsurf_v(getKey_Ref< vector< delphi_integer > >("atsurf")),
    readGlobalVector1D<delphi_integer>("atsurf", 0, sizeofGlobal1D<delphi_integer>("atsurf"), getKey_Ref< vector< delphi_integer > >("atsurf"));
    //atndx_v(getKey_Ref< vector< delphi_integer > >("atndx"))
    readGlobalVector1D<delphi_integer>("atndx", 0, sizeofGlobal1D<delphi_integer>("atndx"), getKey_Ref< vector< delphi_integer > >("atndx"));
    if (debug_data)std::cout << "*    Sync_pre_space on node " << myid() << " done   *\n";
}

void CDelphiData::sync_pre_solver()
{
    //----- uniform parameters
    getKey_Ref<string>("biomodel")          = readGlobalChar("biomodel");
    getKey_Ref<string>("solver")            = readGlobalChar("solver");

    //----- set by Statements
    getKey_Ref<delphi_integer>("natom")     = readGlobalVar<delphi_integer>("natom");
    getKey_Ref<bool>("isph")                = readGlobalVar<bool>("isph");
    getKey_Ref<bool>("iuspec")              = readGlobalVar<bool>("iuspec");
    getKey_Ref<delphi_real>("uspec")        = readGlobalVar<delphi_real>("uspec");
    readGlobalVector1D<char>("iper", 0, sizeofGlobal1D<char>("iper"), getKey_Ref< vector<char> >("iper"));
    getKey_Ref<int>("nnit")                 = readGlobalVar<int>("nnit");
    getKey_Ref<delphi_real>("exrad")        = readGlobalVar<delphi_real>("exrad");
    getKey_Ref<bool>("icheb")               = readGlobalVar<bool>("icheb");
    getKey_Ref<bool>("iautocon")            = readGlobalVar<bool>("iautocon");
    getKey_Ref<delphi_real>("gten")         = readGlobalVar<delphi_real>("gten");
    getKey_Ref<delphi_real>("res1")         = readGlobalVar<delphi_real>("res1");
    getKey_Ref<delphi_real>("res2")         = readGlobalVar<delphi_real>("res2");
    getKey_Ref<bool>("igraph")              = readGlobalVar<bool>("igraph");
    getKey_Ref<bool>("ipoten")              = readGlobalVar<bool>("ipoten");
    getKey_Ref<bool>("imanual")             = readGlobalVar<bool>("imanual");
    getKey_Ref<int>("phiintype")            = readGlobalVar<int>("phiintype");
    getKey_Ref<delphi_real>("repsout")      = readGlobalVar<delphi_real>("repsout");
    getKey_Ref<delphi_real>("repsin")       = readGlobalVar<delphi_real>("repsin");

    //----- io file names
    getKey_Ref<string>("epsnam")            = readGlobalChar("epsnam");
    getKey_Ref<string>("phiinam")           = readGlobalChar("phiinam");

    //----- set by functions
    getKey_Ref<bool>("iwgcrg")              = readGlobalVar<bool>("iwgcrg");
    getKey_Ref<bool>("epswrt")              = readGlobalVar<bool>("epswrt");
    getKey_Ref<bool>("logions")             = readGlobalVar<bool>("logions");
    getKey_Ref<int>("gaussian")             = readGlobalVar<int>("gaussian");
    getKey_Ref<int>("inhomo")               = readGlobalVar<int>("inhomo");

    //----- set by DelPhi
    getKey_Ref<delphi_real>("epsout")       = readGlobalVar<delphi_real>("epsout");
    getKey_Ref<delphi_real>("deblen")       = readGlobalVar<delphi_real>("deblen");
    getKey_Ref<delphi_real>("scale")        = readGlobalVar<delphi_real>("scale");
    getKey_Ref<delphi_real>("epsin")        = readGlobalVar<delphi_real>("epsin");
    getKey_Ref<delphi_real>("rionst")       = readGlobalVar<delphi_real>("rionst");
    getKey_Ref<int>("idirectalg")           = readGlobalVar<int>("idirectalg");
    getKey_Ref<delphi_real>("epkt")         = readGlobalVar<delphi_real>("epkt");
    getKey_Ref<SGrid<delphi_real>>("oldmid")= readGlobalVar<SGrid<delphi_real>>("oldmid");
    getKey_Ref<SGrid<delphi_real>>("vdrop") = readGlobalVar<SGrid<delphi_real>>("vdrop");
    getKey_Ref<delphi_real>("chi2")         = readGlobalVar<delphi_real>("chi2");
    getKey_Ref<delphi_real>("chi3")         = readGlobalVar<delphi_real>("chi3");
    getKey_Ref<delphi_real>("chi4")         = readGlobalVar<delphi_real>("chi4");
    getKey_Ref<delphi_real>("chi5")         = readGlobalVar<delphi_real>("chi5");
    getKey_Ref<bool>("uniformdiel")         = readGlobalVar<bool>("uniformdiel");
    getKey_Ref<delphi_real>("srfcut")       = readGlobalVar<delphi_real>("srfcut");

    //----- set by IO class
    getKey_Ref<delphi_integer>("nmedia")    = readGlobalVar<delphi_integer>("nmedia");
    getKey_Ref<delphi_integer>("nobject")   = readGlobalVar<delphi_integer>("nobject");

    /*
     * solver does not sync medeps(0:nmediamax) to save memory usage
     */
    //readGlobalVector1D<delphi_real>("medeps", 0, sizeofGlobal1D<delphi_real>("medeps"), getKey_Ref< vector<delphi_real> >("medeps"));

    //----- set by Surface class
    getKey_Ref<delphi_integer>("nqass")     = readGlobalVar<delphi_integer>("nqass");
    getKey_Ref<delphi_real>("qmin")         = readGlobalVar<delphi_real>("qmin");
    getKey_Ref<delphi_real>("qplus")        = readGlobalVar<delphi_real>("qplus");
    getKey_Ref<SGrid<delphi_real>>("cqplus")= readGlobalVar<SGrid<delphi_real>>("cqplus");
    getKey_Ref<SGrid<delphi_real>>("cqmin") = readGlobalVar<SGrid<delphi_real>>("cqmin");
    getKey_Ref<delphi_integer>("ibnum")     = readGlobalVar<delphi_integer>("ibnum");

    //readGlobalVector1D<SGrid<delphi_integer>>("ibgrd", 0, sizeofGlobal1D<SGrid<delphi_integer>>("ibgrd"), getKey_Ref< vector< SGrid<delphi_integer> > >("ibgrd"));

    getKey_Ref<delphi_integer>("nqgrd")     = readGlobalVar<delphi_integer>("nqgrd");

    //readGlobalVector1D<SGridValue<delphi_real>>("chrgv2", 0, sizeofGlobal1D<SGridValue<delphi_real>>("chrgv2"), getKey_Ref< vector< SGridValue<delphi_real> > >("chrgv2"));

    //readGlobalVector1D<delphi_integer>("nqgrdtonqass", 0, sizeofGlobal1D<delphi_integer>("nqgrdtonqass"), getKey_Ref< vector<delphi_integer> >("nqgrdtonqass"));

    //readGlobalVector1D<delphi_real>("atmeps", 0, sizeofGlobal1D<delphi_real>("atmeps"), getKey_Ref< vector<delphi_real> >("atmeps"));

    //readGlobalVector1D<SGridValue<delphi_real>>("atmcrg", 0, sizeofGlobal1D<SGridValue<delphi_real>>("atmcrg"), getKey_Ref< vector<SGridValue<delphi_real> > >("atmcrg"));

    //readGlobalVector1D<SGrid<delphi_real>>("scspos", 0, sizeofGlobal1D<SGrid<delphi_real>>("scspos"), prgfgSurfCrgA);

    //readGlobalVector1D<SGrid<delphi_integer>>("iepsmp", 0, sizeofGlobal1D<SGrid<delphi_integer>>("iepsmp"), getKey_Ref< vector< SGrid<delphi_integer> > >("iepsmp"));

    //readGlobalVector1D<SGrid<delphi_real>>("gepsmp2", 0, sizeofGlobal1D<SGrid<delphi_real>>("gepsmp2"), getKey_Ref< vector< SGrid<delphi_real> > >("gepsmp2"));

    //getKey_Ref<vector< char>>("idebmap").resize(sizeofGlobal1D<char>("idebmap"));
    //readGlobalVector1D<char>("idebmap", 0, sizeofGlobal1D<char>("idebmap"), getKey_Ref<vector< char>>("idebmap"));

    //readGlobalVector1D<delphi_real>("gdensity", 0, sizeofGlobal1D<delphi_real>("gdensity"), getKey_Ref<vector<delphi_real>>("gdensity"));

    //-----reference to read-and-write variables from data container
    getKey_Ref<delphi_integer>("igrid")     = readGlobalVar<delphi_integer>("igrid");
    getKey_Ref<delphi_integer>("ibctyp")    = readGlobalVar<delphi_integer>("ibctyp");
    getKey_Ref<delphi_real>("qnet")         = readGlobalVar<delphi_real>("qnet");
    getKey_Ref<int>("nlit")                 = readGlobalVar<int>("nlit");
    getKey_Ref<int>("icon1")                = readGlobalVar<int>("icon1");
    getKey_Ref<int>("icon2")                = readGlobalVar<int>("icon2");
    getKey_Ref<bool>("idbwrt")              = readGlobalVar<bool>("idbwrt");
    getKey_Ref<delphi_real>("relpar")       = readGlobalVar<delphi_real>("relpar");

    //----- out into data container
    getKey_Ref<delphi_integer>("icount2b")  = readGlobalVar<delphi_integer>("icount2b");
    getKey_Ref<delphi_integer>("icount1b")  = readGlobalVar<delphi_integer>("icount1b");

    //readGlobalVector1D<delphi_real>("gchrg", 0, sizeofGlobal1D<delphi_real>("gchrg"), getKey_Ref< vector<delphi_real> >("gchrg"));

    //readGlobalVector1D<SGrid<delphi_integer>>("gchrgp", 0, sizeofGlobal1D<SGrid<delphi_integer>>("gchrgp"), getKey_Ref< vector< SGrid<delphi_integer> > >("gchrgp"));

    getKey_Ref<delphi_integer>("ibc")       = readGlobalVar<delphi_integer>("ibc");

    //readGlobalVector1D<SDoubleGridValue>("cgbp", 0, sizeofGlobal1D<SDoubleGridValue>("cgbp"), getKey_Ref< vector<SDoubleGridValue> >("cgbp"));

    //getKey_Ref< vector<delphi_real> >("phimap").resize(sizeofGlobal1D<delphi_real>("phimap"));
    //readGlobalVector1D<delphi_real>("phimap", 0, sizeofGlobal1D<delphi_real>("phimap"), getKey_Ref< vector<delphi_real> >("phimap"));

    //getKey_Ref< vector<delphi_real> >("phimap_pre").resize(sizeofGlobal1D<delphi_real>("phimap_pre"));
    //readGlobalVector1D<delphi_real>("phimap_pre", 0, sizeofGlobal1D<delphi_real>("phimap_pre"), getKey_Ref< vector<delphi_real> >("phimap_pre"));
}

void CDelphiData::sync_pre_energy()
{
    //iBndyType(getKey_constRef<int>("ibctyp")),
    getKey_Ref<delphi_integer>("ibctyp") = readGlobalVar<delphi_integer>("ibctyp");

    //strEnergyFile(getKey_constRef<string>("nrgnam")),
    getKey_Ref<string>("nrgnam") = readGlobalChar("nrgnam");

    //strScrgFile(getKey_constRef<string>("scrgnam")),
    getKey_Ref<string>("scrgnam") = readGlobalChar("scrgnam");

    //bEngOut(getKey_constRef<bool>("inrgwrt")),
    getKey_Ref<bool>("inrgwrt") = readGlobalVar<bool>("inrgwrt");

    //bSurfEngOut(getKey_constRef<bool>("isen")),
    getKey_Ref<bool>("isen") = readGlobalVar<bool>("isen");

    //bSurfCrgOut(getKey_constRef<bool>("isch")),
    getKey_Ref<bool>("isch") = readGlobalVar<bool>("isch");

    //iSurfCrgFormatOut(getKey_constRef<int>("scrgfrm")),
    getKey_Ref<delphi_integer>("scrgfrm") = readGlobalVar<delphi_integer>("scrgfrm");

    //bGridEng(getKey_constRef<bool>("logg")),
    getKey_Ref<bool>("logg") = readGlobalVar<bool>("logg");

    //bSolvEng(getKey_constRef<bool>("logs")),
    getKey_Ref<bool>("logs") = readGlobalVar<bool>("logs");

    //bAnalySurfEng(getKey_constRef<bool>("logas")),
    getKey_Ref<bool>("logas") = readGlobalVar<bool>("logas");

    //bAnalyEng(getKey_constRef<bool>("loga")),
    getKey_Ref<bool>("loga") = readGlobalVar<bool>("loga");

    //bIonsEng(getKey_constRef<bool>("logions")),
    getKey_Ref<bool>("logions") = readGlobalVar<bool>("logions");

    //bCoulombEng(getKey_constRef<bool>("logc")),
    getKey_Ref<bool>("logc") = readGlobalVar<bool>("logc");

    //bPotentiallnSite(getKey_constRef<bool>("isitpot")),
    getKey_Ref<bool>("isitpot") = readGlobalVar<bool>("isitpot");

    //bReactFieldlnFRC(getKey_constRef<bool>("irea")),
    getKey_Ref<bool>("irea") = readGlobalVar<bool>("irea");

    //bBuffz(getKey_constRef<bool>("ibufz")),
    getKey_Ref<bool>("ibufz") = readGlobalVar<bool>("ibufz");

    //iCrgGridNum(getKey_constRef<delphi_integer>("nqass")),
    getKey_Ref<delphi_integer>("nqass") = readGlobalVar<delphi_integer>("nqass");

    //fEpsOut(getKey_constRef<delphi_real>("epsout")),
    getKey_Ref<delphi_real>("epsout") = readGlobalVar<delphi_real>("epsout");

    //bNonlinearEng(getKey_constRef<bool>("lognl")),
    getKey_Ref<bool>("lognl") = readGlobalVar<bool>("lognl");

    //iMediaNum(getKey_constRef<delphi_integer>("nmedia")),
    getKey_Ref<delphi_integer>("nmedia") = readGlobalVar<delphi_integer>("nmedia");

    //iGaussian(getKey_constRef<int>("gaussian")),
    getKey_Ref<delphi_integer>("gaussian") = readGlobalVar<delphi_integer>("gaussian");

    //inhomo(getKey_constRef<int>("inhomo")),
    getKey_Ref<delphi_integer>("inhomo") = readGlobalVar<delphi_integer>("inhomo");

    //prgrMediaEps(getKey_constRef<delphi_real>("medeps")),
    //getKey_Ref<delphi_real>("medeps") = readGlobalVar<delphi_real>("medeps");
    readGlobalVector1D<delphi_real>("medeps", 0, sizeofGlobal1D<delphi_real>("medeps"), getKey_Ref < vector < delphi_real > >("medeps"));

    //prgiAtomMediaNum(getKey_constRef<delphi_integer>("iatmmed")),
    //getKey_Ref<delphi_integer>("iatmmed") = readGlobalVar<delphi_integer>("iatmmed");
    readGlobalVector1D<delphi_integer>("iatmmed", 0, sizeofGlobal1D<delphi_integer>("iatmmed"), getKey_Ref< vector<delphi_integer> >("iatmmed"));

    //iTotalBdyGridNum(getKey_constRef<delphi_integer>("ibnum")),
    getKey_Ref<delphi_integer>("ibnum") = readGlobalVar<delphi_integer>("ibnum");

    //iCrgBdyGrid(getKey_constRef<delphi_integer>("ibc")),
    getKey_Ref<delphi_integer>("ibc") = readGlobalVar<delphi_integer>("ibc");

    //prgigBndyGrid(getKey_constRef< vector< SGrid<int> > >("ibgrd")),
    readGlobalVector1D<SGrid<delphi_integer>>("ibgrd", 0, sizeofGlobal1D<SGrid<delphi_integer>>("ibgrd"), getKey_Ref< vector< SGrid<delphi_integer> > >("ibgrd"));

    //iGrid(getKey_constRef<delphi_integer>("igrid")),
    getKey_Ref<delphi_integer>("igrid") = readGlobalVar<delphi_integer>("igrid");

    //iCrgedGridB(getKey_constRef<delphi_integer>("icount1b")),
    getKey_Ref<delphi_integer>("icount1b") = readGlobalVar<delphi_integer>("icount1b");

    //fEpsin(getKey_constRef<delphi_real>("epsin")),
    getKey_Ref<delphi_real>("epsin") = readGlobalVar<delphi_real>("epsin");

    ////        DEVELOPER(getKey_constRef<bool>("ideveloper")),
    //prgfPhimap(getKey_constRef< vector<delphi_real> >("phimap")),  // phimap read function to convert to 3d array. //
    readGlobalVector1D<delphi_real>("phimap", 0, sizeofGlobal1D<delphi_real>("phimap"), getKey_Ref< vector<delphi_real> >("phimap"));

    //fScale(getKey_constRef<delphi_real>("scale")),
    getKey_Ref<delphi_real>("scale") = readGlobalVar<delphi_real>("scale");

    //fgPotentialDrop(getKey_constRef< SGrid<delphi_real> >("vdrop")),
    getKey_Ref<SGrid<delphi_real>>("vdrop") = readGlobalVar<SGrid<delphi_real>>("vdrop");

    //prgfgCrgPoseA(getKey_constRef< vector< SGrid<delphi_real> > >("chgpos")),
    readGlobalVector1D<SGrid<delphi_real>>("chgpos", 0, sizeofGlobal1D<SGrid<delphi_real>>("chgpos"), getKey_Ref< vector< SGrid<delphi_real> > >("chgpos"));

    //prgfGridCrg(getKey_Ref< vector<delphi_real> >("gchrg")),
    readGlobalVector1D<delphi_real>("gchrg", 0, sizeofGlobal1D<delphi_real>("gchrg"), getKey_Ref< vector<delphi_real> >("gchrg"));

    //prgigGridCrgPose(getKey_constRef< vector< SGrid<int> > >("gchrgp")),
    readGlobalVector1D<SGrid<delphi_integer>>("gchrgp", 0, sizeofGlobal1D<SGrid<delphi_integer>>("gchrgp"), getKey_Ref< vector< SGrid<delphi_integer> > >("gchrgp"));

    //prggvAtomicCrg(getKey_constRef< vector< SGridValue<delphi_real> > >("atmcrg")),
    readGlobalVector1D<SGridValue<delphi_real>>("atmcrg", 0, sizeofGlobal1D<SGridValue<delphi_real>>("atmcrg"), getKey_Ref< vector< SGridValue<delphi_real> > >("atmcrg"));

    //fIonStrength(getKey_Val<delphi_real>("rionst")),
    getKey_Ref<delphi_real>("rionst") = readGlobalVar<delphi_real>("rionst");

    //fgBoxCenter(getKey_constRef< SGrid<delphi_real> >("oldmid")),
    getKey_Ref<SGrid<delphi_real>>("oldmid") = readGlobalVar<SGrid<delphi_real>>("oldmid");

    //prgbDielecMap(getKey_Ref< vector<char> >("idebmap")),
    readGlobalVector1D<char>("idebmap", 0, sizeofGlobal1D<char>("idebmap"), getKey_Ref< vector<char> >("idebmap"));

    //prgfAtomEps(getKey_constRef< vector<delphi_real> >("atmeps")),
    readGlobalVector1D<delphi_real>("atmeps", 0, sizeofGlobal1D<delphi_real>("atmeps"), getKey_Ref< vector<delphi_real> >("atmeps"));

    //ieBuffz(getKey_constRef< SExtrema<int> >("buffz")),
    getKey_Ref<SExtrema<delphi_integer>>("buffz") = readGlobalVar<SExtrema<delphi_integer>>("buffz");

    //iObjectNum(getKey_constRef<int>("nobject")),
    getKey_Ref<delphi_integer>("nobject") = readGlobalVar<delphi_integer>("nobject");

    //iAtomNum(getKey_constRef<int>("natom")),
    getKey_Ref<delphi_integer>("natom") = readGlobalVar<delphi_integer>("natom");

    //iDielecBndyOdd(getKey_constRef<int>("icount2b")),
    getKey_Ref<delphi_integer>("icount2b") = readGlobalVar<delphi_integer>("icount2b");

    //fEPKT(getKey_constRef<delphi_real>("epkt")),
    getKey_Ref<delphi_real>("epkt") = readGlobalVar<delphi_real>("epkt");

    //prgdgvCrgBndyGrid(getKey_Ref< vector<SDoubleGridValue> >("cgbp")),
    readGlobalVector1D<SDoubleGridValue>("cgbp", 0, sizeofGlobal1D<SDoubleGridValue>("cgbp"), getKey_Ref< vector<SDoubleGridValue> >("cgbp"));

    //prgapAtomPdb(getKey_Ref< vector<CAtomPdb> >("delphipdb")),
    readGlobalVector1D<CAtomPdb>("delphipdb", 0, sizeofGlobal1D<CAtomPdb>("delphipdb"), getKey_Ref< vector<CAtomPdb> >("delphipdb"));

    //prgfgSurfCrgA(getKey_constRef< vector< SGrid<delphi_real> > >("scspos")),
    readGlobalVector1D<SGrid<delphi_real>>("scspos", 0, sizeofGlobal1D<SGrid<delphi_real>>("scspos"), getKey_Ref< vector< SGrid<delphi_real> > >("scspos"));

    //prgiCrgAt(getKey_constRef< vector<delphi_integer> >("crgatn")),
    readGlobalVector1D<delphi_integer>("crgatn", 0, sizeofGlobal1D<delphi_integer>("crgatn"), getKey_Ref< vector<delphi_integer> >("crgatn"));

    //atsurf(getKey_constRef< vector<delphi_integer> >("atsurf")),
    readGlobalVector1D<delphi_integer>("atsurf", 0, sizeofGlobal1D<delphi_integer>("atsurf"), getKey_Ref< vector<delphi_integer> >("atsurf"));

    //fTaylorCoeff1(getKey_constRef<delphi_real>("chi1")),
    getKey_Ref<delphi_real>("chi1") = readGlobalVar<delphi_real>("chi1");

    //fTaylorCoeff2(getKey_constRef<delphi_real>("chi2")),
    getKey_Ref<delphi_real>("chi2") = readGlobalVar<delphi_real>("chi2");

    //fTaylorCoeff3(getKey_constRef<delphi_real>("chi3")),
    getKey_Ref<delphi_real>("chi3") = readGlobalVar<delphi_real>("chi3");

    //fTaylorCoeff4(getKey_constRef<delphi_real>("chi4")),
    getKey_Ref<delphi_real>("chi4") = readGlobalVar<delphi_real>("chi4");

    //fTaylorCoeff5(getKey_constRef<delphi_real>("chi5")),
    getKey_Ref<delphi_real>("chi5") = readGlobalVar<delphi_real>("chi5");

    //schrg(getKey_Ref< vector<delphi_real> >("schrg")),
    readGlobalVector1D<delphi_real>("schrg", 0, sizeofGlobal1D<delphi_real>("schrg"), getKey_Ref< vector<delphi_real> >("schrg"));

    //ergg(getKey_Ref<delphi_real>("ergg")),
    getKey_Ref<delphi_real>("ergg") = readGlobalVar<delphi_real>("ergg");

    //ergc(getKey_Ref<delphi_real>("ergc")),
    getKey_Ref<delphi_real>("ergc") = readGlobalVar<delphi_real>("ergc");

    //ergs(getKey_Ref<delphi_real>("ergs")),
    getKey_Ref<delphi_real>("ergs") = readGlobalVar<delphi_real>("ergs");

    //ergr(getKey_Ref<delphi_real>("ergr")),
    getKey_Ref<delphi_real>("ergr") = readGlobalVar<delphi_real>("ergr");

    //ergions(getKey_Ref<delphi_real>("ergions")),
    getKey_Ref<delphi_real>("ergions") = readGlobalVar<delphi_real>("ergions");

    //ergsgaussian(getKey_Ref<delphi_real>("ergsgaussian"))
    getKey_Ref<delphi_real>("ergsgaussian") = readGlobalVar<delphi_real>("ergsgaussian");
}

void CDelphiData::sync_pre_site()
{
	readGlobalVector1D<delphi_real>("phimap", 0,sizeofGlobal1D<delphi_real>("phimap"),getKey_Ref< vector<delphi_real> >("phimap"));
	readGlobalVector1D<delphi_real>("radprb", 0, sizeofGlobal1D<delphi_real>("radprb"), getKey_Ref< vector< delphi_real > >("radprb"));
readGlobalVector1D<SGrid<delphi_real>>("xn1", 0, sizeofGlobal1D<SGrid<delphi_real>>("xn1"), getKey_Ref< vector< SGrid<delphi_real> > >("xn1"));
readGlobalVector1D<SGrid<delphi_real>>("xn2", 0, sizeofGlobal1D<SGrid<delphi_real>>("xn2"), getKey_Ref< vector< SGrid<delphi_real> > >("xn2"));
readGlobalVector1D<SGridValue<delphi_real>>("atmcrg", 0, sizeofGlobal1D<SGridValue<delphi_real>>("atmcrg"), getKey_Ref< vector< SGridValue<delphi_real> > >("atmcrg"));
readGlobalVector1D<SGrid<delphi_real>>("chgpos", 0, sizeofGlobal1D<SGrid<delphi_real>>("chgpos"), getKey_Ref< vector< SGrid<delphi_real> > >("chgpos"));
readGlobalVector1D<SGrid<delphi_real>>("scspos", 0, sizeofGlobal1D<SGrid<delphi_real>>("scspos"), getKey_Ref< vector< SGrid<delphi_real> > >("scspos"));
readGlobalVector1D<delphi_integer>("crgatn", 0, sizeofGlobal1D<delphi_integer>("crgatn"), getKey_Ref< vector< delphi_integer > >("crgatn"));
readGlobalVector1D<delphi_integer>("atsurf", 0, sizeofGlobal1D<delphi_integer>("atsurf"), getKey_Ref< vector< delphi_integer > >("atsurf"));
readGlobalVector1D<SGrid<delphi_real>>("scsnor", 0, sizeofGlobal1D<SGrid<delphi_real>>("scsnor"), getKey_Ref< vector< SGrid<delphi_real> > >("scsnor"));
readGlobalVector1D<char>("idebmap", 0, sizeofGlobal1D<char>("idebmap"), getKey_Ref< vector< char > >("idebmap"));
readGlobalVector1D<SGrid<delphi_integer>>("ibgrd", 0, sizeofGlobal1D<SGrid<delphi_integer>>("ibgrd"), getKey_Ref< vector< SGrid<delphi_integer> > >("ibgrd"));
readGlobalVector1D<delphi_real>("atmeps", 0, sizeofGlobal1D<delphi_real>("atmeps"), getKey_Ref< vector< delphi_real > >("atmeps"));
readGlobalVector1D<delphi_integer>("atndx", 0, sizeofGlobal1D<delphi_integer>("atndx"), getKey_Ref< vector< delphi_integer > >("atndx"));
readGlobalVector1D<delphi_real>("schrg", 0, sizeofGlobal1D<delphi_real>("schrg"), getKey_Ref< vector< delphi_real > >("schrg"));
readGlobalVector1D<delphi_real>("medeps", 0, sizeofGlobal1D<delphi_real>("medeps"), getKey_Ref< vector< delphi_real > >("medeps"));
}
