/**
 * @file app_delphi.cpp
 * @brief Main function to generate the executable delphicpp
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * delphicpp is the object-oriented C++ version of the program DelPhi, which takes as input a Brookhaven
 * database coordinate file format of a molecule or equivalent data for geometrical objects and/or charge
 * distributions and calculates the electrostatic potential in and around the system, using a finite
 * difference solution to the Poisson-Boltzmann equation. This can be done for any given concentration of
 * any two different salts. The system can be composed of different parts having different dielectrics.
 *
 * =================================================================================================
 *
 * [QUICK REFERENCE TO ACCESS SVN REPOSITORY ON SERVER consus.clemson.edu] \n
 * - TO SHOW LOG HISTORY FOR THIS PROJECT:
 *   	svn log svn+ssh://username\@consus.clemson.edu/home/common/delphi_distr/delphicpp
 * - TO CHECKOUT THE LATEST REVERSION:
 *   	svn co svn+ssh://username\@consus.clemson.edu/home/common/delphi_distr/delphicpp
 * - TO REMOVE/ADD FILES FROM A CHECKED-OUT PROJECT:
 *   	svn rm  <file/directory name>
 *   	svn add <file/directory name>
 * - TO REVIEW LOCAL CHANGES:
 *   	svn status
 * - TO COMMIT CHANGES:
 *   	export SVN_EDITOR=/usr/bin/vim
 *   	svn commit
 * - TO UPDATE YOUR LOCAL FILES:
 *   	svn update
 *
 *   (See "subversion user guide" created by S. Sarkar and C. Li for other svn options)
 *
 * =================================================================================================
 *
 **/

#include "../delphi/delphi_data.h"
#include "../delphi/delphi_datamarshal.h"
#include "../space/space.h"
#include "../solver/solver_fastSOR.h"
#include "../energy/energy.h"
#include "../site/site.h"

#include <time.h>
#include <math.h>
#include <iostream>
#include <cstddef>
#include <sstream>
#include <string>

using namespace std;

#define DEVELOPER
#define VERBOSE

int main(int argc, char *argv[])
{
    /*
     * bool values are inserted/extracted by their textual representation: either true or false, instead of integral values.
     * This flag can be unset with the noboolalpha manipulator.
     */


    size_t MAXWIDTH = 45;
    string info_string = " Info>";

    #ifdef VERBOSE
    cout << boolalpha;
    cerr << boolalpha;
    #endif

    #ifdef DEVELOPER
    cout << fixed << setprecision(7); //cout.precision(15)
    #else
    cout << fixed << setprecision(3); //cout.precision(7)
    #endif

    try
    {
        CTimer * pTester = new CTimer; // record execution time

        pTester->start();

        shared_ptr<CTimer> pTimer(new CTimer); // record execution time

        //********************************************************************************//
        //                                                                                //
        //                  read the parameter, size, charge, PDB files                   //
        //                                                                                //
        //********************************************************************************//

        //cout << "DPLT IDELPHI Creating" << endl;
        //---------- a shared_ptr to an object of IDataContainer
        shared_ptr<IDataContainer> pDataContainer(new CDelphiData(argc, argv, pTimer));
        //cout << "DPLT IDELPHI Created" << endl;

        #ifdef DEBUG_DELPHI_SPACE
        cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SPACE] data container is written into file test_delphicpp_atbeginning.dat---\n\n"
        << "\033[0m";

        pDataContainer->showMap("test_delphicpp_atbeginning.dat");
        #endif

        //pDataContainer->showMap("test_delphicpp_before_space.dat");

        #ifdef DEVELOPER
        cout << "\n\n---------- delphicpp finishes IO on ";
        pTester->showTime();
        cout << "\n\n";
        #endif

        //********************************************************************************//
        //                                                                                //
        //    realize an object of CDelphiSpace class to construct molecular surfaces     //
        //                                                                                //
        //********************************************************************************//
        int& inhomo(pDataContainer->getKey_Ref<int>("inhomo"));
        const int& iGaussian(pDataContainer->getKey_constRef<int>("gaussian"));
        bool& logs(pDataContainer->getKey_Ref<bool>("logs"));
        bool& logg(pDataContainer->getKey_Ref<bool>("logg"));

        //iGaussian=pDataContainer->getKey_constRef<int>("gaussian");
        //inhomo=pDataContainer->getKey_Ref<int>("inhomo");
        //logs=pDataContainer->getKey_Ref<bool>("logs");
        //logg=pDataContainer->getKey_Ref<bool>("logg");
        inhomo = 0;

        if (iGaussian == 1 && logs)
        {
            logg = true; //for gaussian
            inhomo = 1;
        }
		
        unique_ptr<IAbstractModule> pSpace(new CDelphiSpace(pDataContainer, pTimer));

        pSpace->run();

        pSpace.reset();

        #ifdef DEBUG_DELPHI_SPACE
        cout << "---[DEBUG_DELPHI_SPACE] data container is written into file test_delphicpp_aftersurf.dat---\n\n";

        pDataContainer->showMap("test_delphicpp_aftersurf.dat");
        #endif

        //pDataContainer->showMap("test_delphicpp_aftersurf.dat");

        #ifdef DEVELOPER
        cout << "\n\n---------- delphicpp finishes SPACE class on ";
        pTester->showTime();
        cout << "\n\n";

        cout << endl;
        #endif

        #ifdef VERBOSE
        if ( !(iGaussian == 1 && inhomo == 0 && logs) )
            cout << left << setw(MAXWIDTH) << " Number of atom coordinates read" << " : "
                 << pDataContainer->getKey_constRef< delphi_integer > ("natom") << endl;
        #endif

        if (pDataContainer->getKey_constRef<bool>("isolv"))
        {
            if ( !(iGaussian == 1 && inhomo == 0 && logs) )
            {
                cout << left << setw(MAXWIDTH) << " Total number of assigned charges" << " : "
                     << pDataContainer->getKey_constRef< delphi_integer > ("nqass") << endl;
                cout << left << setw(MAXWIDTH) << " Net assigned charge" << " : "
                     << right << setw(10) << pDataContainer->getKey_constRef< delphi_real > ("qnet") << endl;
                cout << left << setw(MAXWIDTH) << " Assigned positive charge" << " : "
                     << right << setw(10) << pDataContainer->getKey_constRef< delphi_real > ("qplus") << endl;
                cout << left << setw(MAXWIDTH) << " Centred at (gu)" << " : "
                     << pDataContainer->getKey_constRef < SGrid<delphi_real> > ("cqplus").nX << " " << right << setw(10)
                     << pDataContainer->getKey_constRef< SGrid<delphi_real> > ("cqplus").nY << " " << right << setw(10)
                     << pDataContainer->getKey_constRef< SGrid<delphi_real> > ("cqplus").nZ << endl;
                cout << left << setw(MAXWIDTH) << " Assigned negative charge"
                     << " : " << pDataContainer->getKey_constRef< delphi_real > ("qmin") << endl;
                cout << left << setw(MAXWIDTH) << " Centred at (gu)" << " : "
                     << pDataContainer->getKey_constRef < SGrid<delphi_real> > ("cqmin").nX << " " << right << setw(10)
                     << pDataContainer->getKey_constRef < SGrid<delphi_real> > ("cqmin").nY << " " << right << setw(10)
                     << pDataContainer->getKey_constRef < SGrid<delphi_real> > ("cqmin").nZ << endl;

                #ifdef VERBOSE
                cout << left << setw(MAXWIDTH) << " Number of dielectric boundary points" << " : "
                     << pDataContainer->getKey_constRef < delphi_integer > ("ibnum") << endl;
                #endif

                cout << endl;
            }

            if (pDataContainer->getKey_constRef<bool>("iexun") && 0 == pDataContainer->getKey_constRef < delphi_integer > ("ibnum"))
                throw CNoBndyAndDielec(pTimer);

            //********************************************************************************//
            //                                                                                //
            //   realize an object of CDelphiFastSOR class to calculate potentials on grids   //
            //                                                                                //
            //********************************************************************************//

            #ifdef DEBUG_DELPHI_SOLVER
            cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SOLVER] data container is read from file test_fromsurf.dat---\n\n" << "\033[0m";

            pDataContainer->reset("test_fromsurf.dat");

            cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SOLVER] data container is written into file test_delphicpp_beforeitr.dat---\n\n"
            << "\033[0m";

            pDataContainer->showMap("test_delphicpp_beforeitr.dat");

            #endif // DEBUG_DELPHI_SOLVER
        }

        unique_ptr<IAbstractModule> pSolver(new CDelphiFastSOR(pDataContainer, pTimer));

        pSolver->run();

        pSolver.reset();

        #ifdef DEBUG_DELPHI_SOLVER
        cout << "---[DEBUG_DELPHI_SOLVER] data container is written into file test_delphicpp_afteritr.dat---\n\n";

        pDataContainer->showMap("test_delphicpp_afteritr.dat");
        #endif // DEBUG_DELPHI_SOLVER

        //pDataContainer->showMap("test_delphicpp_afteritr.dat");

        #ifdef DEVELOPER
        cout << "\n\n---------- delphicpp finishes SOLVER class on ";
        pTester->showTime();
        cout << "\n\n";
        #endif

        //********************************************************************************//
        //                                                                                //
        //          realize an object of CDelphiEnergy class to calculate energies        //
        //                                                                                //
        //********************************************************************************//

        #ifdef DEBUG_DELPHI_ENERGY
        cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_ENERGY] data container is read from file test_fromsurf.dat "
        << "and test_fromsolver.dat---\n\n" << "\033[0m";

        pDataContainer->reset("test_fromsurf.dat");

        pDataContainer->reset("test_fromsolver.dat");

        cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SOLVER] data container is written into file test_delphicpp_beforeenergy.dat---\n\n"
        << "\033[0m";

        pDataContainer->showMap("test_delphicpp_beforeenergy.dat");

        #endif // DEBUG_DELPHI_ENERGY


        unique_ptr<IAbstractModule> pEnergy(new CDelphiEnergy(pDataContainer, pTimer));

        pEnergy->run();

        pEnergy.reset();

        #ifdef DEBUG_DELPHI_ENERGY
        cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_ENERGY] data container is written into file test_delphicpp_aftereng.dat---\n\n"
        << "\033[0m";

        pDataContainer->showMap("test_delphicpp_aftereng.dat");
        #endif // DEBUG_DELPHI_ENERGY

        //pDataContainer->showMap("test_delphicpp_aftereng.dat");

        #ifdef DEVELOPER
        cout << "\n\n---------- delphicpp finishes ENERGY class on ";
        pTester->showTime();
        cout << "\n\n";
        #endif

        if (iGaussian == 1 && inhomo == 1 && logs) //second run for Gaussian
        {

            //pDataContainer.reset();

            //shared_ptr<IDataContainer> pDataContainer( new CDelphiData(argc,argv,pTimer) );
            inhomo = 0;

            unique_ptr < IAbstractModule > pSpace(new CDelphiSpace(pDataContainer, pTimer));
            pSpace->run();
            pSpace.reset();

            #ifdef DEBUG_DELPHI_SPACE
            cout << "---[DEBUG_DELPHI_SPACE] data container is written into file test_delphicpp_aftersurf.dat---\n\n";
            pDataContainer->showMap("test_delphicpp_aftersurf.dat");
            #endif

            unique_ptr < IAbstractModule > pSolver(new CDelphiFastSOR(pDataContainer, pTimer));
            pSolver->run();
            pSolver.reset();

            #ifdef DEBUG_DELPHI_SOLVER
            cout << "---[DEBUG_DELPHI_SOLVER] data container is written into file test_delphicpp_afteritr.dat---\n\n";

            pDataContainer->showMap("test_delphicpp_afteritr.dat");
            #endif
            unique_ptr < IAbstractModule > pEnergy(new CDelphiEnergy(pDataContainer, pTimer));
            pEnergy->run();
            pEnergy.reset();
        }

        //********************************************************************************//
        //                                                                                //
        //               realize an object of CSite class to write site info              //
        //                                                                                //
        //********************************************************************************//

        #ifdef DEBUG_DELPHI_SITE
        cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SITE] data container is read from file test_fromsurf.dat, "
        << "test_fromsolver.dat and test_fromenergy.dat---\n\n" << "\033[0m";

        pDataContainer->reset("test_fromsurf.dat");
        pDataContainer->reset("test_fromsolver.dat");
        pDataContainer->reset("test_fromenergy.dat");

        cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SOLVER] data container is written into file test_delphicpp_beforesite.dat---\n\n"
        << "\033[0m";

        pDataContainer->showMap("test_delphicpp_beforesite.dat");
        #endif // DEBUG_DELPHI_SITE

        unique_ptr<CSite> pSite(new CSite(pDataContainer, pTimer));
        
	if (pDataContainer->getKey_constRef<bool>("isite"))
        {
            int iisitsf = 0;
            if (pDataContainer->getKey_Ref<bool>("isitsf")) iisitsf = 1;
            pSite->writeSite(iisitsf);
        }

        if (pDataContainer->getKey_constRef<bool>("phiwrt")) pSite->writePhi();
        if (pDataContainer->getKey_constRef < delphi_real > ("radipz") > 0.) pSite->wirtePAnalysis();
        
	pSite.reset();
        
#ifdef DEBUG_DELPHI_SITE
        cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SITE] data container is written into file test_delphicpp_atend.dat---\n\n"
        << "\033[0m";

        pDataContainer->showMap("test_delphicpp_atend.dat");
        #endif // DEBUG_DELPHI_SITE

        #ifdef DEVELOPER
        cout << "\n\n---------- delphicpp finishes SITE class on ";
        pTester->showTime();
        cout << "\n\n";
        #endif

        pDataContainer.reset();
        
	pTimer->exit();
        pTimer.reset();

        delete pTester;

    } // ---------- end of try block
    catch (CException&)
    {
        cerr << "\n\n ......... PROGRAM ABORTS WITH AN EXCEPTION AND "
                << CWarning::iWarningNum << " WARNING(S) ........\n\n";

        return 1;
    }

    cout << "\n\n .........  PROGRAM EXITS SUCCESSFULLY : WITH TOTAL "
            << CWarning::iWarningNum << " WARNING(S) ........\n\n";

	int numWarnings = 0;
	if (iMaxWarn >= CWarning::iWarningNum) 
	{
		numWarnings = CWarning::iWarningNum;
	} 
	else 
	{
		numWarnings = iMaxWarn;
	}

	cout << " .....................  SHOWING " << numWarnings << " OUT OF " 
	     <<  CWarning::iWarningNum << " WARNING(S) ..................." << endl;


	// ARGO: modified the method of printing warnings based on the number requested by the user
	// Specifically, the buffer is made to have "|" char which act as delimiters.
	// The buffer is split by that delimiter.

	string warningStr = cwarn.rdbuf()->str()  + "|";
	int warnCount = -1;
	std::size_t current, previous = 0;

	current = warningStr.find_first_of("|");
	while (current != std::string::npos ) 
	{
		if (warnCount == iMaxWarn) break; 
		cout << warningStr.substr(previous, current - previous) << endl;
		previous = current + 1;
		current = warningStr.find_first_of("|", previous);
		warnCount++;
	}

	cerr << "\033[0m";
    cout.unsetf(ios_base::floatfield); // return to cout default notation

    return 0;

}
