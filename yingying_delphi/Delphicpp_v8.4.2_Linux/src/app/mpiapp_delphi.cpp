/**
 * @file mpiapp_delphi.cpp
 * @brief Main function to generate the executable mpi version of delphicpp
 *
 * @author Chuan Li, cli@wcupa.edu
 *
 * delphicpp is the object-oriented C++ version of the program DelPhi, which takes as input a Brookhaven
 * database coordinate file format of a molecule or equivalent data for geometrical objects and/or charge
 * distributions and calculates the electrostatic potential in and around the system, using a finite
 * difference solution to the Poisson-Boltzmann equation. This can be done for any given concentration of
 * any two different salts. The system can be composed of different parts having different dielectrics.
 *
 * [VERSION HISTORY] \n
 * - r01   chuan   2017Feb23   The 1st version of mpi implementation of delphi program based on delphicpp v76.
 */

#include "../delphi/delphi_data.h"
#include "../space/space.h"
#include "../solver/solver_fastSOR.h"
#include "../energy/energy.h"
#include "../site/site.h"
#include "../dplt/dplt.h"
#include "../dplt/mem_usage.h"

#include <mpi.h>

#include <sstream>
#include <iostream>
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstddef>
#include <sstream>

using namespace std;
//#define VERBOSE
//#define DEVELOPER
int run = 1;
int space_pause = 1;

int main(int argc, char *argv[]) 
{
    /*
     * bool values are inserted/extracted by their textual representation: either true or false, 
     * instead of integral values.
     * This flag can be unset with the noboolalpha manipulator.
     */
    std::ostringstream local;
    std::streambuf *cerr_buff = std::cerr.rdbuf(); // save pointer to std::cout buff

    std::cerr.rdbuf(local.rdbuf()); // substitute internal std::cout buffer with

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

        shared_ptr < CTimer > pTimer(new CTimer); // record execution time

        if(argc>2)run = atoi(argv[2]);
        while (run == 0)
        {
            sleep(5);
        }
        
        //********************************************************************************//
        //                                                                                //
        //                  read the parameter, size, charge, PDB files                   //
        //                                                                                //
        //********************************************************************************//

        /*
         * ++++++++++ MPI is initialzied when a new CDelphiData object is realized ++++++++++
         */
        //cout << "DPLT IDELPHI Creating" << endl;
        //---------- a shared_ptr to an object of IDataContainer
        shared_ptr < IDataContainer > pDataContainer(new CDelphiData(argc, argv, pTimer));
        //cout << "DPLT IDELPHI Created" << endl;
        
        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp finishes IO on ";
            pTester->showTime();
            cout << "\n\n";     
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                    cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
        }
        #endif
        
        int& inhomo(pDataContainer->getKey_Ref<int>("inhomo"));
        const int& iGaussian(pDataContainer->getKey_constRef<int>("gaussian"));
        bool& logs(pDataContainer->getKey_Ref<bool>("logs"));
        bool& logg(pDataContainer->getKey_Ref<bool>("logg"));

        inhomo = 0;

        if (iGaussian == 1 && logs)
        {
            logg = true; //for gaussian
            inhomo = 1;
        }

        /*
         * -------------------- an object to Delphi Space class --------------------
         */
        pDataContainer->sync_pre_space();
        ddm::barrier();
        
        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp Prepared data for SPACE on ";
            pTester->showTime();
            cout << "\n\n";
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
        }
        #endif

        unique_ptr<IAbstractModule> pSpace(new CDelphiSpace(pDataContainer, pTimer));

        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp finishes Created SPACE object on ";
            pTester->showTime();
            cout << "\n\n";
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
         }
        #endif

        //if (pDataContainer->myid() == 0) pDataContainer->showMap("rank0_mpiapp_beforeSpace.dat");
        //if (pDataContainer->myid() == 1) pDataContainer->showMap("rank1_mpiapp_beforeSpace.dat");

        pSpace->run();
        ddm::barrier();

        //if (pDataContainer->myid() == 0) pDataContainer->showMap("rank0_mpiapp_afterSpace.dat");
        //if (pDataContainer->myid() == 1) pDataContainer->showMap("rank1_mpiapp_afterSpace.dat");

        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp finishes SPACE Run on ";
            pTester->showTime();
            cout << "\n\n";
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
        }
        #endif
 
        //pSpace->reallocate_global(0);
        //ddm::barrier();
        pSpace->write_to_global();
        ddm::barrier();

        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp writes SPACE DATA to Global on ";
            pTester->showTime();
            cout << "\n\n";
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
        }
        #endif

        if(argc>3)space_pause = atoi(argv[3]);
        while (0 == space_pause)
        {
            sleep(5);
        }

        pSpace.reset();

        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp Reset SPACE Object on ";
            pTester->showTime();
            cout << "\n\n";
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
        }
        #endif

        //if (pDataContainer->myid() == 0) pDataContainer->showMap("rank0_mpiapp_afterSpaceDetroied.dat");
        //if (pDataContainer->myid() == 1) pDataContainer->showMap("rank1_mpiapp_afterSpaceDetroied.dat");

        /*
         * -------------------- sync local and global data containers before iteration --------------------
         */
        pDataContainer->sync_pre_solver();
        ddm::barrier();

        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp Prepared data for SOLVER on ";
            pTester->showTime();
            cout << "\n\n";
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
        }
        #endif

        if (pDataContainer->myid() == 0)
        {
            #ifdef VERBOSE
            if ( !(iGaussian == 1 && inhomo == 0 && logs) )
                cout << left << setw(MAXWIDTH) << " Number of atom coordinates read" << right << setw(10) 
                     << pDataContainer->getKey_constRef< delphi_integer > ("natom") << endl;
            #endif

            if (pDataContainer->getKey_constRef<bool>("isolv")) 
            {
                if ( !(iGaussian == 1 && inhomo == 0 && logs) )
                {
                    cout << left << setw(MAXWIDTH) << " Total number of assigned charges" << " : " 
                            << pDataContainer->getKey_constRef < delphi_integer > ("nqass") << endl;
                    cout << left << setw(MAXWIDTH) << " Net assigned charge" << " : " 
                            << pDataContainer->getKey_constRef < delphi_real > ("qnet") << endl;
                    cout << left << setw(MAXWIDTH) << " Assigned positive charge" << " : " 
                            << pDataContainer->getKey_constRef < delphi_real > ("qplus") << endl;
                    cout << left << setw(MAXWIDTH) << " Centred at (gu)" << " : " 
                         << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqplus").nX << " "
                         << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqplus").nY << " "
                         << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqplus").nZ << endl;
                    cout << left << setw(MAXWIDTH) << " Assigned negative charge" << " : " 
                         << pDataContainer->getKey_constRef < delphi_real > ("qmin") << endl;
                    cout << left << setw(MAXWIDTH) << " Centred at (gu)" << " : " 
                         << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqmin").nX << " "
                         << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqmin").nY << " "
                         << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqmin").nZ << endl;
                    
                    #ifdef VERBOSE
                    cout << left << setw(MAXWIDTH) << " Number of dielectric boundary points" << " : " 
                            << pDataContainer->getKey_constRef < delphi_integer > ("ibnum") << endl;
                    #endif
                    cout << endl;
                }
            
                if (pDataContainer->getKey_constRef<bool>("iexun") && 0 == pDataContainer->getKey_constRef < delphi_integer > ("ibnum")) 
                    throw CNoBndyAndDielec(pTimer);
            }   
        }
        
        /*
         * -------------------- an object to Delphi fastSolver class --------------------
         */        
        unique_ptr < IAbstractModule > pSolver(new CDelphiFastSOR(pDataContainer, pTimer));         
        
        //if (pDataContainer->myid() == 0) pDataContainer->showMap("rank0_mpiapp_beforeSolver.dat");
        //if (pDataContainer->myid() == 1) pDataContainer->showMap("rank1_mpiapp_beforeSolver.dat");
        //if (pDataContainer->myid() == 2) pDataContainer->showMap("rank2_mpiapp_beforeSolver.dat");
        //if (pDataContainer->myid() == 3) pDataContainer->showMap("rank3_mpiapp_beforeSolver.dat");
        
        pSolver->mpi_run();
        ddm::barrier();
        
        //if (pDataContainer->myid() == 0) pDataContainer->showMap("rank0_mpiapp_afterSolver.dat");
        //if (pDataContainer->myid() == 1) pDataContainer->showMap("rank1_mpiapp_afterSolver.dat");
        //if (pDataContainer->myid() == 2) pDataContainer->showMap("rank2_mpiapp_afterSolver.dat");
        //if (pDataContainer->myid() == 3) pDataContainer->showMap("rank3_mpiapp_afterSolver.dat");

        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp finishes SOLVER calculation on ";
            pTester->showTime();
            cout << "\n\n";
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
        }
        #endif
 
        pSolver->reallocate_global(0);
        ddm::barrier();
        if (pDataContainer->myid() == 0)
        {
            pSolver->write_to_global();
        }

        ddm::barrier();
        pSolver.reset();

        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp finishes SOLVER DATA Transfer on ";
            pTester->showTime();
            cout << "\n\n";
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
        }
        #endif
 
        //if (pDataContainer->myid() == 0) pDataContainer->showMap("rank0_mpiapp_afterSolverSync.dat");
        //if (pDataContainer->myid() == 1) pDataContainer->showMap("rank1_mpiapp_afterSolverSync.dat");
        //if (pDataContainer->myid() == 2) pDataContainer->showMap("rank2_mpiapp_afterSolver.dat");
        //if (pDataContainer->myid() == 3) pDataContainer->showMap("rank3_mpiapp_afterSolver.dat");      
        //cout << "rank " << pDataContainer->myid() << ": after solver." << endl;  
        
        /*
         * --------------- sych local and global data containers before calculating energies ---------------
         */
        pDataContainer->sync_pre_energy();
        ddm::barrier();

        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp finishes Energy init on ";
            pTester->showTime();
            cout << "\n\n";
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
        }
        #endif

        //if (pDataContainer->myid() == 0) pDataContainer->showMap("rank0_mpiapp_beforeEnergy.dat");
        //if (pDataContainer->myid() == 1) pDataContainer->showMap("rank1_mpiapp_beforeEnergy.dat");
        //if (pDataContainer->myid() == 2) pDataContainer->showMap("rank2_mpiapp_beforeEnergy.dat");
        //if (pDataContainer->myid() == 3) pDataContainer->showMap("rank3_mpiapp_beforeEnergy.dat");
        //cout << "rank " << pDataContainer->myid() << ": after pre_energy sync." << endl; 
        
        /*
         * -------------------- an object to Delphi energy class --------------------
         */
        
        unique_ptr<IAbstractModule> pEnergy(new CDelphiEnergy(pDataContainer, pTimer));  
        
        pEnergy->mpi_run();
        ddm::barrier();

        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp finishes ENERGY run on ";
            pTester->showTime();
            cout << "\n\n";
        }
        #endif

        //if (pDataContainer->myid() == 0) pDataContainer->showMap("rank0_mpiapp_afterEnergy.dat");
        //if (pDataContainer->myid() == 1) pDataContainer->showMap("rank1_mpiapp_afterEnergy.dat");
        //if (pDataContainer->myid() == 2) pDataContainer->showMap("rank2_mpiapp_afterEnergy.dat");
        //if (pDataContainer->myid() == 3) pDataContainer->showMap("rank3_mpiapp_afterEnergy.dat");
        //cout << "rank " << pDataContainer->myid() << ": after energy mpi_run" << endl;          
       
        pEnergy->reallocate_global(0);
        ddm::barrier();

        if (pDataContainer->myid() == 0)
        {
            pEnergy->write_to_global();
        }
        ddm::barrier();

        #ifdef DEVELOPER
        if (pDataContainer->myid() == 0)
        {
            cout << "\n\n---------- delphicpp finishes ENERGY on ";
            pTester->showTime();
            cout << "\n\n";
        }

        for(int myproc=0; myproc<pDataContainer->num_procs(); myproc++)
        {
            if(myproc == pDataContainer->myid())
            {
                cout << "On Node " << myproc << " Memory usage: " << getCurrentRSS() << "\n";
            }
            ddm::barrier();
        }
        #endif

        pEnergy.reset();
        
        //cout << "rank " << pDataContainer->myid() << ": after energy" << endl; 

        /*
         * -------------------- below is run on master --------------------
         */
        if (pDataContainer->myid() == 0)
        {
            if (iGaussian == 1 && inhomo == 1 && logs) //second run for Gaussian
            {
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
            
            #ifdef DEBUG_DELPHI_SITE
            cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SITE] data container is read from file test_fromsurf.dat, "
                 << "test_fromsolver.dat and test_fromenergy.dat---\n\n" << "\033[0m";
            pDataContainer->reset("test_fromsurf.dat");
            pDataContainer->reset("test_fromsolver.dat");
            pDataContainer->reset("test_fromenergy.dat");
            cerr << "\n\033[1;35m" 
                 << "---[DEBUG_DELPHI_SOLVER] data container is written into file test_delphicpp_beforesite.dat---\n\n"<< "\033[0m";
           // pDataContainer->showMap("test_delphicpp_beforesite.dat");
            #endif // DEBUG_DELPHI_SITE
            
            #ifdef DEVELOPER
            cout << "\n\n---------- delphicpp starting SITE class on ";
            pTester->showTime();
            cout << "\n\n";
            #endif
	 //pDataContainer->readGlobalVector1D<delphi_real>("phimap", 0, pDataContainer->sizeofGlobal1D<delphi_real>("phimap"), pDataContainer->getKey_Ref< vector<delphi_real> >("phimap"));
            pDataContainer->sync_pre_site();
            unique_ptr < CSite > pSite(new CSite(pDataContainer, pTimer));
            #ifdef DEVELOPER
            cout << "\n\n---------- delphicpp SITE class started ";
            cout << "\n\n";
            #endif
            if (pDataContainer->getKey_constRef<bool>("isite")) 
            {
                int iisitsf = 0;
                if (pDataContainer->getKey_Ref<bool>("isitsf")) iisitsf = 1;
                pSite->writeSite(iisitsf);
            }
            if (pDataContainer->getKey_constRef<bool>("phiwrt")) 
                pSite->writePhi();
            if (pDataContainer->getKey_constRef < delphi_real > ("radipz") > 0.)
                pSite->wirtePAnalysis();
            pSite.reset();
            #ifdef DEBUG_DELPHI_SITE
            cerr << "\n\033[1;35m" 
                 << "---[DEBUG_DELPHI_SITE] data container is written into file test_delphicpp_atend.dat---\n\n" << "\033[0m";
            //pDataContainer->showMap("test_delphicpp_atend.dat");
            #endif // DEBUG_DELPHI_SITE
            
            #ifdef DEVELOPER
            cout << "\n\n---------- delphicpp finishes SITE class on ";
            pTester->showTime();
            cout << "\n\n";
            #endif          
        }
	if (pDataContainer->myid() == 0) pTimer->exit();
        pTimer.reset();
        pDataContainer.reset();
        delete pTester;
    } // ---------- end of try block
    catch (CException&) 
    {
        cerr << "\n\n ......... PROGRAM ABORTS WITH AN EXCEPTION AND " << CWarning::iWarningNum << " WARNING(S) ........\n\n";
        return 1;
    }
    //cout << "\n\n .........  PROGRAM EXITS SUCCESSFULLY : WITH TOTAL " << CWarning::iWarningNum << " WARNING(S) ........\n\n";
    std::cerr.rdbuf(cerr_buff);
    // print 'local' content
    std::cerr << local.str() << "\n";
    cout.unsetf(ios_base::floatfield); // return to cout default notation

    return 0;
}
