#include "space.h"


SGrid<delphi_integer> CDelphiSpace::calcCPUDim()
{

    delphi_integer ncpus = 1;

#ifdef PARALLEL_MPI
    ncpus = pdc->num_procs();
#endif //PARALLEL_MPI

    delphi_integer* _cpuDim = new delphi_integer[3];

    bool change_cpu_num = false;

    //loop over 3 dimensions to find the best size
    for (int i(0); i < 3; ++i)
    {
        //for 3d, the optimal blocks on each dimension would be the cubic  root of nprocs
        //for 2d, it would be the square root
        double dim_optimal = round(pow(ncpus, 1.f / (3.f - i)));

        //the best dimension size would be a factor of nprocs that is very close to
        //the optimal dimension
        for (int factor(1); factor <= ncpus; ++factor)
        {
            if (ncpus % factor == 0) // if there is a factor
            {
                if (factor >= dim_optimal) //it is slightly larger than the optimal size
                {
                    if (factor >= dim_optimal * 2) //but not too large
                    {
                        //if it is too large, an error will be popped
                        change_cpu_num = true;
                    }
                    //this is the size of dimension i
                    _cpuDim[i] = factor;
                    break;
                }
            }
        }
        //now for the other two dimensions
        ncpus = ncpus / _cpuDim[i];
    }
    SGrid<delphi_integer> result = { _cpuDim[0], _cpuDim[1], _cpuDim[2] };
    delete _cpuDim;
    return result;
}

//---------------------------------------------------------------------

SGrid<delphi_integer> CDelphiSpace::calcMyCPULocation()
{
    SGrid<delphi_integer> result = {
            my_id / (cpuDim.nY*cpuDim.nZ) ,
            (my_id % (cpuDim.nY*cpuDim.nZ)) / cpuDim.nZ ,
            my_id % cpuDim.nZ };
    return result;
}

//---------------------------------------------------------------------

SGrid<delphi_integer> CDelphiSpace::calcGridDim()
{
    //allocate memory for dim array
    //dim array stores the size of all blocks
    localGridDim = new SGrid<delphi_integer> **[cpuDim.nX];
    for (int x = 0; x < cpuDim.nX; x++)
    {
        localGridDim[x] = new SGrid<delphi_integer>*[cpuDim.nY];
        for (int y = 0; y < cpuDim.nY; y++)
        {
            localGridDim[x][y] = new SGrid<delphi_integer>[cpuDim.nZ];
        }
    }

    //calculate the dimension of grid in each CPU block
    int x, y, z, grid_x, grid_y, grid_z;
    //for each CPU block
    for (x = 0, grid_x = global_iGrid.nX; x < cpuDim.nX; x++)
    {
        for (y = 0, grid_y = global_iGrid.nY; y < cpuDim.nY; y++)
        {
            for (z = 0, grid_z = global_iGrid.nZ; z < cpuDim.nZ; z++)
            {
                //the number of grid on each dimension is is the left over grids divided by the left over CPUs
                localGridDim[x][y][z].nX = round(grid_x / (cpuDim.nX - x));
                localGridDim[x][y][z].nY = round(grid_y / (cpuDim.nY - y));
                localGridDim[x][y][z].nZ = round(grid_z / (cpuDim.nZ - z));

                //these grids are assigned to CPUs on Z dimension
                grid_z -= localGridDim[x][y][z].nZ;
            }
            //these grids are assigned to CPUs on Y dimension
            grid_y -= localGridDim[x][y][0].nY;
        }
        //these grids are assigned to CPUs on X dimension
        grid_x -= localGridDim[x][0][0].nX;
    }

    mySizeNoBuffer = localGridDim[myCPULocation.nX][myCPULocation.nY][myCPULocation.nZ];
    return mySizeNoBuffer;
}

//---------------------------------------------------------------------

SGrid<delphi_integer> CDelphiSpace::calcLocalGridStartEnd()
{
    //allocate memory for dim array
    //dim array stores the size of all blocks
    SGrid<delphi_integer> *** localStart = new SGrid<delphi_integer> **[cpuDim.nX + 1];
    SGrid<delphi_integer> *** localEnd = new SGrid<delphi_integer> **[cpuDim.nX + 1];

    for (int x = 0; x <= cpuDim.nX; x++)
    {
        localStart[x] = new SGrid<delphi_integer>*[cpuDim.nY + 1];

        for (int y = 0; y <= cpuDim.nY; y++)
        {
            localStart[x][y] = new SGrid<delphi_integer>[cpuDim.nZ + 1];
        }
    }

    //evenly divide the atoms into sub divisions
    vector<SGrid<delphi_real>> atomCoor;

    for (std::vector <CAtomPdb>::iterator it = delphipdb.begin(); it != delphipdb.end(); it++)
    {
        atomCoor.push_back(it->getPose());
    }

    size_t xStart;
    size_t xSize;
    delphi_integer xEnd;
    delphi_integer xStartGrid;

    size_t yStart;
    size_t ySize;
    delphi_integer yEnd;
    delphi_integer yStartGrid;

    size_t zStart;
    size_t zSize;
    delphi_integer zEnd;
    delphi_integer zStartGrid;


    std::sort(atomCoor.begin(), atomCoor.end(), compareCoorX);

    vector<SGrid<delphi_real>> XCoorSorted(atomCoor);
    xSize = max(XCoorSorted.size() - 1, size_t(0));

    for (int x = 0; x <= cpuDim.nX; x++)
    {
        xStart = xSize / cpuDim.nX * x;
        xEnd = min(xSize / cpuDim.nX + xStart - 1, xSize);
        xEnd = max(delphi_integer(xStart), xEnd);

        xStartGrid = global_iGrid.nX - 1;
        if (XCoorSorted.size() > 0)
        {
            delphi_integer xStartCoor = XCoorSorted[xStart].nX;
            xStartGrid = (xStartCoor - global_cOldMid.nX)*fScale + (global_iGrid.nX + 1) / 2;
        }

        //Now sort these atoms by Y corrdinate
        vector<SGrid<delphi_real>> YCoorSorted;
        if (XCoorSorted.size() > 0)
        {
            vector<SGrid<delphi_real>> temp(XCoorSorted.begin() + xStart, XCoorSorted.begin() + xEnd);
            YCoorSorted.swap(temp);
        }
        std::sort(YCoorSorted.begin(), YCoorSorted.end(), compareCoorY);
        ySize = max(YCoorSorted.size() - 1, size_t(0));

        for (int y = 0; y <= cpuDim.nY; y++)
        {
            yStart = ySize / cpuDim.nY * y;
            yEnd = min(ySize / cpuDim.nY * +yStart - 1, ySize);
            yEnd = max(delphi_integer(yStart), yEnd);

            yStartGrid = global_iGrid.nY - 1;

            if (YCoorSorted.size() > 0)
            {
                delphi_real yStartCoor = YCoorSorted[yStart].nY;
                yStartGrid = (yStartCoor - global_cOldMid.nY)*fScale + (global_iGrid.nY + 1) / 2;
            }

            //Now sort these atoms bz Z corrdinate
            vector<SGrid<delphi_real>> ZCoorSorted;
            if (YCoorSorted.size() > 0)
            {
                vector<SGrid<delphi_real>>temp(YCoorSorted.begin() + yStart, YCoorSorted.begin() + yEnd);
                ZCoorSorted.swap(temp);
            }
            std::sort(ZCoorSorted.begin(), ZCoorSorted.end(), compareCoorZ);
            zSize = max(ZCoorSorted.size() - 1, size_t(0));

            for (int z = 0; z <= cpuDim.nZ; z++)
            {
                zStart = zSize / cpuDim.nZ * z;
                zEnd = min(zSize / cpuDim.nZ *  +zStart - 1, zSize);
                zEnd = max(delphi_integer(zStart), zEnd);

                zStartGrid = global_iGrid.nZ - 1;

                if (ZCoorSorted.size() > 0)
                {
                    delphi_real zStartCoor = ZCoorSorted[zStart].nZ;
                    zStartGrid = (zStartCoor - global_cOldMid.nZ)*fScale + (global_iGrid.nZ + 1) / 2;
                }

                //for the else, do such calculation
                localStart[x][y][z].nX = xStartGrid;
                localStart[x][y][z].nY = yStartGrid;
                localStart[x][y][z].nZ = zStartGrid;


                //for the first CPU in each dimension, the start is 0
                if (x == 0) localStart[x][y][z].nX = 0;
                if (y == 0) localStart[x][y][z].nY = 0;
                if (z == 0) localStart[x][y][z].nZ = 0;

                //for the last CPU in each dimension, the end is iGrid-1
                if (x == cpuDim.nX) localStart[x][y][z].nX = global_iGrid.nX;
                if (y == cpuDim.nY) localStart[x][y][z].nY = global_iGrid.nY;
                if (z == cpuDim.nZ) localStart[x][y][z].nZ = global_iGrid.nZ;
            }
        }
    }


    myStartNoBuffer.nX = localStart[myCPULocation.nX][myCPULocation.nY][myCPULocation.nZ].nX;
    myEndNoBuffer.nX = localStart[myCPULocation.nX + 1][myCPULocation.nY][myCPULocation.nZ].nX - 1;

    myStartNoBuffer.nY = localStart[myCPULocation.nX][myCPULocation.nY][myCPULocation.nZ].nY;
    myEndNoBuffer.nY = localStart[myCPULocation.nX][myCPULocation.nY + 1][myCPULocation.nZ].nY - 1;

    myStartNoBuffer.nZ = localStart[myCPULocation.nX][myCPULocation.nY][myCPULocation.nZ].nZ;
    myEndNoBuffer.nZ = localStart[myCPULocation.nX][myCPULocation.nY][myCPULocation.nZ + 1].nZ - 1;

    mySizeNoBuffer = myEndNoBuffer - myStartNoBuffer + delphi_integer(1);
    return mySizeNoBuffer;
}

//---------------------------------------------------------------------

SGrid<delphi_integer> CDelphiSpace::calcLocalGridStartEnd_fast()
{
    //allocate memory for dim array
    //dim array stores the size of all blocks
    SGrid<delphi_integer> *** localStart = new SGrid<delphi_integer> **[cpuDim.nX];
    SGrid<delphi_integer> *** localEnd = new SGrid<delphi_integer> **[cpuDim.nX];

    for (int x = 0; x < cpuDim.nX; x++)
    {
        localStart[x] = new SGrid<delphi_integer>*[cpuDim.nY];
        localEnd[x] = new SGrid<delphi_integer>*[cpuDim.nY];
        for (int y = 0; y < cpuDim.nY; y++)
        {
            localStart[x][y] = new SGrid<delphi_integer>[cpuDim.nZ];
            localEnd[x][y] = new SGrid<delphi_integer>[cpuDim.nZ];
        }
    }

    //evenly divide the atoms into sub divisions
    vector<delphi_real> xcoor, ycoor, zcoor;

    for (std::vector <CAtomPdb>::iterator it = delphipdb.begin(); it != delphipdb.end(); it++)
    {
        xcoor.push_back(it->getPose().nX);
        ycoor.push_back(it->getPose().nY);
        zcoor.push_back(it->getPose().nZ);
    }
    std::sort(xcoor.begin(), xcoor.end());
    std::sort(ycoor.begin(), ycoor.end());
    std::sort(zcoor.begin(), zcoor.end());

    vector<delphi_integer> xdiv, ydiv, zdiv;

    xdiv.push_back(0);
    for (int x = 1; x < cpuDim.nX; x++)
    {
        delphi_real div_coor = xcoor[global_iNatom / cpuDim.nX * x];
        delphi_integer div_grid = (div_coor - global_cOldMid.nX)*fScale + (global_iGrid.nX + 1) / 2;
        xdiv.push_back(div_grid);
    }
    xdiv.push_back(global_iGrid.nX);

    ydiv.push_back(0);
    for (int y = 1; y < cpuDim.nY; y++)
    {
        delphi_real div_coor = ycoor[global_iNatom / cpuDim.nY * y];
        delphi_integer div_grid = (div_coor - global_cOldMid.nY)*fScale + (global_iGrid.nY + 1) / 2;
        ydiv.push_back(div_grid);
    }
    ydiv.push_back(global_iGrid.nY);

    zdiv.push_back(0);
    for (int z = 1; z < cpuDim.nZ; z++)
    {
        delphi_real div_coor = zcoor[global_iNatom / cpuDim.nZ * z];
        delphi_integer div_grid = (div_coor - global_cOldMid.nZ)*fScale + (global_iGrid.nZ + 1) / 2;
        zdiv.push_back(div_grid);
    }
    zdiv.push_back(global_iGrid.nZ);

    myStartNoBuffer.nX = xdiv[myCPULocation.nX];
    myEndNoBuffer.nX = xdiv[myCPULocation.nX + 1] - 1;

    myStartNoBuffer.nY = ydiv[myCPULocation.nY];
    myEndNoBuffer.nY = ydiv[myCPULocation.nY + 1] - 1;

    myStartNoBuffer.nZ = zdiv[myCPULocation.nZ];
    myEndNoBuffer.nZ = zdiv[myCPULocation.nZ + 1] - 1;

    mySizeNoBuffer = myEndNoBuffer - myStartNoBuffer + delphi_integer(1);
    return mySizeNoBuffer;
}

//---------------------------------------------------------------------

bool CDelphiSpace::calcLocalGridStartEnd_naive()
{
    //First calculate the starting grid number without buffer
    myStartNoBuffer = { 0, 0, 0 };
    for (int x = 0; x < myCPULocation.nX; x++)
    {
        myStartNoBuffer.nX += localGridDim[x][myCPULocation.nY][myCPULocation.nZ].nX;
    }
    for (int y = 0; y < myCPULocation.nY; y++)
    {
        myStartNoBuffer.nY += localGridDim[myCPULocation.nX][y][myCPULocation.nZ].nY;
    }
    for (int z = 0; z < myCPULocation.nZ; z++)
    {
        myStartNoBuffer.nZ += localGridDim[myCPULocation.nX][myCPULocation.nY][z].nZ;
    }

    //Second, calculate the end grid number without buffer
    myEndNoBuffer = myStartNoBuffer + mySizeNoBuffer - delphi_integer(1);

    return true;
}

//---------------------------------------------------------------------

bool CDelphiSpace::assignBuffer()
{
    //find the largest atom radius
    max_radius = 0;

    for (int i = 0; i < iNatom; i++)
    {
        max_radius = max(max_radius, delphipdb[i].getRadius());
    }
    buffer = ceil((fExternRadius + max_radius * 2 + 1)*fScale);

    myStart = myStartNoBuffer - optMin(myStartNoBuffer, buffer);
    myEnd = myEndNoBuffer + optMin(global_iGrid - myEndNoBuffer - delphi_integer(1), buffer);
    mySize = myEnd - myStart + delphi_integer(1);

    //make grid size odd by changing expanding smaller edge
    if (mySize.nX % 2 == 0)
    {
        myStart.nX = myStart.nX - min(myStartNoBuffer.nX, delphi_integer(1));
    }
    if (mySize.nY % 2 == 0)
    {
        myStart.nY = myStart.nY - min(myStartNoBuffer.nY, delphi_integer(1));
    }
    if (mySize.nZ % 2 == 0)
    {
        myStart.nZ = myStart.nZ - min(myStartNoBuffer.nZ, delphi_integer(1));
    }
    mySize = myEnd - myStart + delphi_integer(1);

    //if still doesn't work, make grid size odd by changing expanding smaller edge
    if (mySize.nX % 2 == 0)
    {
        myEnd.nX = myEnd.nX + min(global_iGrid.nX - myEndNoBuffer.nX - delphi_integer(1), delphi_integer(1));
    }
    if (mySize.nY % 2 == 0)
    {
        myEnd.nY = myEnd.nY + min(global_iGrid.nY - myEndNoBuffer.nY - delphi_integer(1), delphi_integer(1));
    }
    if (mySize.nZ % 2 == 0)
    {
        myEnd.nZ = myEnd.nZ + min(global_iGrid.nZ - myEndNoBuffer.nZ - delphi_integer(1), delphi_integer(1));
    }
    mySize = myEnd - myStart + delphi_integer(1);

    //assign internal grid index without buffer
    myInternalStartNoBuffer = myStartNoBuffer - myStart;
    myInternalEndNoBuffer = myInternalStartNoBuffer + mySizeNoBuffer - delphi_integer(1);

    //update the grid handled by current process
    iGrid = mySize;
    fRMid = (optCast<delphi_real, delphi_integer>(global_iGrid) + 1.0) / 2.0;
    halfSize = optCast<delphi_real, delphi_integer>(iGrid) / fScale / 2.0;

    //the box center is changed to local box
    cOldMid = global_cOldMid + (optCast<delphi_real, delphi_integer>(myStart) * 2.0 + optCast<delphi_real, delphi_integer>(mySize) - optCast<delphi_real, delphi_integer>(global_iGrid)) / fScale / 2.0;
    myStartCoor = cOldMid - halfSize;
    myEndCoor = cOldMid + halfSize;
    myStartCoorNoBuffer = myStartCoor + (optCast<delphi_real, delphi_integer>(myStartNoBuffer - myStart) - 0.5) / fScale;
    myEndCoorNoBuffer = myEndCoor - (optCast<delphi_real, delphi_integer>(myEnd - myEndNoBuffer) + 0.5) / fScale;
}

//---------------------------------------------------------------------

void CDelphiSpace::split()
{
    cpuDim = calcCPUDim();
    myCPULocation = calcMyCPULocation();

    //calcGridDim();
    calcLocalGridStartEnd();
    assignBuffer();
};

#ifdef PARALLEL_MPI

void CDelphiSpace::verify()
{
    cout << "Verify parameters into space module " << endl;

    checkValue<delphi_integer>("natom");
    checkValue<bool>("iautocon");
    checkValue<delphi_real>("scale");
    checkValue<delphi_integer>("igrid");
    checkValue<delphi_integer>("nobject");
    checkValue<delphi_real>("repsout");
    checkValue<delphi_real>("repsin");
    checkValue<bool>("uniformdiel");
    checkValue<bool>("ionlymol");
    checkValue<bool>("isolv");
    checkValue<bool>("irea");
    checkValue<bool>("logs");
    checkValue<bool>("lognl");
    checkValue<bool>("isen");
    checkValue<bool>("isch");
    checkValue<delphi_real>("epsout");
    checkValue<delphi_real>("deblen");
    checkValue<delphi_real>("epsin");
    checkValue<bool>("epswrt");
    checkValue<bool>("isite");
    checkValue<bool>("ibem");
    checkValue<int>("ibctyp");

    checkString("epsnam"); //strEpsFile(pdc->getKey_constRef<string>("epsnam")),

    checkValue<bool>("isitsf");
    checkValue<delphi_real>("cutoff");
    checkValue<delphi_real>("sigma");
    checkValue<int>("inhomo");

    checkValue<delphi_real>("srfcut");
    checkValue<int>("gaussian");

    checkValue<int>("nmedia");
    checkValue<int>("numbmol");
    checkValue<int>("scrgfrm");

    checkValue<delphi_real>("rionst");
    checkValue<delphi_real>("exrad");
    checkValue<delphi_real>("rdmx");

    checkVector<delphi_real>("radprb"); //fRadPrb_v(pdc->getKey_constRef<vector <delphi_real> >("radprb")),
};

//---------------------------------------------------------------------

void CDelphiSpace::reallocate_global(delphi_integer usingCPU)
{
    ddm::Shared<size_t> globalSize;
    dplt::globalVec1D<delphi_integer> distSize(pdc->num_procs());
    ddm::barrier();

    pdc->deallocateGlobal1D<char>("idebmap");
    if (pdc->myid() == usingCPU) globalSize.set(bDebMap_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<char>("idebmap", globalSize.get());
    ddm::barrier();

    //ARGO: doing the 'pdc' thing just like its done with idebmap
    //zetaSurfMap_v(pdc->getKey_Ref< vector< char > >("zetaSurfMap")),
    pdc->deallocateGlobal1D<char>("zetaSurfMap");
    if (pdc->myid() == usingCPU) globalSize.set(zetaSurfMap_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<char>("zetaSurfMap", globalSize.get());
    ddm::barrier();

    //iEpsMap_v(pdc->getKey_Ref< vector< SGrid<delphi_integer> > >("iepsmp")),
    pdc->deallocateGlobal1D<SGrid<delphi_integer>>("iepsmp");
    if (pdc->myid() == usingCPU) globalSize.set(iEpsMap_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_integer>>("iepsmp", globalSize.get());
    ddm::barrier();

    //fGepsMap_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("gepsmp")),
    pdc->deallocateGlobal1D<SGrid<delphi_real>>("gepsmp");
    if (pdc->myid() == usingCPU) globalSize.set(fGepsMap_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_real>>("gepsmp", globalSize.get());
    ddm::barrier();

    //fGepsMap2_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("gepsmp2")),
    pdc->deallocateGlobal1D<SGrid<delphi_real>>("gepsmp2");
    if (pdc->myid() == usingCPU) globalSize.set(fGepsMap2_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_real>>("gepsmp2", globalSize.get());
    ddm::barrier();

    //fGDensityMap_v(pdc->getKey_Ref< vector< delphi_real> >("gdensity")),
    pdc->deallocateGlobal1D<delphi_real>("gdensity");
    if (pdc->myid() == usingCPU) globalSize.set(fGDensityMap_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<delphi_real>("gdensity", globalSize.get());
    ddm::barrier();

    //ibgrd_v(pdc->getKey_Ref< vector< SGrid<delphi_integer> > >("ibgrd")),
    pdc->deallocateGlobal1D<SGrid<delphi_integer>>("ibgrd");
    if (pdc->myid() == usingCPU) globalSize.set(ibgrd_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_integer>>("ibgrd", globalSize.get());
    ddm::barrier();

    //scspos_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("scspos")),
    pdc->deallocateGlobal1D<SGrid<delphi_real>>("scspos");
    if (pdc->myid() == usingCPU) globalSize.set(scspos_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_real>>("scspos", globalSize.get());
    ddm::barrier();

    //chgpos_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("chgpos")),
    pdc->deallocateGlobal1D<SGrid<delphi_real>>("chgpos");
    if (pdc->myid() == usingCPU) globalSize.set(chgpos_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_real>>("chgpos", globalSize.get());
    ddm::barrier();

    //crgatn_v(pdc->getKey_Ref< vector< delphi_integer > >("crgatn")),
    pdc->deallocateGlobal1D<delphi_integer>("crgatn");
    if (pdc->myid() == usingCPU) globalSize.set(crgatn_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<delphi_integer>("crgatn", globalSize.get());
    ddm::barrier();

    //nqgrdtonqass_v(pdc->getKey_Ref< vector< delphi_integer > >("nqgrdtonqass")),
    pdc->deallocateGlobal1D<delphi_integer>("nqgrdtonqass");
    if (pdc->myid() == usingCPU) globalSize.set(nqgrdtonqass_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<delphi_integer>("nqgrdtonqass", globalSize.get());
    ddm::barrier();

    //atmeps_v(pdc->getKey_Ref< vector< delphi_real > >("atmeps")),
    pdc->deallocateGlobal1D<delphi_real>("atmeps");
    if (pdc->myid() == usingCPU) globalSize.set(atmeps_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<delphi_real>("atmeps", globalSize.get());
    ddm::barrier();

    //atmcrg_v(pdc->getKey_Ref< vector< SGridValue<delphi_real> > >("atmcrg")),
    pdc->deallocateGlobal1D<SGridValue<delphi_real>>("atmcrg");
    if (pdc->myid() == usingCPU) globalSize.set(atmcrg_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<SGridValue<delphi_real>>("atmcrg", globalSize.get());
    ddm::barrier();

    //chrgv2_v(pdc->getKey_Ref< vector< SGridValue<delphi_real> > >("chrgv2")),
    pdc->deallocateGlobal1D<SGridValue<delphi_real>>("chrgv2");
    if (pdc->myid() == usingCPU) globalSize.set(chrgv2_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<SGridValue<delphi_real>>("chrgv2", globalSize.get());
    ddm::barrier();

    //scsnor_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("scsnor")),
    pdc->deallocateGlobal1D<SGrid<delphi_real>>("scsnor");
    if (pdc->myid() == usingCPU) globalSize.set(scsnor_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_real>>("scsnor", globalSize.get());
    ddm::barrier();

    //atsurf_v(pdc->getKey_Ref< vector< delphi_integer > >("atsurf")),
    pdc->deallocateGlobal1D<delphi_integer>("atsurf");
    if (pdc->myid() == usingCPU) globalSize.set(atsurf_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<delphi_integer>("atsurf", globalSize.get());
    ddm::barrier();

    //atndx_v(pdc->getKey_Ref< vector< delphi_integer > >("atndx"))
    pdc->deallocateGlobal1D<delphi_integer>("atndx");
    if (pdc->myid() == usingCPU) globalSize.set(atndx_v.size());
    ddm::barrier();
    pdc->allocateGlobal1D<delphi_integer>("atndx", globalSize.get());
    ddm::barrier();
}

//---------------------------------------------------------------------

void CDelphiSpace::write_to_global()
{
    size_t ncpus = pdc->num_procs();
    ddm::Shared<size_t> globalSize;
    ddm::Shared<delphi_real> globalSum;
    dplt::globalVec1D<delphi_integer> distSize(ncpus * ddm::size(), ddm::BLOCKED);
    dplt::globalVec1D<delphi_real> distSum(ncpus* ddm::size(), ddm::BLOCKED);
    dplt::globalVec1D<SGrid<delphi_real>> distGridSum(ncpus* ddm::size(), ddm::BLOCKED);
    size_t localStart, localSize;
    size_t iGridCube = global_iGrid.nX * global_iGrid.nY * global_iGrid.nZ;

    //*****************************************
    //**       Due with boundary points      **
    //*****************************************

    //Recalculate the following arrays on each CPU

    delphi_integer newBoundNum = 0;
    vector< SGrid<delphi_integer> > new_ibgrd_v;
    vector< SGrid<delphi_real> > new_scspos_v;
    vector< delphi_integer > new_atsurf_v;

    for (int i = 0; i < iBoundNum; i++)
    {
        if (checkNotBufferGrid(ibgrd_v[i]))
        {
            newBoundNum++;
            ibgrd_v[i] = convLocalGridToGlobal(ibgrd_v[i]);
            new_ibgrd_v.push_back(ibgrd_v[i]);
            scspos_v[i] = scspos_v[i];
            new_scspos_v.push_back(scspos_v[i]);
            if (!isitsf && !isite && !(isch&&scrgfrm != 0))
            {

            }
            else
            {
                new_atsurf_v.push_back(atsurf_v[i]);
            }
        }
    }

    //Collect the size of the arrays on each CPU
    localSize = newBoundNum;
    distSize[my_id] = localSize;
    ddm::barrier();
    if (my_id == 0)
    {
        delphi_integer sum = 0;
        for (int i = 0; i < ncpus; i++)
        {
            sum += distSize[i];
        }
        globalSize.set(sum);
        //work on "ibnum"
        pdc->writeGlobalVar<delphi_integer>("ibnum", sum);
    }
    ddm::barrier();
    localStart = 0;
    for (int i = 0; i < my_id; i++)
    {
        localStart += distSize[i];
    }

    ddm::barrier();

    //work on "ibgrd"
    pdc->deallocateGlobal1D<SGrid<delphi_integer>>("ibgrd");
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_integer>>("ibgrd", globalSize.get());
    pdc->writeGlobalVector1D<SGrid<delphi_integer>>("ibgrd", localStart, localSize, new_ibgrd_v);

    //work on "scspos"
    pdc->deallocateGlobal1D<SGrid<delphi_real>>("scspos");
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_real>>("scspos", globalSize.get());
    pdc->writeGlobalVector1D<SGrid<delphi_real>>("scspos", localStart, localSize, new_scspos_v);

    //work on "atsurf"
    pdc->deallocateGlobal1D<delphi_integer>("atsurf");
    ddm::barrier();
    if (!isitsf && !isite && !(isch&&scrgfrm != 0))
    {
        pdc->allocateGlobal1D<delphi_integer>("atsurf", 0);
    }
    else
    {
        pdc->allocateGlobal1D<delphi_integer>("atsurf", globalSize.get());
        pdc->writeGlobalVector1D<delphi_integer>("atsurf", localStart, localSize, new_atsurf_v);
    }

    //*****************************************
    //**     Due with charged grid points    **
    //*****************************************

    delphi_integer new_nqgrd = 0;
    vector< delphi_integer > new_nqgrdtonqass_v;
    vector< SGridValue<delphi_real> > new_chrgv2_v;

    for (int i = 0; i < nqgrd; i++)
    {
        if (checkNotBufferGrid(optCast<delphi_integer, delphi_real>(chrgv2_v[i].nGrid)))
        {
            new_nqgrd++;
            new_nqgrdtonqass_v.push_back(nqgrdtonqass_v[i]);
            chrgv2_v[i].nGrid = convLocalGridToGlobal(chrgv2_v[i].nGrid);
            new_chrgv2_v.push_back(chrgv2_v[i]);
        }
    }

    //Collect the size of the arrays on each CPU
    localSize = new_nqgrd;
    distSize[my_id] = localSize;
    ddm::barrier();
    if (my_id == 0)
    {
        delphi_integer sum = 0;
        for (int i = 0; i < ncpus; i++)
        {
            sum += distSize[i];
        }
        globalSize.set(sum);
        //work on "nqgrd"
        pdc->writeGlobalVar<delphi_integer>("nqgrd", sum);
    }
    ddm::barrier();

    localStart = 0;
    for (int i = 0; i < my_id; i++)
    {
        localStart += distSize[i];
    }

    ddm::barrier();

    //work on "nqgrdtonqass"
    pdc->deallocateGlobal1D<delphi_integer>("nqgrdtonqass");
    ddm::barrier();
    pdc->allocateGlobal1D<delphi_integer>("nqgrdtonqass", globalSize.get());
    pdc->writeGlobalVector1D<delphi_integer>("nqgrdtonqass", localStart, localSize, new_nqgrdtonqass_v);

    //work on "chrgv2"
    pdc->deallocateGlobal1D<SGridValue<delphi_real>>("chrgv2");
    ddm::barrier();
    pdc->allocateGlobal1D<SGridValue<delphi_real>>("chrgv2", globalSize.get());
    pdc->writeGlobalVector1D<SGridValue<delphi_real>>("chrgv2", localStart, localSize, new_chrgv2_v);


    //*****************************************
    //**       Due with charged atoms      **
    //*****************************************

    delphi_integer new_nqass = 0;
    vector< SGrid<delphi_real> > new_chgpos_v;
    vector< delphi_integer> new_crgatn_v;
    vector< delphi_real > new_atmeps_v;
    vector< SGridValue<delphi_real> > new_atmcrg_v;

    for (int i = 0; i < nqass; i++)
    {
        if (checkNotBufferCoor(chgpos_v[i]))
        {
            new_nqass++;
            new_chgpos_v.push_back(chgpos_v[i]);
            new_crgatn_v.push_back(crgatn_v[i]);
            new_atmeps_v.push_back(atmeps_v[i]);
            atmcrg_v[i].nGrid = convLocalGridToGlobal(atmcrg_v[i].nGrid);
            new_atmcrg_v.push_back(atmcrg_v[i]);
        }
    }

    //Collect the size of the arrays on each CPU
    localSize = new_nqass;
    distSize[my_id] = localSize;
    ddm::barrier();
    if (my_id == 0)
    {
        delphi_integer sum = 0;
        for (int i = 0; i < ncpus; i++)
        {
            sum += distSize[i];
        }
        globalSize.set(sum);
        //work on "nqass"
        pdc->writeGlobalVar<delphi_integer>("nqass", sum);
    }
    ddm::barrier();

    localStart = 0;
    for (int i = 0; i < my_id; i++)
    {
        localStart += distSize[i];
    }

    ddm::barrier();

    //work on "chgpos"
    pdc->deallocateGlobal1D<SGrid<delphi_real>>("chgpos");
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_real>>("chgpos", globalSize.get());
    pdc->writeGlobalVector1D<SGrid<delphi_real>>("chgpos", localStart, localSize, new_chgpos_v);

    //work on "crgatn"
    pdc->deallocateGlobal1D<delphi_integer>("crgatn");
    ddm::barrier();
    pdc->allocateGlobal1D<delphi_integer>("crgatn", globalSize.get());
    pdc->writeGlobalVector1D<delphi_integer>("crgatn", localStart, localSize, new_crgatn_v);

    //work on "atmeps"
    pdc->deallocateGlobal1D<delphi_real>("atmeps");
    ddm::barrier();
    pdc->allocateGlobal1D<delphi_real>("atmeps", globalSize.get());
    pdc->writeGlobalVector1D<delphi_real>("atmeps", localStart, localSize, new_atmeps_v);

    //work on "atmcrg"
    pdc->deallocateGlobal1D<SGridValue<delphi_real>>("atmcrg");
    ddm::barrier();
    pdc->allocateGlobal1D<SGridValue<delphi_real>>("atmcrg", globalSize.get());
    pdc->writeGlobalVector1D<SGridValue<delphi_real>>("atmcrg", localStart, localSize, new_atmcrg_v);

    //*****************************************
    //**       Due with idebmap and epsmap   **
    //*****************************************

    //work on "idebmap"
    pdc->deallocateGlobal1D<char>("idebmap");
    ddm::barrier();
    pdc->allocateGlobal1D<char>("idebmap", iGridCube);
    ddm::barrier();

    dplt::global1D<char> &global_idebmap = pdc->getGlobalObj<dplt::global1D<char>>("idebmap");
    for (int k = myInternalStartNoBuffer.nZ; k <= myInternalEndNoBuffer.nZ; k++)
    {
        for (int j = myInternalStartNoBuffer.nY; j <= myInternalEndNoBuffer.nY; j++)
        {
            vector<char> temp;
            size_t globalIndex = (k + myStart.nZ)*global_iGrid.nY*global_iGrid.nX + (j + myStart.nY)*global_iGrid.nX + myStart.nX + myInternalStartNoBuffer.nX;
            size_t localIndex = k*iGrid.nY*iGrid.nX + j*iGrid.nX + myInternalStartNoBuffer.nX;

            for (int i = myInternalStartNoBuffer.nX; i <= myInternalEndNoBuffer.nX; i++)
            {
                temp.push_back(bDebMap_v[localIndex]);
                localIndex++;
            }

            pdc->writeGlobalVector1D<char>("idebmap", globalIndex, mySizeNoBuffer.nX, temp);
        }
    }
    ddm::barrier();

    size_t global_start_z = myInternalStartNoBuffer.nZ + myStart.nZ;
    size_t global_start_y = myInternalStartNoBuffer.nY + myStart.nY;
    size_t global_start_x = myInternalStartNoBuffer.nX + myStart.nX;
    dplt::Dim3 idebmap_size = { mySizeNoBuffer.nZ,mySizeNoBuffer.nY,mySizeNoBuffer.nX };
    dplt::Dim3 idebmap_start = { global_start_z ,global_start_y, global_start_x };

    //work on "iepsmp"
    pdc->deallocateGlobal1D<SGrid<delphi_integer>>("iepsmp");
    ddm::barrier();
    pdc->allocateGlobal1D<SGrid<delphi_integer>>("iepsmp", iGridCube);
    ddm::barrier();
    dplt::global1D<SGrid<delphi_integer>> &global_iepsmp = pdc->getGlobalObj<dplt::global1D<SGrid<delphi_integer>>>("iepsmp");
    for (int k = myInternalStartNoBuffer.nZ; k <= myInternalEndNoBuffer.nZ; k++)
    {
        for (int j = myInternalStartNoBuffer.nY; j <= myInternalEndNoBuffer.nY; j++)
        {
            vector<SGrid<delphi_integer>> temp;
            size_t globalIndex = (k + myStart.nZ)*global_iGrid.nY*global_iGrid.nX + (j + myStart.nY)*global_iGrid.nX + myStart.nX + myInternalStartNoBuffer.nX;
            size_t localIndex = k*iGrid.nY*iGrid.nX + j*iGrid.nX + myInternalStartNoBuffer.nX;
            for (int i = myInternalStartNoBuffer.nX; i <= myInternalEndNoBuffer.nX; i++)
            {
                temp.push_back(iEpsMap_v[localIndex]);
                localIndex++;
            }
            pdc->writeGlobalVector1D<SGrid<delphi_integer>>("iepsmp", globalIndex, mySizeNoBuffer.nX, temp);
        }

    }
    ddm::barrier();

    //*****************************************
    //**        Due with variables           **
    //*****************************************

    //work on "qnet"
    distSum[my_id] = qnet;
    ddm::barrier();
    if (my_id == 0)
    {
        delphi_real sum = 0;
        for (int i = 0; i < ncpus; i++)
        {

            sum += distSum[i];
        }
        pdc->writeGlobalVar<delphi_real>("qnet", sum);
    }
    ddm::barrier();

    //work on "qmin" and "cqmin"
    distSum[my_id] = qmin;
    distGridSum[my_id] = convLocalGridToGlobal(cqmin);
    ddm::barrier();
    if (my_id == 0)
    {
        delphi_real sum = 0;
        delphi_real temp = 0;
        SGrid<delphi_real> gridSum = { 0, 0, 0 };
        SGrid<delphi_real> gridTemp = { 0, 0, 0 };
        for (int i = 0; i < ncpus; i++)
        {
            temp = distSum[i];
            gridTemp = distGridSum[i];
            sum += temp;
            gridSum = gridSum + gridTemp * temp;
        }
        if (abs(sum) > 1.e-6) gridSum = gridSum / sum;
        pdc->writeGlobalVar<delphi_real>("qmin", sum);
        pdc->writeGlobalVar<SGrid<delphi_real>>("cqmin", gridSum);
    }
    ddm::barrier();

    //work on "qplus" and "cqplus"
    distSum[my_id] = qplus;
    distGridSum[my_id] = convLocalGridToGlobal(cqplus);
    ddm::barrier();
    if (my_id == 0)
    {
        delphi_real sum = 0;
        delphi_real temp = 0;
        SGrid<delphi_real> gridSum = { 0, 0, 0 };
        SGrid<delphi_real> gridTemp = { 0, 0, 0 };
        for (int i = 0; i < ncpus; i++)
        {
            temp = distSum[i];
            gridTemp = distGridSum[i];
            sum += temp;
            gridSum = gridSum + gridTemp * temp;
        }
        if (sum > 1.e-6) gridSum = gridSum / sum;
        pdc->writeGlobalVar<delphi_real>("qplus", sum);
        pdc->writeGlobalVar<SGrid<delphi_real>>("cqplus", gridSum);
    }
    ddm::barrier();

}

#endif //PARALLEL_MPI
