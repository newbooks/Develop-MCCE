#include <unistd.h>
#include <iostream>
#include <cstddef>
#include <sstream>

#include "ddm/DDM.h"
#include "any.h"
#include <time.h>
#include <math.h>
#include <map>

#define bCheckValue false
#define VERBOSE false
#define PATTERN BLOCKCYCLIC(8)
#define BLOCK_SIZE 8

int ARRAY_SIZE=512;

using namespace std;
using uint = unsigned int;
  typedef std::map<std::string, any> mapType;
using std::cout;
using std::cin;
using std::endl;
using std::vector;


bool fillArray3D_Method_1 (ddm::NArray<double, 3> &matrix) 
{
    //std::vector<ddm::Future<ddm::NArray<double, 3>*>> fut_dest_end;
    using mT = ddm::GlobViewIter<double, ddm::BlockPattern<3, (ddm::MemArrange)1u, long int>, ddm::GlobMem<double, ddm::allocator::CollectiveAllocator<double> >, ddm::GlobPtr<double, ddm::BlockPattern<3, (ddm::MemArrange)1u, long int> >, ddm::GlobRef<double> >;
    vector<ddm::Future<mT>> fut_dest_end;
    for(size_t i = 0; i < matrix.extent(0); i++) 
      for(size_t j = 0; j < matrix.extent(1); ++j)
      {
            for(size_t k = 0; k < matrix.extent(2); ++k)
            {
                double temp=(i * 10000) + (j*100) + (k * 1);
                auto reg = matrix.sub<0>(i,1).sub<1>(j,1).sub<2>(k,1);
                fut_dest_end.push_back(ddm::copy_async(&temp, &temp+1, reg.begin()));
                if(fut_dest_end.size()>matrix.extent(2))
                (fut_dest_end.end()- matrix.extent(2))->get();
            }
            
        }
        
 
    for(vector<ddm::Future<mT>>::iterator it=fut_dest_end.end()- matrix.extent(2); it!=fut_dest_end.end(); it++)
    {
        it->get();
    }
    matrix.barrier();
    return false;
}


bool fillArray3D_Method_2 (ddm::NArray<double, 3> &matrix) 
{
    if(VERBOSE)cout << "filling Matrix size " << matrix.extent(0)<<"  "<< matrix.extent(1)<<"  "<< matrix.extent(2)<<endl;
    for(size_t i = 0; i < matrix.extent(0); i++) 
      for(size_t j = 0; j < matrix.extent(1); ++j){
          std::vector<double> temp (matrix.extent(2));          
          for(size_t k = 0; k < matrix.extent(2); ++k)
          {
              temp[k]=(i * 10000) + (j*100) + (k * 1);
          }
        auto reg = matrix.sub<0>(i,1).sub<1>(j,1).sub<2>(0,matrix.extent(2));
        ddm::copy(temp.data(), temp.data()+ matrix.extent(2), reg.begin());        
      }
    return false;
}

bool fillArray3D_Method_3 (ddm::NArray<double, 3> &matrix) 
{
    //std::vector<ddm::Future<ddm::NArray<double, 3>*>> fut_dest_end;
    for(size_t i = 0; i < matrix.extent(0); i++) {
      for(size_t j = 0; j < matrix.extent(1); ++j)
        for(size_t k = 0; k < matrix.extent(2); ++k){
            matrix[i][j][k] = (i * 10000) + (j*100) + (k * 1);
      }
    }
    matrix.barrier();
    return false;
}

bool fillArray3D_MPI_1 (ddm::NArray<double, 3> &matrix, int width) 
{
    if(VERBOSE)cout << "Filling: width " <<width;
    int offset=ddm::myid()* width;

    //cout<<"My id="<<ddm::myid()<<"   Offset="<<offset<<endl;
    
    for (int i = 0; i<width; ++i)
        for (int j = 0; j<matrix.extent(1); ++j)
            for (int k = 0; k<matrix.extent(1); ++k)
            {
                //auto reg = matrix.sub<0>(i+offset,1).sub<1>(j,1).sub<2>(k,1);
                //ddm::copy_async(reg.begin(), reg.end(),receive[i][j].data()+k).get();
                //cout << "i="<<i<<" j="<<j<<" k="<<k<< " V="<<receive[i][j][k]<<endl;
                matrix[i+offset][j][k] = ((i+offset) * 10000) + (j*100) + (k * 1);
            }
    ddm::barrier();
    if(VERBOSE)cout << "  Done..."<<endl;
    return false;
}

bool fillArray3D_MPI_2 (ddm::NArray<double, 3> &matrix, int width) 
{
    if(VERBOSE)cout << "Filling: width " <<width;
    int offset=ddm::myid()* width;
    
    for(size_t i = 0; i < width; i++) 
      for(size_t j = 0; j < matrix.extent(1); ++j){
          std::vector<double> temp (matrix.extent(2));          
          for(size_t k = 0; k < matrix.extent(2); ++k)
          {
              temp[k]=((i+offset) * 10000) + (j*100) + (k * 1);
          }
        auto reg = matrix.sub<0>(i+offset,1).sub<1>(j,1).sub<2>(0,matrix.extent(2));
        ddm::copy(temp.data(), temp.data()+ matrix.extent(2), reg.begin());        
      }
    return false;
}

bool copyArray3D_Method1 (ddm::NArray<double, 3> &matrix,std::vector<vector<vector<double>>> &receive, int globalSize, int localSize, int width) 
{
    int offset=ddm::myid()* width;

    //cout<<"My id="<<ddm::myid()<<"   Offset="<<offset<<endl;
    
    for (int i = 0; i<width; ++i)
        for (int j = 0; j<localSize; ++j)
            //for (int k = 0; k<copy_z; ++k)
            {
                // slice in all three dimensions
                auto reg = matrix.sub<0>(offset+i,1).sub<1>(j,1).sub<2>(0,localSize);
                ddm::copy(reg.begin(), reg.end(), receive[i][j].data());
            }
    return false;
}  

bool copyArray3D_Method2 (ddm::NArray<double, 3> &matrix,std::vector<vector<vector<double>>> &receive, int globalSize, int localSize, int width) 
{
    if(VERBOSE)cout << "Copying Method 2: width " <<width;
    int offset=ddm::myid()* width;

    //cout<<"My id="<<ddm::myid()<<"   Offset="<<offset<<endl;
    
    for (int i = 0; i<width; ++i)
        for (int j = 0; j<localSize; ++j)
            for (int k = 0; k<localSize; ++k)
            {
                //cout << "i="<<i<<" j="<<j<<" k="<<k<<endl;
                receive[i][j][k]=matrix[offset+i][j][k];
            }
    if(VERBOSE)cout << "  Done..."<<endl;
    return false;
}  

bool copyArray3D_Method3 (ddm::NArray<double, 3> &matrix,std::vector<vector<vector<double>>> &receive, int globalSize, int localSize, int width) 
{
    if(VERBOSE)cout << "Copying Method 2: width " <<width;
    int offset=ddm::myid()* width;

    //cout<<"My id="<<ddm::myid()<<"   Offset="<<offset<<endl;
    
    for (int i = 0; i<width; ++i)
        for (int j = 0; j<localSize; ++j)
            for (int k = 0; k<localSize; ++k)
            {
                auto reg = matrix.sub<0>(i+offset,1).sub<1>(j,1).sub<2>(k,1);
                ddm::copy(reg.begin(), reg.end(),receive[i][j].data()+k);
                //cout << "i="<<i<<" j="<<j<<" k="<<k<< " V="<<receive[i][j][k]<<endl;
            }
    ddm::barrier();
    if(VERBOSE)cout << "  Done..."<<endl;
    return false;
}

bool copyArray3D_Method4 (ddm::NArray<double, 3> &matrix,std::vector<vector<vector<double>>> &receive, int globalSize, int localSize, int width) 
{
    if(VERBOSE)cout << "Copying Method 2: width " <<width;
    int offset=ddm::myid()* width;

    //cout<<"My id="<<ddm::myid()<<"   Offset="<<offset<<endl;
    
    for (int i = 0; i<width; ++i)
        for (int j = 0; j<localSize; ++j)
            for (int k = 0; k<localSize; k=k+ARRAY_SIZE)
            {
                auto reg = matrix.sub<0>(i+offset,1).sub<1>(j,1).sub<2>(k,ARRAY_SIZE);
                ddm::copy(reg.begin(), reg.end(),receive[i][j].data()+k);
                //cout << "i="<<i<<" j="<<j<<" k="<<k<< " V="<<receive[i][j][k]<<endl;
            }
    ddm::barrier();
    if(VERBOSE)cout << "  Done..."<<endl;
    return false;
}

bool copyArray3D_Method5 (ddm::NArray<double, 3> &matrix,std::vector<vector<vector<double>>> &receive, int globalSize, int localSize, int width) 
{
    if(VERBOSE)cout << "Copying Method 2: width " <<width;
    int offset=ddm::myid()* width;

    //cout<<"My id="<<ddm::myid()<<"   Offset="<<offset<<endl;
    
    for (int i = 0; i<width; ++i)
        for (int j = 0; j<localSize; ++j)
            {                
                auto reg = matrix.sub<0>(i+offset,1).sub<1>(j,1).sub<2>(0,localSize);
                ddm::copy(reg.begin(), reg.end(),receive[i][j].data());
                //cout << "i="<<i<<" j="<<j<<" k="<<k<< " V="<<receive[i][j][k]<<endl;
            }
    ddm::barrier();
    if(VERBOSE)cout << "  Done..."<<endl;
    return false;
}


bool checkArray3D (ddm::NArray<double, 3> &matrix,std::vector<vector<vector<double>>> &receive, int localSize, int width)
{
    bool sucess=true;
    int offset=ddm::myid()* width;
      for (int i = 0; i<width; ++i)
        for (int j = 0; j<localSize; ++j)
            for (int k = 0; k<localSize; ++k){
                int value    = receive[i][j][k];
                int expected = ((i+offset) * 10000) + (j*100) + k;
                //cout << i << " " << j << " "<< k << " " << "E=" << expected << "  V=" << value << "Origin="<<matrix[i+offset][j][k]<<endl;
                if(expected != value)sucess=false;
            }
    ddm::Team::All().barrier();
    return sucess;
}

double calcAverTime (ddm::Array<double> &matrix_time, int num_units)
{
    // allocate local data structure
        std::vector<double> rec_t(num_units);
        
    // copy from region to local data structure
        ddm::copy(matrix_time.begin(), matrix_time.end(), rec_t.data());

        double sum_t=0;

        for (int i=0; i<num_units; ++i)
        {
            sum_t+=rec_t[i];
        }
        sum_t/=num_units;
        return sum_t;
}

template <typename T>
inline bool insert(mapType& objMap, std::string objStr,  const T & obj)
{
    return objMap.insert(std::make_pair(objStr, obj)).second;
    //objMap[objStr] = obj;
    //(*((objMap.insert(std::make_pair("A", (obj)))).first)).second;
    return false;
}

template <typename T>
inline bool insert1(mapType& objMap, std::string objStr,  const T & obj)
{
    return objMap.insert(std::make_pair(objStr, std::make_unique<T>())).second;
    //objMap[objStr] = obj;
    //(*((objMap.insert(std::make_pair("A", (obj)))).first)).second;
    return false;
}

/*template <typename T>
inline bool remove(mapType& objMap, std::string objStr)
{    
    mapType::iterator it = objMap.find(objStr);
    T objPtr = *(any_cast<T*>((*it).second));
    delete objPtr;
    return false;
}*/


template <typename T>
inline T& read(mapType& objMap, std::string objStr)
{
   mapType::iterator it = objMap.find("A");
   
   return *(any_cast<T*>((*it).second));
}

template<typename T> T * createInstance() 
{ 
    static T reg;
    return &reg;
}

template<typename T> T * deleteInstance(T& reg) 
{ 
    delete reg;
}

/*ddm::NArray<double, 3> * createInstance() 
{ 
    static ddm::NArray<double, 3> reg (ddm::Team::Null());
    return &reg;
}*/


int main(int argc, char* argv[])
{
    // Check the number of parameters
    if (argc < 2) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " Array_Size" << std::endl;
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        return 1;
    }
    ARRAY_SIZE=atoi(argv[1]);

  time_t start_s,init_s,filled_s;  
  ddm::init(&argc, &argv);
  size_t team_size    = ddm::Team::All().size();
  ddm::TeamSpec<3> teamspec_3d(team_size, 1,1);
  teamspec_3d.balance_extents();

  ddm::global_unit_t myid   = ddm::myid();
  if(0 == myid) 
  {
    start_s=clock();
    if(VERBOSE)cout <<"Array size="<<ARRAY_SIZE<<endl;
    }
  size_t num_units   = ddm::Team::All().size();

  
  size_t rows = ARRAY_SIZE ;
  size_t cols = ARRAY_SIZE ;
  size_t wids = ARRAY_SIZE ;

  size_t tilesize_x  = rows;
  size_t tilesize_y  = cols;
  size_t tilesize_z  = wids;

  //ddm::NArray<double, 3> matrix( ddm::SizeSpec<3>( rows,cols,wids), ddm::DistributionSpec<3>(ddm::CYCLIC,ddm::CYCLIC,ddm::NONE),ddm::Team::All(),teamspec_3d);                         
                         
  size_t matrix_size = rows * cols * wids;
  if(VERBOSE)cout<<"Asserting..."<<endl;

  if(VERBOSE)cout<<"Assert done..."<<endl;
  ddm::Array<double> matrix_t1(num_units);                         
  ddm::Array<double> matrix_t2(num_units);
  ddm::Array<double> matrix_t3(num_units);                         
  ddm::Array<double> matrix_t4(num_units);
  ddm::Array<double> matrix_t5(num_units);
  

  mapType dataMap;

  insert(dataMap,"A", createInstance <ddm::NArray<double, 3>>());
  
   //mapType::iterator it = dataMap.find("A");
   
   //ddm::NArray<double, 3> &matrix = *(any_cast<ddm::NArray<double, 3>*>((*it).second));
   
   ddm::NArray<double, 3> &matrix1= read<ddm::NArray<double, 3>>(dataMap,"A");
   //delete &matrix1;
   
   ddm::NArray<double, 3> &matrix= read<ddm::NArray<double, 3>>(dataMap,"A");
   
   matrix.allocate( ddm::SizeSpec<3>( rows,cols,wids), ddm::DistributionSpec<3>(ddm::CYCLIC,ddm::CYCLIC,ddm::NONE),teamspec_3d,ddm::Team::All());
  
  DDM_ASSERT(matrix_size == matrix.size());
  DDM_ASSERT(rows == matrix.extent(0));
  DDM_ASSERT(cols == matrix.extent(1));
  DDM_ASSERT(wids == matrix.extent(2));

  if (0 == myid) {
    init_s=clock();
    if(VERBOSE)cout << "Initialized: " << (init_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
    if(VERBOSE)cout << "Matrix size: " << rows << " x " << cols << " x " << wids << " == " << matrix_size << endl;
  }

  // Fill matrix
  if (0 == myid && false) // Filling method 2
  {
    //if(VERBOSE)
    cout << "Assigning matrix values..." << endl;
    fillArray3D_Method_1(matrix);      
    filled_s=clock();
    //if(VERBOSE)
    cout << "Filled: " << (filled_s-init_s)/double(CLOCKS_PER_SEC)*1000 << endl;

 }
 
     int width=int(ARRAY_SIZE/ddm::Team::All().size());
    
    if(true)
    {
        if (0 == myid ) 
      {
        //if(VERBOSE)
        cout << "Assigning matrix values..." << endl;
        }
        fillArray3D_MPI_2(matrix,width);
        
        if (0 == myid ) 
      {
            filled_s=clock();
        //if(VERBOSE)
        cout << "Filled: " << (filled_s-init_s)/double(CLOCKS_PER_SEC)*1000 << endl;
        }
    }
   
  // Units waiting for value initialization
  ddm::Team::All().barrier();
  
    // Testing Transfer Method 4
    if(VERBOSE)
    {
        if (0 == myid )cout << "Now start transfer method 4..."<<endl;
    }
            //Now doing the transfer
            int copy_size=ARRAY_SIZE;
    
            time_t check_s,transfer_s,valid_s;
            check_s=clock();  

            // allocate local data structure
            std::vector<vector<vector<double>>> receive(width,vector<vector<double>>(copy_size,vector<double>(copy_size)));
          
            //copyArray3D (matrix,receive, ARRAY_SIZE, copy_size,width);
            copyArray3D_Method5 (matrix,receive, ARRAY_SIZE, copy_size,width);
            ddm::barrier();
            
            transfer_s=clock();
            matrix_t1[myid]=(transfer_s-check_s)/double(CLOCKS_PER_SEC)*1000;

            bool sucess=checkArray3D(matrix,receive,copy_size,width);
            
            valid_s=clock();
            matrix_t2[myid]=(valid_s-transfer_s)/double(CLOCKS_PER_SEC)*1000;

            if(!sucess)cout << "Node " << myid << " Copy Failed!" << endl;    
                      
            ddm::Team::All().barrier();
          
            if(0==myid)
            {
              // calculate average time
                cout <<calcAverTime(matrix_t1,num_units)<< endl;          
            }
            ddm::Team::All().barrier();

  
  ddm::barrier();
  ddm::finalize();
  
  return 0;
  
}
