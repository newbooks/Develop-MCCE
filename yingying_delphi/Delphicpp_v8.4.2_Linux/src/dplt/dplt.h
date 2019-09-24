
#ifndef _DPLT_H_
#define _DPLT_H_

#include "ddm/DDM.h"
#include "any.h"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <algorithm>

#define DEBUG_DPLT false

#define BLOCKSIZE 1024

namespace dplt {
    typedef double real;
    typedef unsigned int uint;

    /**
    * Internal definition of 1 dimensional vectors, not for external use
    */
    template <typename T> using _vector1D = std::vector<T>;

    /**
    * Internal definition of 2 dimensional vectors, not for external use
    */
    template <typename T> using _vector2D = std::vector<std::vector<T>>;

    /**
    * Internal definition of 3 dimensional vectors, not for external use
    */
    template <typename T> using _vector3D = std::vector<std::vector<std::vector<T>>> ;

    /**
    * 1 dimensional arrays store in the distributed memory
    */
        template <typename T> using globalVec1D = ddm::Array<T>;

    template <typename T> using global1D = ddm::NArray<T,2>;

    /**
    * 2 dimensional arrays store in the distributed memory
    */
    template <typename T> using global2D = ddm::NArray<T, 2>;

    /**
    * 3 dimensional arrays store in the distributed memory
    */
    template <typename T> using global3D = ddm::NArray<T, 3>;

    /**
    * Global variables
    */
    template <typename T> using globalVar = ddm::GlobRef<T>;

    /**
    * 1 dimensional vectors for use in local memory
    */
    template <typename T> using vector1D = _vector1D<T>;

    /**
    * 2 dimensional vectors for use in local memory
    */
    template <typename T> struct vector2D : public std::vector<vector1D<T>>
    {
        vector2D(size_t i, size_t j) :_vector2D<T>(i, vector1D<T>(j)) {}
    };

    /**
    * 3 dimensional vectors for use in local memory
    */
    template <typename T> struct vector3D : public std::vector<vector2D<T>>
    {
        vector3D(size_t i, size_t j, size_t k) :std::vector<vector2D<T>>(i, vector2D<T>(j, k)) {}
    };

    /**
    * 2D coordinate type for convenience. Dim2
    */
    typedef struct Dim2
    {
        size_t x;
        size_t y;
        Dim2(size_t i, size_t j) : x(i), y(j) {};
    } Dim2;

    /**
    * 3D coordinate type for convenience. Dim3
    */
    typedef struct Dim3
    {
        size_t x;
        size_t y;
        size_t z;
        Dim3 (size_t i, size_t j, size_t k) : x(i), y(j), z(k){};
    } Dim3;

    /**
    * Return ID of the MPI process
    */
    inline size_t myid() { return ddm::myid(); };

    /**
    * Return the total number of the MPI processes
    */
    inline size_t num_procs() { return ddm::Team::All().size();; };

    /**
    * MPI barrier
    */
    inline void barrier() { ddm::barrier(); };

};



class DPLT 
{

public:


private:
    typedef std::map<std::string, any> mapType;

    mapType dataMap;
    size_t team_size;
    ddm::TeamSpec<3> teamspec_3d;
    ddm::TeamSpec<2> teamspec_2d;
    ddm::TeamSpec<1> teamspec_1d;


    /**
    * Create an instance of object type T and return the pointer to the new object
    */
    template<typename T> T * createInstance()
    {
        return new T();
    }

    /**
    * Create an instance of object type T with size and return the pointer to the new object
    */
    template<typename T> ddm::GlobMem<T> * createGlobalMem(size_t size)
    {
        return new ddm::GlobMem<T> (size);
    }


    /**
    * MPI finalize
    */
    inline void finalize() { ddm::finalize(); };



public:
    /**
    * Constructor
    */
    DPLT(int argc, char* argv[])
    {
        if(DEBUG_DPLT)std::cout << "DPLT init called" << std::endl;
        ddm::init(&argc, &argv);
        if (DEBUG_DPLT)std::cout << "DDM init done" << std::endl;
        team_size = ddm::Team::All().size();
        teamspec_3d = ddm::TeamSpec<3>(team_size, 1, 1);
        teamspec_2d = ddm::TeamSpec<2>(team_size, 1);
        teamspec_1d = ddm::TeamSpec<1>(team_size);
        teamspec_3d.balance_extents();
        teamspec_2d.balance_extents();
        teamspec_1d.balance_extents();
        if (DEBUG_DPLT)std::cout << "DPLT init done" << std::endl;
    }

    DPLT() { if (DEBUG_DPLT)std::cout << "DPLT default constructor called" << std::endl; } //char* messages[] ={"delphicpp_release"};  init(1, messages); }


    ~DPLT(void)
    {
        finalize();
    }


    /**
    * Create an object in the dataMap indexed by given string
    * Input: the name of the globale array as objKey
    */
    template <typename T>
    bool createGlobalObj(const std::string objKey);

    /**
    * Get reference to the object by given string
    * Input: the name of the globale array as objKey
    */
    template <typename T>
    T& getGlobalObj(const std::string objKey);

    /**
    * Get reference to the object by given string, without safety checks
    * Input: the name of the globale array as objKey
    */
    template <typename T>
    T& getGlobalObjUnsafe(const std::string objKey);

    /**
    * Create an object in the dataMap indexed by given string
    * Input: the name of the globale memory as objKey
    */
    template <typename N>
    bool createGlobalMemObj(const std::string objKey, size_t size);


    /**
    * Get reference to the globalMem object by given string
    * Input: the name of the globale memory as objKey
    */
    template <typename T>
    ddm::GlobMem<T> & getGlobalMemObj(const std::string objKey);

    /**
    * Create an object in the dataMap indexed by given string
    * Input: the name of the globale memory as objKey, size is one
    */
    template <typename T>
    bool createGlobalVarObj(const std::string objKey);

    /**
    * Get size of the globalMem object by given string
    * Input: the name of the globale memory as objKey
    */
    template <typename T>
    size_t sizeofGlobalMemObj(const std::string objKey);

    /**
    * Write the content of a vector to global memory
    * Input: the name of the globale memory as objKey, 1D vector
    */
    template <typename T>
    bool writeGlobalMem(const std::string objKey, dplt::vector1D<T> &input);

    /**
    * Read the content of global memory to a vector 
    * Input: the name of the globale memory as objKey, 1D vector
    */
    template <typename T>
    bool readGlobalMem(const std::string objKey, dplt::vector1D<T> &output);


    /**
    * Allocate 3 dimensional array in distributed memory
    * Input:<T> type of the elements in the 3D array
    * Input objKey: the name of the global array as objKey
    * Input size: the 3D size of the global array
    */
    template <typename T>
    bool allocateGlobal3D(const std::string objKey, dplt::Dim3 size);

    /**
    * Deallocate 3 dimensional array in distributed memory
    * Input:<T> type of the elements in the 3D array
    * Input objKey: the name of the global array as objKey
    */
    template <typename T>
    bool deallocateGlobal3D(const std::string objKey);

    /**
    * Allocate 23 dimensional array in distributed memory
    * Input:<T> type of the elements in the 2D array
    * Input objKey: the name of the global array as objKey
    * Input size: the 2D size of the global array
    */
    template <typename T>
    bool allocateGlobal2D(const std::string objKey, dplt::Dim2 size);


    /**
    * Allocate 1 dimensional array in distributed memory
    * Input:<T> type of the elements in the 1D array
    * Input objKey: the name of the global array as objKey
    * Input size: the length of the global array
    */
    template <typename T>
    bool allocateGlobal1D(const std::string objKey, size_t size);
    template <typename T>
    bool allocateGlobal1D_None(const std::string objKey, size_t size);

    /**
    * Allocate 1 dimensional array in distributed memory
    * Input:<T> type of the elements in the 1D array
    * Input objKey: the name of the global array as objKey
    */
    template <typename T>
    bool deallocateGlobal1D(const std::string objKey);


    /**
    * Allocate a shared variable in distributed memory
    * Input:<T> type of the elements 
    * Input objKey: the name of the global array as objKey
    */
    template <typename T>
    bool allocateGlobalVar(const std::string objKey);



    /**
    * Get the size of a 1D matrix in distributed memory
    * Input:<T> type of the elements
    * Input objKey: the name of the global array as objKey
    */
    template <typename T>
    inline size_t sizeofGlobal1D(const std::string objKey);

    /**
    * Get the size of a 3Dmatrix in distributed memory
    * Input:<T> type of the elements
    * Input objKey: the name of the global array as objKey
    */
    template <typename T>
    inline size_t sizeofGlobal3D(const std::string objKey);

    /**
    * Read from distributed 3D array to a 3D vector
    * Input:<T> type of the elements in the 3D array
    * Input: the name of the globale array as objKey
    * Input: the upper left (smallest) corner of the 3D global array
    * Input: the size of inquiry in 3 dimensions
    */
    template <typename T>
    bool readGlobalVector3D(const std::string objKey, dplt::Dim3 start, dplt::Dim3 size, dplt::vector3D<T> &receive);


    /**
    * Write to distributed 3D array from a local 3D array
    * Input:<T> type of the elements in the 3D array
    * Input: the name of the globale array as objKey
    * Input: the upper left (smallest) corner of the 3D global array
    * Input: the size of inquiry in 3 dimensions
    */
    template <typename T>
    bool writeGlobalVector3D(const std::string objKey, dplt::Dim3 start, dplt::Dim3 size, dplt::vector3D<T> &data);

    /**
    * Read from distributed 1D array to a 1D vector
    * Input:<T> type of the elements in the 1D array
    * Input: the name of the globale array as objKey
    * Input: the starting element of the 1D global array
    * Input: the size of inquiry in 1 dimension
    */
    template <typename T>
    bool readGlobalVector1D(const std::string objKey, size_t start, size_t size, dplt::vector1D<T> &receive);
    template <typename T>
    bool readGlobalVector1D_None(const std::string objKey, size_t start, size_t size, dplt::vector1D<T> &receive);
        
        template <typename T>
        T readGlobalVector1D(const std::string objKey, size_t index);

    /**
    * Write to distributed 1 array from a local 1D array
    * Input:<T> type of the elements in the 1D array
    * Input: the name of the globale array as objKey
    * Input: the starting element of the 1D global array
    * Input: the size of inquiry in 1 dimension
    */
    template <typename T>
    bool writeGlobalVector1D(const std::string objKey, size_t start, size_t size, dplt::vector1D<T> &data);
    template <typename T>
    bool writeGlobalVector1D_None(const std::string objKey, size_t start, size_t size, dplt::vector1D<T> &data);
    

    /**
    * Read from distributed memory variable to local variable
    * Input:<T> type of the elements
    * Input: the name of the globale array as objKey
    * Input: the reference of local variable
    */
    template <typename T>
    T readGlobalVar(const std::string objKey);


    /**
    * Write to distributed memory variable from local variable
    * Input:<T> type of the elements
    * Input: the name of the globale array as objKey
    * Input: the reference of local variable
    */
    template <typename T>
    bool writeGlobalVar(const std::string objKey, T &send);

    /**
    * Create a char array in the dataMap indexed by given string
    * Input: the name of the globale memory as objKey
    */
    inline bool createGlobalChar(const std::string objKey, size_t size);

    /**
    * Write to distributed memory variable from local string
    * Input: the name of the globale string as objKey
    * Input: the reference of local string
    */

    inline bool writeGlobalChar(const std::string objKey, const std::string &send);

    /**
    * Read from distributed memory variable to local string
    * Input: the name of the globale string as objKey
    * Input: the reference of local string
    */

    inline std::string readGlobalChar(const std::string objKey);

    /**
    * Create a string vector in the dataMap indexed by given string
    * Input: the name of the globale memory as objKey
    */
    inline bool createGlobal1DString(const std::string objKey, std::vector<std::string> & str);



    /**
    * Write to distributed memory variable from local string vector
    * Input: the name of the globale string as objKey
    * Input: the reference of local string
    */
    inline bool writeGlobal1DString(const std::string objKey, std::vector<std::string> & str);

    /**
    * Read from distributed memory variable to local string vector
    * Input: the name of the globale string as objKey
    * Input: the reference of local string
    */
    inline bool readGlobal1DString(const std::string objKey, std::vector<std::string> & str);


    /**
    * Return ID of the MPI process
    */
    inline size_t myid() { return ddm::myid(); };

    /**
    * Return the total number of the MPI processes
    */
    inline size_t num_procs() { return team_size; };

    /**
    * MPI barrier
    */
    inline void barrier() { ddm::barrier(); };
    
};



inline bool DPLT::createGlobalChar(const std::string objKey, size_t size)
{
    ddm::Shared<size_t> globalSize;
    if (myid() == 0) globalSize.set(size);
    ddm::barrier();
    bool success = dataMap.insert(std::make_pair(objKey, createGlobalMem<char>(globalSize.get()))).second;
    DDM_ASSERT_MSG(success, ("Failed to create string object " + objKey + " in global memory map!"));
    ddm::barrier();
    return true;
}

inline bool DPLT::writeGlobalChar(const std::string objKey, const std::string &send)
{
    ddm::GlobMem<char>& gString = getGlobalMemObj<char>(objKey);

    size_t i = 0;
    for (std::string::const_iterator it = send.begin(); it != send.end(); it++)
    {
        gString.put_value(*it, i);
        i++;
    }
    gString.flush_all();
    return true;
}


inline std::string DPLT::readGlobalChar(const std::string objKey)
{
    ddm::GlobMem<char>& gString = getGlobalMemObj<char>(objKey);
    std::string receive;
    size_t j = gString.size();
    char temp_char;
    for (size_t i = 0; i<j; i++)
    {
        gString.get_value(&temp_char, i);
        receive += temp_char;
    }
    return receive;
}


template <typename T>
bool DPLT::createGlobalObj(const std::string objKey)
{
    bool success = dataMap.insert(std::make_pair(objKey, createInstance <T>())).second;
    //DDM_ASSERT_MSG(success, ("Failed to create object " + objKey + " in global memory map!"));
    ddm::barrier();
    return success;
}


template <typename T>
inline T& DPLT::getGlobalObj(const std::string objKey)
{
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("object " + objKey + " not found in global memory map!"));
    return *(any_cast<T*>((*it).second));
}

template <typename T>
inline T& DPLT::getGlobalObjUnsafe(const std::string objKey)
{
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("object " + objKey + " not found in global memory map!"));
    return *(any_cast_unsafe<T*>((*it).second));
}


template <typename T>
bool DPLT::createGlobalMemObj(const std::string objKey, size_t size)
{

    bool success = dataMap.insert(std::make_pair(objKey, createGlobalMem<T>(size))).second;
    DDM_ASSERT_MSG(success, ("Failed to create object " + objKey + " in global memory map!"));
    ddm::barrier();
    return success;
}


template <typename T>
inline ddm::GlobMem<T> & DPLT::getGlobalMemObj(const std::string objKey)
{
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("object " + objKey + " not found in global memory map!"));
    return *(any_cast_unsafe<ddm::GlobMem<T>*>((*it).second));
}

template <typename T>
inline size_t DPLT::sizeofGlobalMemObj(const std::string objKey)
{
    ddm::GlobMem<T>& gMem = getGlobalMemObj<T>(objKey);
    return gMem.size();

}

template <typename T>
bool DPLT::writeGlobalMem(const std::string objKey, dplt::vector1D<T> &input)
{
    ddm::GlobMem<T>& gMem= getGlobalMemObj<T>(objKey);
    DDM_ASSERT_EQ(input.size(), gMem.size(), ("Size of GlobalMemobj " + objKey + " does not match the input size!"));
    size_t i = 0;
    for (typename dplt::vector1D<T>::iterator it = input.begin(); it != input.end(); it++)
    {
        gMem.put_value(*it, i);
        i++;
    }
    return true;
}


template <typename T>
bool DPLT::readGlobalMem(const std::string objKey, dplt::vector1D<T> &output)
{
    ddm::GlobMem<T>& gMem = getGlobalMemObj<T>(objKey);
    
    DDM_ASSERT_EQ(output.size(), gMem.size(), ("Size of GlobalMemobj " + objKey + " does not match the output size!"));
    size_t i = 0;
    for (typename dplt::vector1D<T>::iterator it = output.begin(); it != output.end(); it++)
    {
        gMem.get_value(it, i);
        i++;
    }
    return true;
}

template <typename T>
bool DPLT::allocateGlobal3D(const std::string objKey, dplt::Dim3 size)
{
    if (DEBUG_DPLT)std::cout << "*    allocating 3D matrix " << objKey << " \n";
    if (size.x==0 && size.y==0 && size.z==0)return true;
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("allocateGlobal3D " + objKey + " key not found!"));
    ddm::NArray<T, 3> &matrix = *(any_cast<ddm::NArray<T, 3>*>((*it).second));
    matrix.allocate(ddm::SizeSpec<3>(size.x, size.y, size.z), ddm::DistributionSpec<3>(ddm::CYCLIC, ddm::CYCLIC, ddm::NONE), teamspec_3d, ddm::Team::All());
    ddm::barrier();
    return true;
}

template <typename T>
bool DPLT::deallocateGlobal3D(const std::string objKey)
{
    if (DEBUG_DPLT)std::cout << "*    deallocating 3D matrix " << objKey << " \n";
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("deallocateGlobal3D " + objKey + " key not found!"));
    ddm::NArray<T, 3> &matrix = *(any_cast<ddm::NArray<T, 3>*>((*it).second));
    if (matrix.size() == 0)return true;
    matrix.deallocate();
    ddm::barrier();
    return true;
}

template <typename T>
inline size_t DPLT::sizeofGlobal3D(const std::string objKey)
{
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("allocateGlobal1D " + objKey + " key not found!"));
    dplt::global3D<T> &matrix = *(any_cast<dplt::global1D<T>*>((*it).second));
    return matrix.size();

}


template <typename T>
bool DPLT::writeGlobalVector3D(const std::string objKey, dplt::Dim3 start, dplt::Dim3 size, dplt::vector3D<T> &data)
{
    if (size.x == 0 && size.y == 0 && size.z == 0)return true;
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("writeGlobalVector3D " + objKey + " key not found!"));
    dplt::global3D<T> &matrix = *(any_cast<dplt::global3D<T>*>((*it).second));

    //REPORT_ERROR
    DDM_ASSERT_LE(size.x + start.x, matrix.extent(0), ("When writting " + objKey + ", X Range too large!"));
    DDM_ASSERT_LE(size.y + start.y, matrix.extent(1), ("When writting " + objKey + ", Y Range too large!"));
    DDM_ASSERT_LE(size.z + start.z, matrix.extent(2), ("When writting " + objKey + ", Z Range too large!"));

    for (size_t i = 0; i < size.x; ++i)
        for (size_t j = 0; j < size.y; ++j)
        {
            auto reg = ((matrix.template sub<0>(i+start.x,1)).template sub<1>(j+start.y,1)).template sub<2>(start.z, size.z);
            ddm::copy(data[i][j].data(), data[i][j].data() + size.z, reg.begin());
        }
    //ddm::barrier();
    return true;
}


template <typename T>
bool DPLT::readGlobalVector3D(const std::string objKey, dplt::Dim3 start, dplt::Dim3 size, dplt::vector3D<T> &receive)
{
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("readGlobalVector3D " + objKey + " key not found!"));
    dplt::global3D<T> &matrix = *(any_cast<dplt::global3D<T>*>((*it).second));

    DDM_ASSERT_LE(size.x + start.x, matrix.extent(0), ("When reading " + objKey + ", X Range too large!"));
    DDM_ASSERT_LE(size.y + start.y, matrix.extent(1), ("When reading " + objKey + ", Y Range too large!"));
    DDM_ASSERT_LE(size.z + start.z, matrix.extent(2), ("When reading " + objKey + ", Z Range too large!"));

    for (size_t i = 0; i < size.x; i++)
        for (size_t j = 0; j < size.y; ++j)
        {
            auto reg = ((matrix.template sub<0>(i + start.x, 1)).template sub<1>(j + start.y, 1)).template sub<2>(start.z, size.z);
            //if ( reg.is_local(ddm::Team::All().size()) )
            {
                //receive[i][j]=std::vector<T>(reg.begin()+ start.z,sizeof(T)*size.z);
            }
            //else
            {                
                ddm::copy(reg.begin(), reg.end(), receive[i][j].data());
            }
        }
    return true;
}

template <typename T>
bool DPLT::allocateGlobal1D_None(const std::string objKey, size_t size)
{
    if (size == 0)
    {
        return true;
    }
    //if (size == 1) size++;
    if (DEBUG_DPLT)std::cout << "Allocating " <<objKey<< " on node " << myid() << "  size="<<size<<std::endl;
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("allocateGlobal1D " + objKey + " key not found!"));
    dplt::global1D<T> &matrix = *(any_cast<dplt::global1D<T>*>((*it).second));
    matrix.allocate(size, ddm::DistributionSpec<1>(ddm::NONE), ddm::Team::All());
    return true;
}

template <typename T>
bool DPLT::allocateGlobal1D(const std::string objKey, size_t size)
{
    if (size == 0)
    {
        return true;
    }
    //if (size == 1) size++;
    if (DEBUG_DPLT)std::cout << "Allocating " << objKey << " on node " << myid() << "  size=" << size << std::endl;
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("allocateGlobal1D " + objKey + " key not found!"));
    dplt::global2D<T> &matrix = *(any_cast<dplt::global2D<T>*>((*it).second));
    matrix.allocate(ddm::SizeSpec<2>(size/BLOCKSIZE +1, BLOCKSIZE), ddm::DistributionSpec<2>(ddm::BLOCKCYCLIC(64), ddm::NONE), teamspec_2d, ddm::Team::All());
    std::string sizeStr = objKey;
    sizeStr += "Size";
    createGlobalObj<ddm::Shared<size_t>>(sizeStr);
    writeGlobalVar<size_t>(sizeStr, size);
    return true;
}

template <typename T>
bool DPLT::deallocateGlobal1D(const std::string objKey)
{
    if (DEBUG_DPLT)std::cout << "Deallocating " << objKey << " on node " << myid() << std::endl;
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("allocateGlobal1D " + objKey + " key not found!"));
    dplt::global2D<T> &matrix = *(any_cast<dplt::global2D<T>*>((*it).second));
    if (matrix.size() == 0)return true;
    matrix.deallocate();
    std::string sizeStr = objKey;
    sizeStr += "Size";
    size_t zero = 0;
    writeGlobalVar<size_t>(sizeStr, zero);
    
    return true;
}

template <typename T>
inline size_t DPLT::sizeofGlobal1D(const std::string objKey)
{
    std::string sizeStr = objKey;
    sizeStr += "Size";
    mapType::iterator it = dataMap.find(sizeStr);
    bool success = (it != dataMap.end()) ? true : false;
    //DDM_ASSERT_MSG(success, ("sizeofGlobal1D " + objKey + " key not found!"));
    if (success) 
    {
        return readGlobalVar<size_t>(sizeStr);
    }
    else
    {
        if(DEBUG_DPLT)
        {
            std::cout << "Size of " << objKey << " not found!" << std::endl;
        }
    }
    return size_t(0);
    
    
       /* mapType::iterator it = dataMap.find(objKey);
        bool success = (it != dataMap.end()) ? true : false;
        DDM_ASSERT_MSG(success, ("sizeofGlobal1D " + objKey + " key not found!"));
        dplt::global2D<T> &matrix = *(any_cast<dplt::global2D<T>*>((*it).second));
        return matrix.size();
    */
}

template <typename T>
bool DPLT::readGlobalVector1D_None(const std::string objKey, size_t start, size_t size, dplt::vector1D<T> &receive)
{

    return true;
}

template <typename T>
bool DPLT::readGlobalVector1D(const std::string objKey, size_t start, size_t size, dplt::vector1D<T> &receive)
{
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("readGlobalVector1D " + objKey + " key not found!"));
    dplt::global2D<T> &matrix = *(any_cast<dplt::global2D<T>*>((*it).second));

    //DDM_ASSERT_LE(start + size, matrix.size(), ("When reading 1D array " + objKey + ", Range too large!"));
    if (size == 0)
    {
        if (DEBUG_DPLT) std::cout << "skip reading " << objKey << " array size is zero! " << std::endl;
        return true;
    }

    if (receive.size() < size)
    {
        if (DEBUG_DPLT) std::cout << "Resizing " << objKey << " local array size is too small! " << std::endl;
        receive.resize(size);
    }
        size_t receive_start = 0;
        for (long i = start/BLOCKSIZE; i < (start+size-1)/BLOCKSIZE + 1; i++)
        {
            size_t x_start = std::max((long(start)-i*BLOCKSIZE), long(0));
            size_t x_end = std::min((start+size-i*BLOCKSIZE), size_t(BLOCKSIZE));
            size_t x_size = x_end - x_start;
            auto reg = ((matrix.template sub<0>(i, 1)).template sub<1>(x_start, x_size));
            
            ddm::copy(reg.begin(), reg.end(), receive.data() + receive_start);
            receive_start += x_size;
        }
    if (DEBUG_DPLT)std::cout << "*    read global1D " << objKey << " on node " << myid() << " size = " << size << "    *\n";
    return true;
}

template <typename T>
T DPLT::readGlobalVector1D(const std::string objKey, size_t index)
{
    T return_val;
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    //if (DEBUG_DPLT)std::cout << "*    read global1D by each " << objKey << " on node " << myid() << " index = " << index << "    *\n";
    DDM_ASSERT_MSG(success, ("readGlobalVector1D " + objKey + " key not found!"));
    dplt::global1D<T> &matrix = *(any_cast<dplt::global1D<T>*>((*it).second));

    DDM_ASSERT_LE(index, matrix.size(), ("When reading 1D array " + objKey + ", Range too large!"));
    
    size_t x = index/BLOCKSIZE;
    size_t y = index - x*BLOCKSIZE;

    return_val = matrix[x][y];

    return return_val;
}


template <typename T>
bool DPLT::writeGlobalVector1D_None(const std::string objKey, size_t start, size_t size, dplt::vector1D<T> &data)
{
    if (size == 0)return true;
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("writeGlobalVector1D " + objKey + " key not found!"));
    dplt::global1D<T> &matrix = *(any_cast<dplt::global1D<T>*>((*it).second));

    DDM_ASSERT_LE(start + size, matrix.size(), ("When writing 1D array " + objKey + ", Range too large!"));


    ddm::copy(data.data(), data.data() + size, matrix.begin()+start);
    if (DEBUG_DPLT)std::cout << "*    write global1D " << objKey << " on node " << myid() << " size = "<< size<< "    *\n";
    //ddm::barrier();
    return true;
}

template <typename T>
bool DPLT::writeGlobalVector1D(const std::string objKey, size_t start, size_t size, dplt::vector1D<T> &data)
{
    if (size == 0)return true;
    mapType::iterator it = dataMap.find(objKey);
    bool success = (it != dataMap.end()) ? true : false;
    DDM_ASSERT_MSG(success, ("writeGlobalVector1D " + objKey + " key not found!"));
    dplt::global2D<T> &matrix = *(any_cast<dplt::global2D<T>*>((*it).second));

    //DDM_ASSERT_LE(start + size, matrix.size(), ("When writing 1D array " + objKey + ", Range too large!"));


    size_t data_start = 0;
    for (long i = start/BLOCKSIZE; i < (start+size-1)/BLOCKSIZE + 1; i++)
    {
        size_t x_start = std::max((long(start)-i*BLOCKSIZE), long(0));
        size_t x_end = size_t(std::min( long((start+size)-long(i*BLOCKSIZE)), long(BLOCKSIZE)));
        size_t x_size = x_end - x_start;
    
                if(DEBUG_DPLT) std::cout << objKey << "  start:" << start << " size: " << size << " x_start: " << x_start << 
                          " x_end:" << x_end << " x_size: " << x_size << " i: " << i << " node: " << myid() <<
              "ArraySize: " << DPLT::sizeofGlobal1D<T>(objKey) << std::endl;
              auto reg = ((matrix.template sub<0>(i, 1)).template sub<1>(x_start, x_size));
        
        ddm::copy(data.data()+data_start, data.data()+data_start+x_size, reg.begin());
        data_start += x_size;
    }

    if (DEBUG_DPLT)std::cout << "*    write global1D " << objKey << " on node " << myid() << " size = " << size << "    *\n";
    //ddm::barrier();
    return true;
}

template <typename T>
bool DPLT::createGlobalVarObj(const std::string objKey)
{

}


template <typename T>
T DPLT::readGlobalVar(const std::string objKey)
{

    ddm::Shared<T> &sVar = getGlobalObj <ddm::Shared<T>>(objKey);
    T receive = sVar.get();
    if (DEBUG_DPLT)std::cout << "*    readGlobalVar " << objKey << " = " << receive << " on node " << myid() << "    *\n";
    return receive;

}


template <typename T>
bool DPLT::writeGlobalVar(const std::string objKey, T &send) 
{
    ddm::Shared<T> &sVar = getGlobalObj <ddm::Shared<T>>(objKey);
    sVar.set(send);
    if (DEBUG_DPLT)std::cout << "*    writeGlobalVar " << objKey << " = " << send << " on node " << myid() << "    *\n";
    return true;

}


inline bool DPLT::createGlobal1DString(const std::string objKey, std::vector<std::string> & str)
{
    ddm::Shared<size_t> globalSize;
    if (myid() == 0) globalSize.set(str.size());
    ddm::barrier();
    size_t str_size = static_cast<size_t>(globalSize.get());
    ddm::barrier();
    if (DEBUG_DPLT)std::cout << "Allocated " << objKey << " on node " << myid() << "  size=" << str_size << std::endl;
    if (str_size == 0)
    {
        return true;
    }
    
    for (size_t i = 0; i < str_size; i++)
    {
        if (myid() == 0) globalSize.set(str[i].size());
        ddm::barrier();
        createGlobalChar(objKey+ std::to_string(i), globalSize.get());
        ddm::barrier();
    }
        
    return true;
}

inline bool DPLT::writeGlobal1DString(const std::string objKey, std::vector<std::string> & str)
{
    if (str.size() == 0)
    {
        return true;
    }
    if (DEBUG_DPLT)std::cout << "*    writeGlobal1DString " << objKey << " on node " << myid() << "\n";

    for (size_t i = 0; i < str.size(); i++)
    {
        writeGlobalChar(objKey + std::to_string(i), str[i]);
    }
    
    return true;
}

inline bool DPLT::readGlobal1DString(const std::string objKey, std::vector<std::string> & str)
{
    std::vector<std::string> result;
    size_t i = 0;
    bool success = true;
    while (true)
    {
        std::string objKeyi = objKey + std::to_string(i);
        mapType::iterator it = dataMap.find(objKeyi);
        bool success = (it != dataMap.end()) ? true : false;
        if (!success) break;
        result.push_back(readGlobalChar(objKeyi));
        i++;
    }
    str = result;
    if (DEBUG_DPLT)std::cout << "*    readGlobal1DString " << objKey << " on node " << myid() << "\n";
    return true;
}


#endif  //DPLT_H_
