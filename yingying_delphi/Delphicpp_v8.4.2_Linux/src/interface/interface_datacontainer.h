/**
 * @file interface_datacontainer.h
 * @brief interface IDataContainer
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * This file declares the interface IDataContainer, the base class for classes such as CDelphiData.
 * Direct realizing an object of this base class is forbidden, but a point to this class can be used to point to
 * an instance of its derived class via polymorphism to provide unified access to various data containers.
 */

#ifndef IDATACONTAINER_H_
#define IDATACONTAINER_H_

#include <iostream> 
#include <string.h>      // STL::string 
#include <vector>        // STL::vector
#include <map>           // STL::map
#include <boost/any.hpp> // boost::any

#include "environment.h"
#include "../io/io_datatype.h"
#include "interface_exceptions.h"

#ifdef PARALLEL_MPI
#include "../dplt/dplt.h"
#include <cctype>
#endif

using namespace std;
using boost::any_cast;

typedef map<string, boost::any> DataMap;
//typedef map<string, any> DataMap;

//-----------------------------------------------------------------------//
#ifdef PARALLEL_MPI
class IDataContainer:virtual public DPLT
{
   protected:   

      /**
       * data map containing variables to be read/modified among many other classes. Each entry contains a
       * key (string type) and associated value (boost::any type)
       */
      DataMap myData;
	               	                  
      /**
       * member function to set above data map. must be virtual since only the derived class knows the
       * contents to be written in the data container.
       */
	   virtual void setMap() = 0;

	public:

	   /**
	    * constructor
	    */
	  IDataContainer()
	  {
          #ifdef DEBUG_OBJECT
	      cout << endl;
	      cout << "****************************************************************\n";
	      cout << "*               IDataContainer is constructed                  *\n";
	      cout << "****************************************************************\n";
          #endif
	  };
	   
      /**
       * virtual destructor allowing an instance of a derived class can be deleted properly
       * through a pointer to this base class
       */
	   virtual ~IDataContainer()
	   {
	      DataMap().swap(myData);

          #ifdef DEBUG_OBJECT
          cout << endl;
          cout << "****************************************************************\n";
          cout << "*                IDataContainer is destroyed                   *\n";
          cout << "****************************************************************\n";
          #endif
	   };

      /**
       * Function to check if the user-specified key exists or not in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    true if the key exists, false otherwise
       */
	   bool keyExists(const string &strKey);
         
	   //----------read-only key content in the map
	   /**
	    * Template function to get a constant reference to an entry of the data container
	    *
	    * @param[in] strKey The key to be searched in the data container
	    * @return    a constant reference to an entry if the key is found. Otherwise, an exception is thrown.
	    */
	   template <class T> const T& getKey_constRef(const string& strKey);
	   	   
	   /**
	    * Template function to get a constant pointer pointing to the data of a vector-type entry in the data container
	    *
	    * @param[in] strKey The key to be searched in the data container
	    * @return    a constant pointer if the key is found. Otherwise, an exception is thrown.
	    */
	   template <class T> const T* getKey_constPtr(const string& strKey);
	   
	   /**
	    * Template function to get a constant 2D pointer pointing to the data of a vector-type entry in the data container
	    *
	    * @param[in] strKey The key to be searched in the data container
	    * @param[in] iRows Number of rows
	    * @param[in] iColumns Number of columns
	    * @return    a constant 2D pointer if the key is found. Otherwise, an exception is thrown.
	    */
	   template <class T> const T** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns);

      /**
       * Template function to get a constant 3D pointer pointing to the data of a vector-type entry in the data container
       *
	    * @param[in] strKey The key to be searched in the data container
       * @param[in] iRows Number of rows
       * @param[in] iColumns Number of columns
       * @param[in] iPages Number of pages
       * @return    a constant 3D pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> const T*** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages);

      /**
       * Template function to get a constant 4D pointer pointing to the data of a vector-type entry in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @param[in] iRows Number of rows
       * @param[in] iColumns Number of columns
       * @param[in] iPages Number of pages
       * @param[in] iSects Number of sections
       * @return    a constant 4D pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> const T**** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages,const int& iSects);

	   //----------read and write key content in the map
      /**
       * Template function to get a reference to an entry of the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    a reference to an entry if the key is found. Otherwise, an exception is thrown.
       */
      template <class T> T& getKey_Ref(const string& strKey);

      /**
       * Template function to get the value to an entry of the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    Value of an entry if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> T getKey_Val(const string& strKey);
	   
      /**
       * Template function to get a pointer pointing to the data of a vector-type entry in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    a pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> T* getKey_Ptr(const string& strKey);

      /**
       * Template function to get a 3D pointer pointing to the data of a vector-type entry in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @param[in] iRows Number of rows
       * @param[in] iColumns Number of columns
       * @param[in] iPages Number of pages
       * @return    a 3D pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> T*** getKey_Ptr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages);

      /**
       * Function to show contents in the data container
       *
       * @param[in] strMapFile The file name that the data container to be written into
       */
	   virtual void showMap(const string& strMapFile) = 0;
	   
      /**
       * Function to reset data container by the values obtained from FORTRAN program
       *
       * @param[in] strF95File The file name storing values obtained from FORTRAN program
       */
	   virtual void reset(const string& strF95File) = 0;

	   virtual void sync_pre_space() = 0;
	   virtual void sync_pre_solver() = 0;
	   virtual void sync_pre_energy() = 0;
           virtual void sync_pre_site() = 0;
};
#else
class IDataContainer
{
   protected:   

      /**
       * data map containing variables to be read/modified among many other classes. Each entry contains a
       * key (string type) and associated value (boost::any type)
       */
      DataMap myData;
	               	                  
      /**
       * member function to set above data map. must be virtual since only the derived class knows the
       * contents to be written in the data container.
       */
	   virtual void setMap() = 0;

	public:

	   /**
	    * constructor
	    */
	   IDataContainer()
	   {
          #ifdef DEBUG_OBJECT
	      cout << endl;
	      cout << "****************************************************************\n";
	      cout << "*               IDataContainer is constructed                  *\n";
	      cout << "****************************************************************\n";
          #endif
	   };
	   
      /**
       * virtual destructor allowing an instance of a derived class can be deleted properly
       * through a pointer to this base class
       */
	   virtual ~IDataContainer()
	   {
	      DataMap().swap(myData);

          #ifdef DEBUG_OBJECT
          cout << endl;
          cout << "****************************************************************\n";
          cout << "*                IDataContainer is destroyed                   *\n";
          cout << "****************************************************************\n";
          #endif
	   };

      /**
       * Function to check if the user-specified key exists or not in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    true if the key exists, false otherwise
       */
	   bool keyExists(const string &strKey);
         
	   //----------read-only key content in the map
	   /**
	    * Template function to get a constant reference to an entry of the data container
	    *
	    * @param[in] strKey The key to be searched in the data container
	    * @return    a constant reference to an entry if the key is found. Otherwise, an exception is thrown.
	    */
	   template <class T> const T& getKey_constRef(const string& strKey);
	   	   
	   /**
	    * Template function to get a constant pointer pointing to the data of a vector-type entry in the data container
	    *
	    * @param[in] strKey The key to be searched in the data container
	    * @return    a constant pointer if the key is found. Otherwise, an exception is thrown.
	    */
	   template <class T> const T* getKey_constPtr(const string& strKey);
	   
	   /**
	    * Template function to get a constant 2D pointer pointing to the data of a vector-type entry in the data container
	    *
	    * @param[in] strKey The key to be searched in the data container
	    * @param[in] iRows Number of rows
	    * @param[in] iColumns Number of columns
	    * @return    a constant 2D pointer if the key is found. Otherwise, an exception is thrown.
	    */
	   template <class T> const T** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns);

      /**
       * Template function to get a constant 3D pointer pointing to the data of a vector-type entry in the data container
       *
	    * @param[in] strKey The key to be searched in the data container
       * @param[in] iRows Number of rows
       * @param[in] iColumns Number of columns
       * @param[in] iPages Number of pages
       * @return    a constant 3D pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> const T*** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages);

      /**
       * Template function to get a constant 4D pointer pointing to the data of a vector-type entry in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @param[in] iRows Number of rows
       * @param[in] iColumns Number of columns
       * @param[in] iPages Number of pages
       * @param[in] iSects Number of sections
       * @return    a constant 4D pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> const T**** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages,const int& iSects);

	   //----------read and write key content in the map
      /**
       * Template function to get a reference to an entry of the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    a reference to an entry if the key is found. Otherwise, an exception is thrown.
       */
      template <class T> T& getKey_Ref(const string& strKey);

      /**
       * Template function to get the value to an entry of the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    Value of an entry if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> T getKey_Val(const string& strKey);
	   
      /**
       * Template function to get a pointer pointing to the data of a vector-type entry in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    a pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> T* getKey_Ptr(const string& strKey);

      /**
       * Template function to get a 3D pointer pointing to the data of a vector-type entry in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @param[in] iRows Number of rows
       * @param[in] iColumns Number of columns
       * @param[in] iPages Number of pages
       * @return    a 3D pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> T*** getKey_Ptr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages);

      /**
       * Function to show contents in the data container
       *
       * @param[in] strMapFile The file name that the data container to be written into
       */
	   virtual void showMap(const string& strMapFile) = 0;

      /**
       * Function to reset data container by the values obtained from FORTRAN program
       *
       * @param[in] strF95File The file name storing values obtained from FORTRAN program
       */
	   virtual void reset(const string& strF95File) = 0;
};
#endif

//-----------------------------------------------------------------------//
template <class T> const T& IDataContainer::getKey_constRef(const string& strKey)
{
	//cout << "Getting Key_constref " << strKey << endl;

	if (!keyExists(strKey)) throw CInexistentKey(strKey);

	const T * nRetConstPtr = any_cast<const T>(&myData[strKey]);

	const T& nConstRetRef = *nRetConstPtr;

	return nConstRetRef;
}

//-----------------------------------------------------------------------//
template <class T> const T * IDataContainer::getKey_constPtr(const string& strKey)
{
	if (!keyExists(strKey)) throw CInexistentKey(strKey);

	const T * nConstRetPtr = any_cast<const T>(&myData[strKey]);

	return nConstRetPtr;
}

//-----------------------------------------------------------------------//
template <class T> const T ** IDataContainer::getKey_constPtr(const string& strKey, const int& iRows, const int& iColumns)
{
	if (!keyExists(strKey)) throw CInexistentKey(strKey);

	vector<T> * nConstVectorPtr = any_cast< vector<T> >(&myData[strKey]);

	if (nConstVectorPtr->size() != (unsigned int)iRows*iColumns) throw CMismatchSize(strKey);

	const T * nConstDataPtr = nConstVectorPtr->data();

	const T** prg2D = new const T *[iRows];

	for (int i = 0; i < iRows; i++)
		prg2D[i] = &nConstDataPtr[i*iColumns];

	return prg2D;
}

//-----------------------------------------------------------------------//
template <class T> const T *** IDataContainer::getKey_constPtr(const string& strKey, const int& iRows, const int& iColumns, const int& iPages)
{
	if (!keyExists(strKey)) throw CInexistentKey(strKey);

	vector<T> * nConstVectorPtr = any_cast< vector<T> >(&myData[strKey]);

	if (nConstVectorPtr->size() != (unsigned int)iRows*iColumns*iPages) throw CMismatchSize(strKey);

	T * nConstDataPtr = nConstVectorPtr->data();

	const T *** prg3D = new const T **[iRows];

	for (int i = 0; i < iRows; i++)
	{
		prg3D[i] = new const T *[iColumns];

		for (int j = 0; j < iColumns; j++)
			prg3D[i][j] = &nConstDataPtr[i*iColumns*iPages + j*iPages];
	}

	return prg3D;
}

//-----------------------------------------------------------------------//
template <class T> T *** IDataContainer::getKey_Ptr(const string& strKey, const int& iRows, const int& iColumns, const int& iPages)
{
	if (!keyExists(strKey)) throw CInexistentKey(strKey);

	vector<T> * nConstVectorPtr = any_cast< vector<T> >(&myData[strKey]);

	if (nConstVectorPtr->size() != (unsigned int)iRows*iColumns*iPages) 
	{
		throw CMismatchSize(strKey);
	}

	T * nConstDataPtr = nConstVectorPtr->data();

	T *** prg3D = new T **[iRows];
	for (int i = 0; i < iRows; i++)
	{
		prg3D[i] = new T *[iColumns];

		for (int j = 0; j < iColumns; j++)
			prg3D[i][j] = &nConstDataPtr[i*iColumns*iPages + j*iPages];
	}

	return prg3D;
}

//-----------------------------------------------------------------------//
template <class T> const T **** IDataContainer::getKey_constPtr(const string& strKey, const int& iRows, const int& iColumns, const int& iPages, const int& iSects)
{
	if (!keyExists(strKey)) throw CInexistentKey(strKey);

	vector<T> * nConstVectorPtr = any_cast< vector<T> >(&myData[strKey]);

	if (nConstVectorPtr->size() != (unsigned int)iRows*iColumns*iPages*iSects) throw CMismatchSize(strKey);

	const T * nConstDataPtr = nConstVectorPtr->data();

	const T **** prg4D = new const T ***[iRows];

	for (int i = 0; i < iRows; i++)
	{
		prg4D[i] = new const T **[iColumns];

		for (int j = 0; j < iColumns; j++)
		{
			prg4D[i][j] = new const T *[iPages];

			for (int k = 0; k < iPages; k++)

				prg4D[i][j][k] = &nConstDataPtr[i*iColumns*iPages*iSects + j*iPages*iSects + k*iSects];
		}
	}

	return prg4D;
}

//-----------------------------------------------------------------------//
template <class T> T& IDataContainer::getKey_Ref(const string& strKey)
{
	//cout << "Getting Key_ref " << strKey << endl;

	if (!keyExists(strKey)) throw CInexistentKey(strKey);

	T * nRetPtr = any_cast<T>(&myData[strKey]);

	T& nRetRef = *nRetPtr;

	return nRetRef;
}

//-----------------------------------------------------------------------//
template <class T> T IDataContainer::getKey_Val(const string& strKey)
{
	if (!keyExists(strKey)) throw CInexistentKey(strKey);

	T nRetVal = any_cast<T>(myData[strKey]);

	return nRetVal;
}

//-----------------------------------------------------------------------//
template <class T> T * IDataContainer::getKey_Ptr(const string& strKey)
{
	if (!keyExists(strKey)) throw CInexistentKey(strKey);

	T * nRetPtr = any_cast<T>(&myData[strKey]);

	return nRetPtr;
}


#endif // IDATACONTAINER_H_
