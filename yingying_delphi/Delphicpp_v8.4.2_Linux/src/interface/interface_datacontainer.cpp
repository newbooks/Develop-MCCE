/*
 * interface_datacontainer.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: chuan
 */

#include "interface_datacontainer.h"
#include "interface_datacontainer_impls.cpp"
#include <iostream>

//-----------------------------------------------------------------------//
bool IDataContainer::keyExists(const string& strKey)
{
   DataMap::iterator it = myData.find(strKey);
   
   if (myData.end() == it) return false;
   
   return true; 
}  

