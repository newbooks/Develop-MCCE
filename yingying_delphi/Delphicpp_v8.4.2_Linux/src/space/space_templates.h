#ifndef SPACE_TEMPLATES_H
#define SPACE_TEMPLATES_H

#include "../interface/environment.h"

template <class T> void get_pt2d(T ** & pt2d, const delphi_integer& length, const delphi_integer& width)
{
    delphi_integer i;

    pt2d= new T *[length];
    for (i=0; i<length; i++)
        pt2d[i] = new T [width]();
}

template <class T> void get_pt3d(T *** & pt3d, const delphi_integer& length, const delphi_integer& width, const delphi_integer& hight)
{
    delphi_integer i,j;

    pt3d= new T **[length];
    for (i=0; i<length; i++)
    {
        pt3d[i] = new T* [width];
        for (j=0; j<width; j++)
            pt3d[i][j]= new T[hight]();
    }
}

template <class T> void get_pt4d(T **** & pt4d, const delphi_integer& length, const delphi_integer& width, const delphi_integer& hight, const delphi_integer& page)
{
    delphi_integer i,j,k;

    pt4d= new T ***[length];
    for (i=0; i<length; i++)
    {
        pt4d[i] = new T** [width];
        for (j=0; j<width; j++)
        {
            pt4d[i][j]= new T* [hight];

            for(k=0; k<hight; k++)
                pt4d[i][j][k]= new T [page]();
        }
    }
}

template <class T> void free_pt2d(T ** & pt2d, const delphi_integer& length, const delphi_integer& width)
{
    delphi_integer i;

    pt2d= new T *[length];
    for (i=0; i<length; i++)
        delete [] pt2d[i];

    delete [] pt2d;
    pt2d=NULL;
}

template <class T> void free_pt3d(T *** & pt3d, const delphi_integer& length, const delphi_integer& width, const delphi_integer& hight)
{
    delphi_integer i,j;

    for (i=0; i<length; i++)
    {
        for (j=0; j<width; j++)
            delete [] pt3d[i][j];

        delete [] pt3d[i];
    }
    delete [] pt3d;
    pt3d=NULL;
}

template <class T> void free_pt4d(T **** & pt4d, const delphi_integer& length, const delphi_integer& width, const delphi_integer& hight, const delphi_integer& page)
{
    delphi_integer i,j,k;

    for (i=0; i<length; i++)
    {
        for (j=0; j<width; j++)
        {
            for(k=0; k<hight; k++)
                delete [] pt4d[i][j][k];

            delete []  pt4d[i][j];
        }
        delete [] pt4d[i];
    }
    delete [] pt4d;
    pt4d=NULL;
}

#endif // SPACE_TEMPLATES_H
