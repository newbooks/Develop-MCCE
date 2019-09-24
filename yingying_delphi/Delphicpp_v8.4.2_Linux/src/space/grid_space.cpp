#include "space.h"

/**
 * int_coord(): covert 3 delphi_integer numbers to be a SGrid <delphi_integer> structure
 */
SGrid <delphi_integer> CDelphiSpace::int_coord(const delphi_integer& a, const delphi_integer& b, const delphi_integer& c)
{
    SGrid <delphi_integer>  d;
    d.nX=a;
    d.nY=b;
    d.nZ=c;

    return(d);
}

/**
 * coord(): covert 3 float numbers to be a SGrid <float> structure
 */
SGrid <delphi_real> CDelphiSpace::coord(const delphi_real& a,  const delphi_real& b,  const delphi_real& c)
{
    SGrid <delphi_real>  d;
    d.nX=a;
    d.nY=b;
    d.nZ=c;

    return(d);
}

/**
 * Float2Int(): 1. convert SGrid <float> structure to SGrid <int> structure; 2. convert float to int
 */
SGrid <int> CDelphiSpace::Float2Int( const SGrid <float>& a )
{
    SGrid <int> b;

    b.nX=int(a.nX);
    b.nY=int(a.nY);
    b.nZ=int(a.nZ);

    return(b);
}

int CDelphiSpace::Float2Int(const float& a)
{
    int b;
    b=int(a);

    return(b);
}

/**
 * Int2Float(): 1. convert SGrid <int> structure to SGrid <float> structure; 2. convert int to float
 */
SGrid <float> CDelphiSpace::Int2Float( const SGrid <int>& a )
{
    SGrid <float> b;

    b.nX=float(a.nX);
    b.nY=float(a.nY);
    b.nZ=float(a.nZ);

    return(b);
}

float CDelphiSpace::Int2Float(const int& a)
{
    float b;
    b=float(a);

    return(b);
}




