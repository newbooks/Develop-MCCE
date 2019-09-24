#ifndef SPACE_EXCEPTIONS_H_
#define SPACE_EXCEPTIONS_H_

#include "../interface/exceptions.h"

class CZeroChargeRadius : public CWarning 
{
public:
    CZeroChargeRadius(const int &cnt)
    {
        cwarn << cnt << " CHARGED ATOMS WERE FOUND TO HAVE ZERO RADIUS." << endl;

    }
};


#endif /* SPACE_EXCEPTIONS_H_ */
