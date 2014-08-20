#ifndef __DISTANCERESTRAINT_H__
#define __DISTANCERESTRAINT_H__

#include <vector>
using namespace std;

#include "Restraint.h"

class DistanceRestraint : public Restraint {
public :
    DistanceRestraint(vector<int>& pis, const char *description, float mindist, float maxdist);
    bool satisfied(VVF & pts);
    void describe();
private :
    float minsq, maxsq;
};

#endif //__DISTANCERESTRAINT_H__
