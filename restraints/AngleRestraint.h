#ifndef __ANGLERESTRAINT_H__
#define __ANGLERESTRAINT_H__

#include <vector>
using namespace std;

#include "Restraint.h"

class AngleRestraint : public Restraint {
public :
    AngleRestraint(vector<int>& pis, const char *description, float min, float max);
    bool satisfied(VVF & pts);
    void describe();
private :
    float minAng, maxAng;
};

#endif //__ANGLERESTRAINT_H__
