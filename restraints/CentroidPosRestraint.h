#ifndef __CENTROIDPOSRESTRAINT_H__
#define __CENTROIDPOSRESTRAINT_H__

#include <vector>
using namespace std;

#include "Restraint.h"

class CentroidPosRestraint : public Restraint {
public :
    CentroidPosRestraint(vector<int>& pis, const char *desc, vector<float>& c, float r);
    bool satisfied(VVF & pts);
    void describe();
private:
    vector<float> cen;
    float rad;
};


#endif // __CENTROIDPOSRESTRAINT_H__
