#ifndef __DIHEDRESTRAINT_H__
#define __DIHEDRESTRAINT_H__

#include <vector>
using namespace std;

#include "Restraint.h"

class DihedRestraint : public Restraint {
public :
    DihedRestraint(vector<int>& pis, const char *description, vector<float> expectedVals, vector<float> errorRanges);
    bool satisfied(VVF & pts);
    void describe();
private :
    vector<float> expected, error;
};

#endif //__DIHEDRESTRAINT_H__
