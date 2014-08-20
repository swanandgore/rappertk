#ifndef __RATRESTRAINT_H__
#define __RATRESTRAINT_H__

#include <string>
using namespace std;

#include "Restraint.h"
class RATdata;

// this class checks whether for a consecutive C,N,CA triplet, there are phi-psi-omega
// triplets to make next CA within a certain sphere. useful for restraint propagation as well as debugging.
class RATrestraint : public Restraint {
public:
    RATrestraint(vector<int> & pis, const char *desc, vector<float> & p, float r, const char *resname, RATdata *ratdata);
    bool satisfied(VVF & pts);
private:
    vector<float> cen;
    float rad;
    string resn;
    RATdata * ratdata;
};

#endif // __RATRESTRAINT_H__
