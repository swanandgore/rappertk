#ifndef __SPHPOSRESTR_H__
#define __SPHPOSRESTR_H__

#include "Restraint.h"

#include <vector>
using namespace std;

class SphPosRestr : public Restraint {
public:    
    SphPosRestr(vector<int>& pis, const char *desc, vector<float>& pos, float rad0, float rad1);
    bool satisfied(VVF & pts);
    void describe();
private:
    float radMin, radMax, radMinsq, radMaxsq; // spherical shell. both radMin, radMax shd be positive
    vector<float> pos;
};
#endif // __SPHPOSRESTR_H__
