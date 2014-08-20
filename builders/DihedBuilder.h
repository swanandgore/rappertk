#ifndef __DIHEDBUILDER_H__
#define __DIHEDBUILDER_H__

#include "Builder.h"

class DihedBuilder : public Builder {
public:
    DihedBuilder(vector<int>& ipInds, vector<int>& opInds, const char *desc, float length, float angle, float dihed);
    bool build(VVF & pts);
    DihedBuilder makeCopy();
private:
    float len, ang, dih;
};

#endif //__DIHEDBUILDER_H__
