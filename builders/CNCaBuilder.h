#ifndef __CNCABUILDER_H__
#define __CNCABUILDER_H__

#include "Builder.h"

class Constants;

class CNCaBuilder : public Builder {
public:
    CNCaBuilder(vector<int>& ipInds, vector<int>& opInds,
        Constants *con, const char* desc, vector<float> & capos, float capostol);
    bool build(VVF & pts);
protected:
    vector<float> CApos;
    float CAposTol;
};

#endif // __CNCABUILDER_H__
