#ifndef __NANCHORBUILDER_H__
#define __NANCHORBUILDER_H__

#include "Builder.h"

class NanchorBuilder : public Builder {
public:
    NanchorBuilder(vector<int> & ipInds, vector<int> & opInds, Constants *con, const char *desc,
        vector<float>& CA0, float R0, vector<float>& CA1, float R1);
    bool build(VVF & pts);
    NanchorBuilder makeCopy();
private:
    vector<float> ca0, ca1;
    float r0, r1;
};

#endif // __NANCHORBUILDER_H__
