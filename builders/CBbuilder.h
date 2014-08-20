#ifndef __CBBUILDER_H__
#define __CBBUILDER_H__

#include "Builder.h"

class CBbuilder : public Builder {
public:
    CBbuilder(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char *desc, const char *resn);
    bool build(VVF & pts);
    bool buildSample(VVF & pts, int sampleIndex);
    CBbuilder makeCopy();
private :
    float CA_CB, C_CA_CB, N_CA_CB, N_C_CA_CB, C_N_CA_CB;
};

#endif //__CBBUILDER_H__
