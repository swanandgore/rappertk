#ifndef __INFORMEDPEPTIDEBUILDER_H__
#define __INFORMEDPEPTIDEBUILDER_H__

#include "BasicPepBuilder.h"

class InformedPhipsiSampler;
class PhipsiSampler;
class OmegaSampler;

class InformedPeptideBuilder : public BasicPepBuilder {
public:
    InformedPeptideBuilder(vector<int>& ipInds, vector<int>& opInds,
	    Constants* con, const char *resn, const char* desc, PhipsiSampler *pss, OmegaSampler *os, InformedPhipsiSampler *ipps, bool fwdbkwd);
    void sample(VVF & pts, float& phi, float& psi, float& omega);
    InformedPeptideBuilder makeCopy();
    bool build(VVF & pts);
private:
    InformedPhipsiSampler *ipps;
    PhipsiSampler *pss;
    OmegaSampler *os;

    float r, a, t;
    set<long long int> sessionRATs;
    bool checkAndAddIfNew(float, float, float);
    void deleteSessionInfo();
};


#endif // __INFORMEDPEPTIDEBUILDER_H__
