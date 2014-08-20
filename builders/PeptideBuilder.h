#ifndef __PEPTIDEBUILDER_H__
#define __PEPTIDEBUILDER_H__

#include "BasicPepBuilder.h"

class PhipsiSampler;
class OmegaSampler;

class PeptideBuilder : public BasicPepBuilder {
public :
    PeptideBuilder(vector<int>& ipInds, vector<int>& opInds,
	    Constants* con, const char *resn, const char* desc, PhipsiSampler* p, OmegaSampler *o, bool fwdbkwd);
    void sample(VVF & pts, float& phi, float& psi, float& omega);
    PeptideBuilder makeCopy();
private :
    PhipsiSampler* pps;
    OmegaSampler* os;
};

#endif //__PEPTIDEBUILDER_H__
