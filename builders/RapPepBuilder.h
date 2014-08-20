#ifndef __RAPBUILDER_H__
#define __RAPBUILDER_H__

#include "BasicPepBuilder.h"

class RapSampler;
class OmegaSampler;

class RapPepBuilder : public BasicPepBuilder {
public :
    RapPepBuilder(vector<int>& ipInds, vector<int>& opInds,
	    Constants* con, const char *resn, const char* desc, RapSampler* rs, OmegaSampler *o, bool fwdbkwd);
    void sample(VVF & pts, float& phi, float& psi, float& omega);
private :
    RapSampler* raps;
    OmegaSampler* os;
};

#endif // __RAPBUILDER_H__
