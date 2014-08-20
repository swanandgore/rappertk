#ifndef __BASICPEPBUILDER_H__
#define __BASICPEPBUILDER_H__

#include "Builder.h"
#include <math.h>
#include <vector>
#include <map>
#include <set>
using namespace std;

// this builder shd do everything expected of peptide-builder,
// except the task of sampling of phi, psi, omega

class BasicPepBuilder : public Builder {
public :
    BasicPepBuilder(vector<int>& ipInds, vector<int>& opInds,
	    Constants* con, const char *resn, const char* desc, bool fwdbkwd);
    bool build(VVF & pts);
    virtual void sample(VVF & pts, float& phi, float& psi, float& omega);
    bool build1(VVF & pts, float phi, float psi, float omega);
    virtual ~BasicPepBuilder() {}
protected :
    float CA_C, N_CA_C, C_N, CA_C_N, N_CA, C_N_CA, C_O, CA_C_O, N_C_O;
    float CA_CB, N_CA_CB, C_N_CA_CB;

    bool fwd;
    bool buildBkwd(VVF & pts, float phi, float psi, float omega);
    bool buildFwd(VVF & pts, float phi, float psi, float omega);

    virtual bool checkAndAddIfNew(float phi, float psi, float omega);
    virtual void deleteSessionInfo();
    set<long long int> sessionPSO;
};

#endif // __BASICPEPBUILDER_H__
