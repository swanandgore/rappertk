#ifndef __PROTRESBUILDER_H__
#define __PROTRESBUILDER_H__

class PsoChis {
public:
    PsoChis() {}
    PsoChis(int p, int s, int o, int c1, int c2, int c3, int c4);
    bool operator()(const PsoChis & psoca, const PsoChis & psocb);
    int phi,psi,omega,chi1,chi2,chi3,chi4;
};

#include <set>
using namespace std;

#include "Builder.h"

class BasicPepBuilder;
class CBbuilder;
class ChiBuilder;
class PhipsiSampler;
class BBdepChiSampler;
class OmegaSampler;

class ProtResBuilder : public Builder {
public:
    ProtResBuilder() {}
    ProtResBuilder(vector<int>& ipInds, vector<int>& opInds, Constants *consts, const char* desc, const char *resname,
        PhipsiSampler *p, OmegaSampler *o, BBdepChiSampler *cs, vector<float>& CA, float CArad, bool fwdbkwd);
    ~ProtResBuilder();
    bool build(VVF & pts);
    bool checkAndAddIfNew(int, int, int);
    ProtResBuilder makeCopy();
private:
    PhipsiSampler *pps; 
    OmegaSampler *os;
    BBdepChiSampler *chisampler;
    vector<float> ca; float carad; // CA restraint

    BasicPepBuilder *pb;
    CBbuilder *cbb;
    ChiBuilder *chib;

    void deleteSessionInfo();
    set<long> sessionIds, psoBlacklist;
};

#endif // __PROTRESBUILDER_H__
