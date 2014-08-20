#ifndef __CHIBUILDER_H__
#define __CHIBUILDER_H__

#include "Builder.h"
#include <set>
using std::set;

class BBdepChiSampler;

class ChiBuilder : public Builder {
public :
    ChiBuilder(vector<int>& ipInds, vector<int>& opInds,
	    Constants *consts, const char* desc, const char *resname, BBdepChiSampler *cs);
    bool build(VVF & pts);
    bool buildSample(VVF & pts, int sampleIndex);
    bool build(VVF & pts, float *sample);
    ChiBuilder makeCopy();

protected :
    void buildSC_VAL(VVF & pts, float *sample);
    void buildSC_ARG(VVF & pts, float *sample);
    void buildSC_LYS(VVF & pts, float *sample);
    void buildSC_GLU(VVF & pts, float *chis);
    void buildSC_ILE(VVF & pts, float *chis);
    void buildSC_ASN(VVF & pts, float *chis);
    void buildSC_LEU(VVF & pts, float *chis);
    void buildSC_HIS(VVF & pts, float *chis);
    void buildSC_THR(VVF & pts, float *chis);
    void buildSC_MET(VVF & pts, float *chis);
    void buildSC_GLN(VVF & pts, float *chis);
    void buildSC_ASP(VVF & pts, float *chis);
    void buildSC_PRO(VVF & pts, float *chis);
    void buildSC_SER(VVF & pts, float *chis);
    void buildSC_PHE(VVF & pts, float *chis);
    void buildSC_TYR(VVF & pts, float *chis);
    void buildSC_CYS(VVF & pts, float *chis);
    void buildSC_TRP(VVF & pts, float *chis);

private :
    BBdepChiSampler *chisampler;
    string resname;

    bool checkAndAddIfNew(int ci);
    set<int> sessionChis;
    void deleteSessionInfo();
};

#endif //__CHIBUILDER_H__
