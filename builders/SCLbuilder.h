#ifndef __SCLBUILDER_H__
#define __SCLBUILDER_H__

#include "Builder.h"
#include <set>
using std::set;

class SCLsampler;

class SCLbuilder : public Builder {
public :
    SCLbuilder(vector<int>& ipInds, vector<int>& opInds,
	    Constants *consts, const char* desc, SCLsampler *cs);
    bool build(VVF & pts);
    bool buildSample(VVF & pts, int sampleIndex);
    SCLbuilder makeCopy();

private :
    void getPhipsi(float & phi, float & psi, VVF & pts);
    bool basicbuild(VVF *rotpts, VVF & pts);

    SCLsampler *sclsampler;
    vector<float> t1; vector<vector<float> > rot; vector<float> t2; // dont create again and again

    bool checkAndAddIfNew(int ci);
    set<int> sessionChis;
    void deleteSessionInfo();
};

#endif //__SCLBUILDER_H__
