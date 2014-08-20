#ifndef __INFORMEDPHIPSISAMPLER_H__
#define __INFORMEDPHIPSISAMPLER_H__

#include <vector>
#include <string>
using namespace std;

class InformedPhipsiSampler {
public:
    InformedPhipsiSampler(vector<float> & target, float targetTol, const char* resn, const char *ratdatapath, bool fwdbwd);
    void sample(vector<float>& c, vector<float>& n, vector<float>& ca, float & phi, float & psi, float & omega,
                float &r, float &a, float &t);
protected:
    string resn, RATdataPath;
    vector<float> target; float targetTol;
    bool fwd;
};

#endif // __INFORMEDPHIPSISAMPLER_H__
