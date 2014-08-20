#ifndef __NASUITESAMPLER_H__
#define __NASUITESAMPLER_H__

#include <map>
#include <vector>
#include <string>
using namespace std;

class SuiteVal {
    public:
        SuiteVal(int deltaPrev, int epsilon, int zeta, int alpha, int beta, int gamma, int deltaNext);
        SuiteVal();
        int diff(SuiteVal & sv, bool verbose=false);
        void describe(ostream & os);

        int d0, e, z, a, b, g, d1;
};

class NAsuiteSampler {
public:
    NAsuiteSampler(const char *filename, bool randomize);
    void readRNAsuite(const char *filename);
    void sample(int puckerPrev, int puckerNext, SuiteVal & sv);
    SuiteVal closestSuiteVal(SuiteVal & sv);
    string closestSuiteKey(SuiteVal & sv);
protected:
    map<string,SuiteVal> keySuite;
    
    vector<float> probs22; vector<string> keys22;
    vector<float> probs23; vector<string> keys23;
    vector<float> probs32; vector<string> keys32;
    vector<float> probs33; vector<string> keys33;

    bool randomize;
};


#endif // __NASUITESAMPLER_H__
