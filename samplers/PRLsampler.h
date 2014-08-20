#ifndef __PRLSAMPLE_H__
#define __PRLSAMPLE_H__

#include "BBdepChiSampler.h"

class PRLsampler : public BBdepChiSampler {
public :
    PRLsampler(const char *filename, char aa);
    int sample(float phi, float psi, float* sample);
    int sample(float phi, float psi, float* sample, int sampleIndex);
    void addSample(float prop, vector<float> chis);
    PRLsampler makeCopy();
    void describe();
public :
    map< char,vector<CHI> > chis; // map from state(Alpha,Beta,Other) to chi
    map< char,vector<int> > r2ind; // for each state, random number to chi-index
};

#endif // __PRLSAMPLE_H__
