#ifndef __FRAGSAMPLER_H__
#define __FRAGSAMPLER_H__

#include "Sampler.h"

#include <iostream>
#include <vector>
using namespace std;

class FragSampler : public Sampler {
public :
    FragSampler(int numfrags, vector<vector<float> > & pts);
    int sample(vector<vector<float> > & pts);
    int maxUniqueSamples();
protected :
    vector<vector<vector<float> > > samples; // assume equiprobable samples for the time being
};

#endif//__FRAGSAMPLER_H__
