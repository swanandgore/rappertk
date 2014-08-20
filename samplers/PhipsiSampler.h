#ifndef __PHIPSISAMPLER_H__
#define __PHIPSISAMPLER_H__

#include<vector>
#include<map>
using namespace std;

#include "Sampler.h"

class PhipsiSampler : public Sampler {
public :
    PhipsiSampler(const char *filename, int coilAlphaBeta);
    void readPhipsiMap(const char *filename);
    float findProb(float phi, float psi);
    int sample(int& phi, int& psi);
    void printSample();
    int maxUniqueSamples();
protected :
    static int magfactor;
    vector<int> phis, psis;
    vector<int> toPhipsiIndex;
    int coilORalphaORbeta;
};

#endif //__PHIPSISAMPLER_H__
