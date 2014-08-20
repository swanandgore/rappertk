#ifndef __OMEGASAMPLER_H__
#define __OMEGASAMPLER_H__

#include "Sampler.h"

#include <string>
using std::string;


#include <map>
using std::map;
#include <string>
using std::string;

#include<iostream>; 
using namespace std ; 

class OmegaSampler : public Sampler {
public :
    OmegaSampler(float transprobability);
    int sample(int& omega);
    int maxUniqueSamples();
private :
    float transprop; // probablibility to be in trans, else cis
};

#endif //__OMEGASAMPLER_H__
