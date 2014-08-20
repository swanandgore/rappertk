#ifndef __SAMPLE_H__
#define __SAMPLE_H__

class Sampler {
public :
    // gives max unique samples that a sampler can provide
    virtual int maxUniqueSamples() = 0;
    virtual ~Sampler() {}
};

#endif//__SAMPLE_H__
