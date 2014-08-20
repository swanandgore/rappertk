#ifndef __RAPSAMPLER_H__
#define __RAPSAMPLER_H__

class RapSampler {
public :
    RapSampler(const char *fn, const char *resn);
    void sample(float &phi, float &psi);
private :
    vector<int> phiseq, psiseq;
    int seqnum;
};

#endif // __RAPSAMPLER_H__
