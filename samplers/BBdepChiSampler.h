#ifndef __BBDEPCHISAMPLER_H__
#define __BBDEPCHISAMPLER_H__

#include <map>
#include <vector>
using namespace std;

typedef pair<int,int> PII;
typedef struct { float chi1,chi2,chi3,chi4; } CHI;

class BBdepChiSampler {
public:
    BBdepChiSampler(const char *filename);
    virtual int sample(float phi, float psi, float* sample);
    virtual int sample(float phi, float psi, float* sample, int sampleIndex);

    BBdepChiSampler(); // for sake of inheritance
    virtual ~BBdepChiSampler();
private:
    void findOccupiedNbr(float phi, float psi, PII & i_phipsi);
    map<PII, vector<CHI> > psChis;
    map<PII, vector<float> > cumProbs;
    map<PII, PII> closestNbr;
};

#endif //__BBDEPCHISAMPLER_H__
