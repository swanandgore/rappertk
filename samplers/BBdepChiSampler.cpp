#include <assert.h>
#include <math.h>
#include "BBdepChiSampler.h"
#include "misc/RanGen.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
#include <iostream>
using namespace std;

void cumulate(vector<float> & vf) {
    if(vf.size() == 0) return;

    //for(unsigned int i=0; i < vf.size(); i++) cout << vf[i] << endl; cout << endl;
    float sum = 0.;
    for(unsigned int i=0; i < vf.size(); i++) sum += vf[i];

    for(unsigned int i=vf.size()-1; i >= 1; i--) vf[i] = vf[i-1];
    vf[0] = 0.;

    for(unsigned int i=1; i < vf.size(); i++) vf[i] = vf[i] + vf[i-1];

    for(unsigned int i=0; i < vf.size(); i++) vf[i] = vf[i]/sum;
    //for(unsigned int i=0; i < vf.size(); i++) cout << vf[i] << endl; cout << endl;
}

BBdepChiSampler::BBdepChiSampler() {}
BBdepChiSampler::~BBdepChiSampler() {}

BBdepChiSampler::BBdepChiSampler(const char *filename) {
    float prob; PII phipsi; CHI chi; char line[500];
    FILE *fp = fopen(filename, "r");
    while(fgets(line, 500, fp) != NULL) {
        sscanf(line, "%d %d %f %f %f %f %f", &phipsi.first, &phipsi.second, &prob, &chi.chi1, &chi.chi2, &chi.chi3, &chi.chi4);
        if(prob <= 1e-10) continue;
        assert(phipsi.first % 10 == 0); assert(phipsi.second % 10 == 0);
        phipsi.first /= 10; phipsi.second /= 10;
        if(phipsi.first == 18) phipsi.first = -18;
        if(phipsi.second == 18) phipsi.second = -18;
        if(cumProbs.find(phipsi) == cumProbs.end()) cumProbs[phipsi] = vector<float>();
        if(psChis.find(phipsi) == psChis.end()) psChis[phipsi] = vector<CHI>();
        cumProbs[phipsi].push_back(prob);
        psChis[phipsi].push_back(chi);
    }
    // for each phi-psi, arrange probabilities in descending order
    for(map<PII,vector<float> >::iterator it = cumProbs.begin(); it != cumProbs.end(); ++it) {
        vector<float> & p = it->second;
        vector<CHI> & c = psChis[it->first];
        float tp; CHI tc;
        for(int i=0; i < p.size(); i++)
            for(int k=i+1; k < p.size(); k++)
                if(p[i] < p[k]) {
                    tp = p[i]; p[i] = p[k]; p[k] = tp;
                    tc = c[i]; c[i] = c[k]; c[k] = tc;
                }
    }
    // for each phipsi, cumulate the probabilities
    for(map<PII,vector<CHI> >::iterator it=psChis.begin(); it != psChis.end(); it++) {
        phipsi = it->first;
        assert(psChis[phipsi].size() == cumProbs[phipsi].size());
        cumulate(cumProbs[phipsi]);
    }
    fclose(fp);
}

const char *spiral = "\
ul\
ddrr\
uuulll\
ddddrrrr\
uuuuulllll\
ddddddrrrrrr\
uuuuuuulllllll\
ddddddddrrrrrrrr\
uuuuuuuuulllllllll\
ddddddddddrrrrrrrrrr\
";

void BBdepChiSampler::findOccupiedNbr(float phi, float psi, PII & i_phipsi) {
    if(closestNbr.find(i_phipsi) != closestNbr.end()) {
        i_phipsi = closestNbr[i_phipsi]; return;
    }
    PII givenPhipsi = i_phipsi;
    // find phipsi closest to given one
    //cout << phi <<" phipsi " << psi << endl;
    assert(-180. <= phi and phi <= 180.); assert(-180. <= psi and psi <= 180.);
    i_phipsi.first = (int)(round(phi/10.)); i_phipsi.second = (int)(round(psi/10.));
    // spiral-walk till an occupied nbr is reached
    for(int step=0; step < strlen(spiral); step++) {
        if(i_phipsi.first == 18) i_phipsi.first = -18;
        else if(i_phipsi.first == -19) i_phipsi.first = 17;
        if(i_phipsi.second == 18) i_phipsi.second = -18;
        else if(i_phipsi.second == -19) i_phipsi.second = 17;

        if(psChis.find(i_phipsi) != psChis.end()) break;

        if(     spiral[step] == 'u') i_phipsi.second++;
        else if(spiral[step] == 'd') i_phipsi.second--;
        else if(spiral[step] == 'l') i_phipsi.first--;
        else if(spiral[step] == 'r') i_phipsi.first++;
    }
    assert(psChis.find(i_phipsi) != psChis.end());
    closestNbr[givenPhipsi] = i_phipsi;
}

int BBdepChiSampler::sample(float phi, float psi, float *sample, int sampleIndex) {
    PII i_phipsi;
    assert(-180. <= phi and phi <= 180.); assert(-180. <= psi and psi <= 180.);
    i_phipsi.first = (int) round(phi / 10.);
    i_phipsi.second = (int) round(psi / 10.);
    if(i_phipsi.first == 18) i_phipsi.first = -18;
    if(i_phipsi.second == 18) i_phipsi.second = -18;
    if(psChis.find(i_phipsi) == psChis.end()) {
        sample[0] = -9999; sample[1] = -9999; sample[2] = -9999; sample[3] = -9999;
        return -999;
    }

    int size = cumProbs[i_phipsi].size();
    if(sampleIndex >= size) return -999;
    sample[0] = psChis[i_phipsi][sampleIndex].chi1;
    sample[1] = psChis[i_phipsi][sampleIndex].chi2;
    sample[2] = psChis[i_phipsi][sampleIndex].chi3;
    sample[3] = psChis[i_phipsi][sampleIndex].chi4;
    return sampleIndex;
}

int BBdepChiSampler::sample(float phi, float psi, float *sample) {
    PII i_phipsi;
    //findOccupiedNbr(phi, psi, i_phipsi);
    //assert(psChis.find(i_phipsi) != psChis.end());
    assert(-180. <= phi and phi <= 180.); assert(-180. <= psi and psi <= 180.);
    i_phipsi.first = (int) round(phi / 10.);
    i_phipsi.second = (int) round(psi / 10.);
    if(i_phipsi.first == 18) i_phipsi.first = -18;
    if(i_phipsi.second == 18) i_phipsi.second = -18;
    if(psChis.find(i_phipsi) == psChis.end()) {
        sample[0] = -9999; sample[1] = -9999; sample[2] = -9999; sample[3] = -9999;
        return -999;
    }

    float rnum = ran01();
    //cout << rnum << endl;
    bool found = false;
    int size = cumProbs[i_phipsi].size();
    for(int i=0; i < size; i++) {
        //cout << cumProbs[i_phipsi][i] << endl;
        if( cumProbs[i_phipsi][i] <= rnum )
            if ( i+1 == size || rnum <= cumProbs[i_phipsi][i+1] ) {
                sample[0] = psChis[i_phipsi][i].chi1;
                sample[1] = psChis[i_phipsi][i].chi2;
                sample[2] = psChis[i_phipsi][i].chi3;
                sample[3] = psChis[i_phipsi][i].chi4;
                found = true;
                return i;
                //return psChis[i_phipsi].size()*((18+i_phipsi.first)*36 + 18+i_phipsi.second) + i;
            }
    }
    assert(found == true);
}
