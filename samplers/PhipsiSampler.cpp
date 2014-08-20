
#include <math.h>
#include <assert.h>
#include<cstring>
#include<stdio.h>
#include "misc/RanGen.h"
#include "PhipsiSampler.h"

#include<iostream>
using namespace std;

int PhipsiSampler::magfactor = 1000000;

void add2map(map<pair<int,int>, float> & themap, int phi, int psi, float prob) {
    while(phi >= 180) phi -= 360; while(phi < -180) phi += 360;
    while(psi >= 180) psi -= 360; while(psi < -180) psi += 360;

    if(themap.find(pair<int,int>(phi,psi)) == themap.end())
        themap[pair<int,int>(phi,psi)] = prob;
    else
        themap[pair<int,int>(phi,psi)] += prob;
}

// put a mask around each phi-psi state, such that it retains x% of its value,
// and 100-x % is distributed equally among its 8 nbrs. do this just once.
// assume that phi-psi are on 5-degree grid and vary between -180 to +180
void smudge(vector<int> & phis, vector<int> & psis, vector<float> & probs) {
    float fudgeFactor = 0.8;
    map< pair<int,int>, float > themap, newmap;
    for(int i=0; i < phis.size(); i++)
        themap[ pair<int,int>(phis[i],psis[i]) ] = probs[i];
    for(int i=0; i < phis.size(); i++) {
        int phi = phis[i], psi = psis[i];
        float prob = themap[ pair<int,int>(phi,psi) ];
        float newprob = prob * fudgeFactor, dprob = (1.-fudgeFactor)/8 * prob;
        for(int p = phi-5; p < phi+6; p+=5)
            for(int s = psi-5; s < psi+6; s+=5)
                if(p == phi and s == psi) add2map(newmap, p, s, newprob);
                else add2map(newmap, p, s, dprob);
    }
    phis.clear(); psis.clear(); probs.clear();
    cout << "SIZES " << themap.size() <<" "<< newmap.size() << endl;
    for(map<pair<int,int>,float>::iterator it = newmap.begin(); it != newmap.end(); ++it) {
        if(themap.find(it->first) == themap.end()) continue; // dont create new states! XXX
        phis.push_back(it->first.first);
        psis.push_back(it->first.second);
        probs.push_back(it->second);
    }
}

void PhipsiSampler::readPhipsiMap(const char *filename) {
    float phi, psi, prob;
    vector<float> cumProb;
    char line[500];
    FILE *fp = fopen(filename, "r");
    while(fgets(line, 500, fp) != NULL) {
        sscanf(line, "%f %f %f", &phi, &psi, &prob);
        if(prob <= 1e-10) continue;
        assert(5*round(phi/5.) == phi);
        assert(5*round(psi/5.) == psi);
        cumProb.push_back(prob); //XXX
        phis.push_back(int(phi));
        psis.push_back(int(psi));
    }

    //for(int i=0; i < phis.size(); i++) cout << "before PSPROP " << phis[i] <<" "<< psis[i] <<" "<< cumProb[i] << endl;

    //for(int smudging=0; smudging < 1; smudging++) smudge(phis, psis, cumProb);

    //for(int i=0; i < phis.size(); i++) cout << "after PSPROP " << phis[i] <<" "<< psis[i] <<" "<< cumProb[i] << endl;

    //for(int i=0; i < cumProb.size(); i++) cout << cumProb[i] <<" "; cout << endl;
    for(unsigned int i=1; i < cumProb.size(); i++) cumProb[i] = cumProb[i] + cumProb[i-1];
    //for(int i=0; i < cumProb.size(); i++) cout << cumProb[i] <<" "; cout << endl;
    float tot = cumProb[cumProb.size()-1];
    //cout << "cum num " << tot << endl;
    for(unsigned int i=0; i < cumProb.size(); i++) cumProb[i] = cumProb[i]/tot;
    //for(int i=0; i < cumProb.size(); i++) cout << cumProb[i] <<" "; cout << endl;

    int startProb = 0, endProb = 0;
    toPhipsiIndex.resize(magfactor);
    for(unsigned int i=0; i < cumProb.size(); i++) {
        //cout << i << endl;
        endProb = (int) round (cumProb[i] * magfactor);
        for(int k=startProb; k < endProb; k++) toPhipsiIndex[k] = i;
        startProb = endProb;
    }
    //for(int i=0; i < toPhipsiIndex.size(); i++) cout << "PS " << i <<" "<< toPhipsiIndex[i] << endl;
    //for(int i=0; i < toPhipsiIndex.size(); i++) cout << i <<" "<< toPhipsiIndex[i] <<" "<< phis[toPhipsiIndex[i]] <<" "<< psis[toPhipsiIndex[i]] << endl;
    //exit(0);
}

float PhipsiSampler::findProb(float phi, float psi) {
    phi = round(phi / 5.) * 5;
    psi = round(psi / 5.) * 5;
    while(phi >= 180) phi -= 360; while(phi < -180) phi += 360;
    while(psi >= 180) psi -= 360; while(psi < -180) psi += 360;
    int psindex = -100;
    for(int i=0; i < phis.size(); i++) {
        //cout << phi <<" "<< psi <<" "<< phis[i] <<" "<< psis[i] << endl;
        if(phi == phis[i] && psi == psis[i])
            { psindex = i; break; }
    }
    if(psindex < 0) return 0;
    int prob = 0;
    for(int i=0; i < toPhipsiIndex.size(); i++)
        if(psindex == toPhipsiIndex[i]) prob++;
    return (prob+0.) / toPhipsiIndex.size();
}

PhipsiSampler::PhipsiSampler(const char *filename, int coilAlphaBeta) {
    if(phis.size() == 0) readPhipsiMap(filename);
    coilORalphaORbeta = coilAlphaBeta;
    if(coilAlphaBeta != 0 && coilAlphaBeta != 1 && coilAlphaBeta != 2) {
        cout << "Secondary structure mask on PhipsiSampler must be 0, 1 or 2 not " << coilAlphaBeta << endl;
        exit(-1);
    }
}

void PhipsiSampler::printSample() {
    int phi, psi;
    sample(phi,psi);
    cout << phi <<" "<< psi << endl;
}

int PhipsiSampler::sample(int& phi, int& psi) {
    while(1) {
        int myrand = (int) floor ( magfactor * ran01() );
        int psindex = toPhipsiIndex[myrand];
        //cout << myrand <<":"<< psindex << endl;
        phi = phis[psindex];
        psi = psis[psindex];
        if(coilORalphaORbeta == 0) return psindex;
        else if(coilORalphaORbeta == 1) { // alpha
            if((-67 <= phi && phi <= -47)   && (-57 <= psi && psi <= -37)) return psindex;
        }
        else { // beta
            if((-149 <= phi && phi <= -109) && (103 <= psi && psi <= 145)) return psindex;
        }
    }
}

int PhipsiSampler::maxUniqueSamples() {
    return phis.size();
}

#ifdef DONT_DEFINE_ME
main() {
    PhipsiSampler pss("test.ps"); // ("../data/PhipsiWeightedProp/psALA");
    cout << "prob " << pss.findProb(31,31) << endl;
    for(int i=0; i < 1000; i++) {
        float a, b;
        pss.sample(a,b);
        cout << "sampled " << a << " " << b << endl;
    }
    //PhipsiSampler pss("../data/PhipsiWeightedProp/psALA");
}
#endif
