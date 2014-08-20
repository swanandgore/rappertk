#include <assert.h>

#include "NAsuiteSampler.h"
#include "misc/RanGen.h"
#include "misc/verbosity.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
#include<stdlib.h>

#include<iostream>
using namespace std;

SuiteVal::SuiteVal() : d0(-999), e(-999), z(-999), a(-999), b(-999), g(-999), d1(-999)
{}

SuiteVal::SuiteVal(int deltaPrev, int epsilon, int zeta, int alpha, int beta, int gamma, int deltaNext)
    : d0(deltaPrev), e(epsilon), z(zeta), a(alpha), b(beta), g(gamma), d1(deltaNext)
{}

float dihedDiff(float theta1, float theta2) {
    while(theta1 <   0) theta1 += 360;
    while(theta1 > 360) theta1 -= 360;
    while(theta2 <   0) theta2 += 360;
    while(theta2 > 360) theta2 -= 360;
    float diff = fabs(theta1-theta2);
    if(diff > 180) diff = 360-diff;
    return diff;
}

void SuiteVal::describe(ostream & os) {
    os << "\t\t" << d0 <<"\t"<< e <<"\t"<< z <<"\t"<< a <<"\t"<< b <<"\t"<< g <<"\t"<< d1 << endl;
}

int SuiteVal::diff(SuiteVal & sv, bool verbo) {
    int d = 0;
    int d0diff = dihedDiff(d0 , sv.d0);  d += d0diff;
    int d1diff = dihedDiff(d1 , sv.d1);  d += d1diff;
    int ediff = dihedDiff(e , sv.e);    d += ediff;
    int zdiff = dihedDiff(z , sv.z);    d += zdiff;
    int adiff = dihedDiff(a , sv.a);    d += adiff;
    int bdiff = dihedDiff(b , sv.b);    d += bdiff;
    int gdiff = dihedDiff(g , sv.g);    d += gdiff;
    if(verbo)
	    cout << "diffs "
	    << d0diff << " " << ediff << " " << zdiff << " " << adiff << " "
        << bdiff << " " << gdiff << " " << d1diff << " " << endl;
    return d;
}

void NAsuiteSampler::readRNAsuite(const char *filename) {
    FILE *fp = fopen(filename, "r");
    int d0, e, z, a, b, g, d1;
    char key[10], line[500];
    while(fgets(line, 500, fp) != NULL) {
        if(line[0] == '#') continue;
        sscanf(line, "%s %d %d %d %d %d %d %d", key, &d0, &e, &z, &a, &b, &g, &d1);
        keySuite[key] = SuiteVal(d0, e, z, a, b, g, d1);
    }
}

void fillProbs(vector<string>& keys, vector<float>& probs) {
    probs.resize( keys.size(), 0 );
    for(int i=1; i < keys.size(); i++) probs[i] = (i+0.)/keys.size();
}

NAsuiteSampler::NAsuiteSampler(const char *filename, bool _randomize) {
    randomize = _randomize;

    readRNAsuite(filename);
    for(map<string,SuiteVal>::iterator it=keySuite.begin(); it != keySuite.end(); ++it) {
        const string & key = it->first;
        if(key[0]=='2' && key[6]=='2') keys22.push_back(key);
        else if(key[0]=='2' && key[6]=='3') keys23.push_back(key);
        else if(key[0]=='3' && key[6]=='2') keys32.push_back(key);
        else if(key[0]=='3' && key[6]=='3') keys33.push_back(key);
        else assert(0);
    }
    fillProbs(keys22, probs22);
    fillProbs(keys23, probs23);
    fillProbs(keys32, probs32);
    fillProbs(keys33, probs33);
    //    if(verbose(7)) cout << "KEY22 23 32 33 " << keys22.size() <<" "<< keys23.size() <<" "<< keys32.size() <<" "<< keys33.size() << endl;
}

void NAsuiteSampler::sample(int puckerPrev, int puckerNext, SuiteVal & sv) {
    vector<float> & probs = probs22;
    vector<string> & keys = keys22;
    if(puckerPrev == 2 && puckerNext == 2) { probs = probs22; keys = keys22; }
    else if(puckerPrev == 2 && puckerNext == 3) { probs = probs23; keys = keys23; }
    else if(puckerPrev == 3 && puckerNext == 2) { probs = probs32; keys = keys32; }
    else if(puckerPrev == 3 && puckerNext == 3) { probs = probs33; keys = keys33; }
    else assert(0);
    //    if(verbose(7)) cout << "PROBS SIZE " << probs.size() << endl;
    float dice = ran01();
    for(int i=0; i < probs.size(); i++) {
        if(probs[i] < dice) {
            if(i==probs.size()-1 || dice < probs[i+1]) {
                sv = SuiteVal(keySuite[ keys[i] ]);
                // change all values a bit by adding +- 20 noise
                if(randomize) {
                    sv.e += (ran01() * 20 - 10);
                    sv.z += (ran01() * 20 - 10);
                    sv.a += (ran01() * 20 - 10);
                    sv.b += (ran01() * 20 - 10);
                    sv.g += (ran01() * 20 - 10);
                }
                return;
            }
        }
    }
}

SuiteVal NAsuiteSampler::closestSuiteVal(SuiteVal & sv) { return keySuite[ closestSuiteKey(sv) ]; }

string NAsuiteSampler::closestSuiteKey(SuiteVal & sv) {
    int maxdiff = 10000, diff;
    string retkey;
    for(map<string,SuiteVal>::iterator it = keySuite.begin(); it != keySuite.end(); ++it) {
        SuiteVal & sv1 = it->second;
        diff = sv.diff(sv1);
        //cout << it->first << " : " << diff << endl;
        if(diff < maxdiff) { maxdiff = diff; retkey = it->first; }
    }
    cout << retkey << " ------- " << maxdiff << endl;
    keySuite[retkey].diff(sv, true);
    keySuite[retkey].describe(cout);
    return retkey;
}
