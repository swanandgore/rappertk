#include <assert.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "PRLsampler.h"
#include "misc/RanGen.h"
#include "misc/verbosity.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
// since min val in richardson.lib is 0.005
#define MAGFAC (1000)

PRLsampler::PRLsampler(const char *filename, char AA) {
    float prob;
    char aa, ss, line[500];
    CHI chi;
    int dontcare;

    map<char, vector<float> > ssChiProb;
    ssChiProb['A'] = vector<float>();
    ssChiProb['B'] = vector<float>();
    ssChiProb['O'] = vector<float>();

    FILE *fp = fopen(filename, "r");
    while(fgets(line, 500, fp) != NULL) {
        if(line[0] == '#') continue;
        sscanf(line, "%c %d %f %f %f %f %f %d %d %c", &aa, &dontcare, &chi.chi1, &chi.chi2, &chi.chi3, &chi.chi4, &prob, &dontcare, &dontcare, &ss);
        if(aa != AA) continue;
        if(prob <= 1e-10) continue;
        assert(ss=='A' || ss=='B' || ss=='O');
        chis[ss].push_back(chi);
        ssChiProb[ss].push_back(prob);
    }

    char *ssStates = "ABO";
    for(int ssi=0; ssi<3; ssi++) { // normalize and make jump-map
        ss = ssStates[ssi];
        float sum = 0.;
        for(int i=0; i < ssChiProb[ss].size(); i++) sum += ssChiProb[ss][i];
        for(int i=0; i < ssChiProb[ss].size(); i++) ssChiProb[ss][i] /= sum;
        int start = 0, end;
        for(int i=0; i < ssChiProb[ss].size(); i++) {
            end = (int) (floor( ssChiProb[ss][i] * MAGFAC )) + start;
            if(i == ssChiProb[ss].size()-1) end = MAGFAC;
            for(; start < end; start++) r2ind[ss].push_back(i);
        }
    }

    /*        if(verbose(7)) {
	      cout << "read PRL lib for residue " << AA << endl;
        for(int ssi=0; ssi < 3; ssi++) {
	ss = ssStates[ssi];
            cout << ss << " " << ssChiProb[ss].size() << " ";
            for(int i=0; i < ssChiProb[ss].size(); i++) cout << ssChiProb[ss][i] << " ";
             cout << endl;
	    }
	    }*/
}

int PRLsampler::sample(float phi, float psi, float* sample, int sampleIndex) {
    char ss = 'O';
    if((-67 <= phi && phi <= -47)   && (-57 <= psi && psi <= -37))  ss = 'A';
    if((-149 <= phi && phi <= -109) && (103 <= psi && psi <= 145))  ss = 'B';

    if(sampleIndex >= chis[ss].size()) return -999;

    sample[0] = chis[ss][sampleIndex].chi1;
    sample[1] = chis[ss][sampleIndex].chi2;
    sample[2] = chis[ss][sampleIndex].chi3;
    sample[3] = chis[ss][sampleIndex].chi4;
    return sampleIndex;
}

int PRLsampler::sample(float phi, float psi, float* sample) {
    char ss = 'O';
    if((-67 <= phi && phi <= -47)   && (-57 <= psi && psi <= -37))  ss = 'A';
    if((-149 <= phi && phi <= -109) && (103 <= psi && psi <= 145))  ss = 'B';

    int index = (int) floor(ran01() * r2ind[ss].size());
    int chiIndex = r2ind[ss][index];
    //for(int i=0; i < r2ind[ss].size(); i++) cout << i <<" : "<< r2ind[ss][i] << endl;
    //cout << "SEEEEEEE " << ss << " " << index << " " << chiIndex << endl;
    sample[0] = chis[ss][chiIndex].chi1;
    sample[1] = chis[ss][chiIndex].chi2;
    sample[2] = chis[ss][chiIndex].chi3;
    sample[3] = chis[ss][chiIndex].chi4;
    return chiIndex;
    //return ((int)ss) * chis[ss].size() + chiIndex;
}

void PRLsampler::addSample(float prop, vector<float> newchis) {
    CHI chi;
    chi.chi1 = newchis[0]; chi.chi2 = newchis[1]; chi.chi3 = newchis[2]; chi.chi4 = newchis[3];
    const char *states = "OAB";
    for(int i=0; i < 3; i++) {
        char s = states[i];
        chis[s].push_back(chi);
        int newplaces = int(floor( prop * r2ind[s].size() ));
        for(int k=0; k < newplaces; k++)
            r2ind[s].push_back(chis[s].size()-1);
    }
}

void PRLsampler::describe() {
    cout << chis.size() << " " << chis['A'].size() << " " << r2ind['A'].size() << endl;
}
PRLsampler PRLsampler::makeCopy() { return PRLsampler(*this); }