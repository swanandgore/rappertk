#include <math.h>

#include <vector>
using namespace std;

#include "RapSampler.h"
#include "misc/RanGen.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>

RapSampler::RapSampler(const char *filename, const char *resname) {
    int phi, psi;
    char line[500], ignme[100], resn[100];
    FILE *fp = fopen(filename, "r");
    while(fgets(line, 500, fp) != NULL) {
        sscanf(line, "%s %d %d", resn, &phi, &psi);
        if(strcmp(resn, resname) != 0) continue;
        phiseq.push_back(phi);
        psiseq.push_back(psi);
        if(phiseq.size() == 100000) break;
    }
    seqnum = 0;
}

void RapSampler::sample(float & phi, float & psi) {
    int si = (int) floor(ran01() * phiseq.size());
    phi = phiseq[si];
    psi = psiseq[si];
    //seqnum = (seqnum+1) % phiseq.size();
}
