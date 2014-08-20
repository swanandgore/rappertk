#include <iostream>
using namespace std;

#include "InformedPeptideBuilder.h"
#include "samplers/InformedPhipsiSampler.h"
#include "samplers/PhipsiSampler.h"
#include "samplers/OmegaSampler.h"
#include "misc/verbosity.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>

InformedPeptideBuilder::InformedPeptideBuilder(vector<int>& ipInds, vector<int>& opInds,
	    Constants* con, const char *resn, const char* desc, PhipsiSampler *PSS, OmegaSampler *OS, InformedPhipsiSampler *infpps, bool fwdbkwd)
        : BasicPepBuilder(ipInds, opInds, con, resn, desc, fwdbkwd)
{
    ipps = infpps;
    pss = PSS;
    os = OS;
}

void InformedPeptideBuilder::sample(VVF & pts, float &phi, float &psi, float &omega) {
    phi = -999; psi = -999; omega = -999;
    ipps->sample ( pts[ip[0]], pts[ip[1]], pts[ip[2]], phi, psi, omega, r, a, t);
}

InformedPeptideBuilder InformedPeptideBuilder::makeCopy() {
    return InformedPeptideBuilder(*this);
}

bool InformedPeptideBuilder::checkAndAddIfNew(float r, float a, float t) {
    if(sessionRATs.size() > sessionSize) return false;
    r *= 100; a *= 10; t = (t+180)*1;
    long long int key = round(a) + 180*10*round(t) + 180*10*360*round(r);
    //r = round(r*20); a = round(a*10); t = round((t+180)*10);
    //long long int key = a + 180*10*t + 180*10*360*10*r;
    if(key < 0) {
        cout << "key " << r <<" "<< a <<" "<< t <<" "<< key <<endl;
        exit(0);
    }
    sessionRATs.insert(key);
    if(sessionRATs.size() == sessionSize) return false;
    sessionSize ++; return true;
}

void InformedPeptideBuilder::deleteSessionInfo() {
    sessionRATs.clear();   
}

bool InformedPeptideBuilder::build(VVF & pts) {
    float phi=-999, psi=-999, omega=-999;
    sample(pts, phi, psi, omega); // this function can be overridden for sampling differently

    //    if(verbose(7)) cout << "phipsi sample " << phi <<" "<< psi <<" "<< omega << endl;

    if(session) {
        for(int count=0; count < 10000; count++) // dont look too hard for a novel sample
            if( checkAndAddIfNew(r,a,t) ) break;
    }
    if(phi < -200) { //no phi-psi to get to r,a,t. copy ip onto op coordinates
        pts[op[1]] = pts[ip[0]];
        pts[op[2]] = pts[ip[1]];
        pts[op[3]] = pts[ip[2]];
        return true;
    }

    return build1(pts, phi, psi, omega); // this can be used for using build functionality w/o sampling functionality
}
