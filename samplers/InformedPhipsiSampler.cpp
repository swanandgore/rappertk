#include <math.h>
#include <assert.h>

#include "InformedPhipsiSampler.h"
#include "RATdata.h"
#include "geometry/functions.h"
#include "misc/RanGen.h"
#include "misc/verbosity.h"

InformedPhipsiSampler::InformedPhipsiSampler(vector<float> & tp, float tol, const char* rn, const char *ratdatapath, bool fwdbwd)
        : resn(rn), targetTol(tol) {
    target = tp;
    targetTol = tol;
    resn = rn;
    RATdataPath = ratdatapath;
    fwd = fwdbwd;
}

bool targetSatisfied(vector<float>& p, vector<float>& cen, float rad) {
    float dx = p[0] - cen[0];
    if(fabs(dx) > rad) return false;
    float dy = p[1] - cen[1];
    if(fabs(dy) > rad) return false;
    float dz = p[2] - cen[2];
    if(fabs(dz) > rad) return false;

    if(dx*dx + dy*dy + dz*dz > rad*rad) return false;
    return true;
}

// if fwd, C,N,CA are C,N,CA, else assumed to be N,C,CA
void InformedPhipsiSampler::sample(
        vector<float>& curC, vector<float>& curN, vector<float>& curCA,
        float & phi, float & psi, float & omega, float & r, float & a, float & t) {
    float d0 = calcDist(target, curCA);

    float cisdist = 2.8, transdist = 3.81;

    bool cis = false, trans = false;
    if(cisdist > d0 - targetTol && cisdist < d0 + targetTol) cis = true;
    if(transdist > d0 - targetTol && transdist < d0 + targetTol) trans = true;

    if(cis && trans) { // if both cis and trans are possible, choose one
        if(ran01() > 0.95) trans = false;
    } else if(!cis && !trans) {
        cout << "target unreachable "
            <<curCA[0]<<" "<<curCA[1]<<" "<<curCA[2] << " : " << d0 <<" : "<< targetTol <<" : "<< target[0]<<" "<<target[1]<<" "<<target[2]
            << endl;
        exit(0);
    }
    // now either cis or trans is true

    float interCA;
    if(cis) { omega = 0; interCA = cisdist; }
    if(trans) { omega = -180; interCA = transdist; }

    phi = -999; psi = -999;

    // find a,t which take u into restraint sphere with this interCA
    vector<float> tar(3,0), noi(3,0);
    while(1) {
        // add some noise to target, proportional to noise
        randomNormalVector(noi);
        linear_combination(tar, 1, target, targetTol*ran01(), noi);
    
        // for this interCA, find a,t that takes you closest to target
        r = calcDist(curCA, tar);
        a = calcAngle(curN, curCA, tar);
        t = calcDihed(curC, curN, curCA, tar);
        find4thPoint(tar, curC, curN, curCA, interCA, a, t);
        if( calcDist(target, tar) <= targetTol ) break;
    }

    // return a phi-psi value for this a,t if available
    if(fwd) {
        RATdata & ratdata = RATdata::fwdinstance(RATdataPath.c_str());
        vector<PSO>::iterator bi, ei;
        ratdata.range(resn, interCA, a, t, bi, ei);
        int rs = ei - bi;
        assert(rs >= 0);
        if(rs == 0) { return; }
        int incr = (int) floor ( (0.+rs) * ran01() );
        bi += incr;
        phi = bi->r; psi = bi->a; omega = bi->t;
    } else {
        RATdata & ratdata = RATdata::bwdinstance(RATdataPath.c_str());
        vector<PSO>::iterator bi, ei;
        ratdata.range(resn, interCA, a, t, bi, ei);
        int rs = ei - bi;
        assert(rs >= 0);
        if(rs == 0) { return; }
        int incr = (int) floor ( (0.+rs) * ran01() );
        bi += incr;
        psi = bi->r; phi = bi->a; omega = bi->t;
    }
    //cout << phi <<" "<< psi <<" "<< omega << " shd take u from ("<< curCA[0]<<" "<< curCA[1]<<" "<< curCA[2] <<" ";
    //vector<float> exppt(3,0);
    //find4thPoint(exppt, curC, curN, curCA, interCA, a, t);
    //cout << ") to ("<< exppt[0]<<" "<< exppt[1]<<" "<< exppt[2] <<") " << endl;;
}

