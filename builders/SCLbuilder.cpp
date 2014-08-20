#include "SCLbuilder.h"
#include "samplers/SCLsampler.h"
#include "geometry/functions.h"
#include "misc/verbosity.h"

SCLbuilder::SCLbuilder(vector<int>& ipInds, vector<int>& opInds, Constants *consts, const char* desc, SCLsampler *cs) 
        : Builder(ipInds, opInds, consts, desc), sclsampler(cs) {
    t1.resize(3); t2.resize(3); rot.resize(3);
    rot[0].resize(3); rot[1].resize(3); rot[2].resize(3);
}

void SCLbuilder::deleteSessionInfo() { sessionChis.clear(); }

bool SCLbuilder::checkAndAddIfNew(int ci) {
    if(sessionChis.size() >= maxSessionSize) return false;
    if(sessionChis.find(ci) != sessionChis.end()) return false;
    sessionChis.insert(ci); return true;
}

void SCLbuilder::getPhipsi(float & phi, float & psi, VVF & pts) {
    int iCp = ip[0], iN = ip[1], iCA = ip[2], iC = ip[3], iNn = ip[4];
    phi = calcDihed(pts[iCp], pts[iN], pts[iCA], pts[iC]);
    psi = calcDihed(pts[iN], pts[iCA], pts[iC], pts[iNn]);
}

// find phi,psi and get the sample
// find superposition operator for N,CA,C and copy CB + onwards
bool SCLbuilder::build(VVF & pts) {
    float phi, psi; getPhipsi(phi,psi,pts);
    vector<vector<float> > *rotpts = NULL; // pointer to rotamer points
    int ci = sclsampler->sample(int(phi), int(psi), &rotpts);

    if(session && !checkAndAddIfNew(ci)) {
      //if(verbose(7)) cout << "this SCL rotamer already sampled" << endl;
        return false; // sampled in this session already
    }
    return basicbuild(rotpts, pts);
}

bool SCLbuilder::basicbuild(VVF *rotpts, VVF & pts) {
    int iN = ip[1], iCA = ip[2], iC = ip[3], iCB = ip[5];

    vector<vector<float> > from, onto;
    from.push_back( (*rotpts)[0] ); onto.push_back( pts[iN] );
    from.push_back( (*rotpts)[1] ); onto.push_back( pts[iCA] );
    from.push_back( (*rotpts)[2] ); onto.push_back( pts[iC] );
    from.push_back( (*rotpts)[4] ); onto.push_back( pts[iCB] );

    findSuperpositionTransform(from, onto, t1, rot, t2);

    int opi;
    for(int i=0; i < op.size(); i++) {
        //cout << opi <<" "<< pts.size()<<" "<< 5+i <<" "<< (*rotpts).size() << endl;
        opi = op[i];
        pts[opi] = (*rotpts)[5+i]; // skip initial N,CA,C,O,CB crds from copying
        TRTtransform1(pts[opi], t1,rot,t2);
    }
    return true;
}

bool SCLbuilder::buildSample(VVF & pts, int sampleIndex) {
    float phi, psi; getPhipsi(phi,psi,pts);
    vector<vector<float> > *rotpts = NULL; // pointer to rotamer points
    int ci = sclsampler->sample(int(phi), int(psi), &rotpts, sampleIndex);
    if(ci < 0 || ci != sampleIndex) return false;
    return basicbuild(rotpts, pts);
}

SCLbuilder SCLbuilder::makeCopy() {
    return SCLbuilder(*this);
}
