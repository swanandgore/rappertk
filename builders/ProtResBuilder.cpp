#include "ProtResBuilder.h"
#include "BasicPepBuilder.h"
#include "CBbuilder.h"
#include "ChiBuilder.h"
#include "samplers/PhipsiSampler.h"
#include "samplers/OmegaSampler.h"
#include "samplers/BBdepChiSampler.h"
#include "misc/verbosity.h"
#include "geometry/functions.h"

#include <iostream>
using namespace std;

//fwd : ip - C N CA, op - C O N CA CB sidechain atoms
//bwd : ip - N C CA, op - N C O CA CB sidechain atoms
ProtResBuilder::ProtResBuilder(vector<int>& ipInds, vector<int>& opInds, Constants *consts, const char* desc, const char *resname,
        PhipsiSampler *p, OmegaSampler *o, BBdepChiSampler *cs, vector<float>& CA, float RAD, bool fwdbkwd)
        : Builder(ipInds, opInds, consts, desc), pps(p), os(o), chisampler(cs), ca(CA), carad(RAD) {
    vector<int> pepIP, pepOP;
    pepIP.push_back(ipInds[0]); pepIP.push_back(ipInds[1]); pepIP.push_back(ipInds[2]);
    pepOP.push_back(opInds[0]); pepOP.push_back(opInds[1]); pepOP.push_back(opInds[2]); pepOP.push_back(opInds[3]);
    pb = new BasicPepBuilder(pepIP, pepOP, consts, resname, desc, fwdbkwd);

    vector<int> cbIP, cbOP;
    if(fwdbkwd) { cbIP.push_back(ipInds[1]); cbIP.push_back(opInds[0]); cbIP.push_back(ipInds[2]); }
    else { cbIP.push_back(opInds[0]); cbIP.push_back(ipInds[1]); cbIP.push_back(ipInds[2]); }
    cbOP.push_back(opInds[4]);
    cbb = new CBbuilder(cbIP, cbOP, consts, desc, resname);

    vector<int> chiIP, chiOP;
    if(fwdbkwd) { chiIP.push_back(ipInds[0]); chiIP.push_back(ipInds[1]); chiIP.push_back(ipInds[2]); chiIP.push_back(opInds[0]); chiIP.push_back(opInds[2]); }
    else { chiIP.push_back(opInds[1]); chiIP.push_back(opInds[0]); chiIP.push_back(ipInds[2]); chiIP.push_back(ipInds[1]); chiIP.push_back(ipInds[0]); }
    chiIP.push_back(opInds[4]);
    for(int i=5; i < opInds.size(); i++) chiOP.push_back(opInds[i]);
    chib = new ChiBuilder(chiIP, chiOP, consts, desc, resname, NULL); // null sampler
}

ProtResBuilder::~ProtResBuilder() {}

bool ProtResBuilder::build(VVF & pts) {
    int phi = -9999, psi = -9999, omega = -9999; float chis[4];
    chis[0] = -9999; chis[1] = -9999; chis[2] = -9999; chis[3] = -9999;
    int pi = pps->sample(phi,psi);
    int oi = os->sample(omega);
    if(psoBlacklist.find(pi + 5184*oi) != psoBlacklist.end()) return false;
    int ci = chisampler->sample(phi,psi,chis);
    if(session && !checkAndAddIfNew(pi,oi,ci)) {
      //if(verbose(7)) cout << "phi-psi-chi combination already sampled or sessionSize exceeded" << endl;
        return false;
    }
    //    if(verbose(7)) cout << "PSOC " << phi <<" "<< psi <<" "<< omega <<" "<< chis[0] <<" "<< chis[1] <<" "<< chis[2] <<" "<< chis[3] << endl;

    pb->build1(pts, phi, psi, omega);
    if(calcDist(ca,pts[op[3]]) > carad) {
      //if(verbose(7)) cout << "Protres CA failure (" << ca[0]<<","<<ca[1]<<","<<ca[2]<<") " << carad<< endl;
      psoBlacklist.insert(pi + 5184*oi);
        return true; // no need to sample further
    }
    cbb->build(pts);
    chib->build(pts, chis);
    return true;
}

ProtResBuilder ProtResBuilder::makeCopy() {
    return ProtResBuilder(*this);
}

void ProtResBuilder::deleteSessionInfo() { sessionIds.clear(); psoBlacklist.clear(); }

PsoChis::PsoChis(int p, int s, int o, int c1, int c2, int c3, int c4)
    : phi(p), psi(s), omega(o), chi1(c1), chi2(c2), chi3(c3), chi4(c4) {}
    
bool ProtResBuilder::checkAndAddIfNew(int psi, int osi, int csi) {
    if(sessionSize >= maxSessionSize) return false;
    long id = psi + 5184*osi + 5184*2*csi; 
    sessionIds.insert(id);
    if(sessionSize == sessionIds.size()) return false;
    sessionSize ++;
    return true;
}

bool PsoChis::operator()(const PsoChis & psoca, const PsoChis & psocb) { // return true if <
    if(psoca.phi < psocb.phi) return true;
    else if(psoca.phi > psocb.phi) return false;
    if(psoca.psi < psocb.psi) return true;
    else if(psoca.psi > psocb.psi) return false;
    if(psoca.omega < psocb.omega) return true;
    else if(psoca.omega > psocb.omega) return false;
    if(psoca.chi1 < psocb.chi1) return true;
    else if(psoca.chi1 > psocb.chi1) return false;
    if(psoca.chi2 < psocb.chi2) return true;
    else if(psoca.chi2 > psocb.chi2) return false;
    if(psoca.chi3 < psocb.chi3) return true;
    else if(psoca.chi3 > psocb.chi3) return false;
    if(psoca.chi4 < psocb.chi4) return true;
    else if(psoca.chi4 > psocb.chi4) return false;
    return false;
}
