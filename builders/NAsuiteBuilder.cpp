#include "NAsuiteBuilder.h"
#include "Constants.h"
#include "NAbuilder.h"
#include "samplers/NAsuiteSampler.h"
#include "geometry/functions.h"
#include "misc/verbosity.h"

#include <math.h>
#include <iostream>
#include <vector>
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>

using namespace std;

#include <assert.h>

NAsuiteBuilder::NAsuiteBuilder(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc,
        NAsuiteSampler* nas, char _base, bool _deoxy, char* _pucker, char* _puckerNext) : Builder(ipInds, opInds, con, desc) {
    naSampler = nas;

    base = _base;
    assert(base == 'A' || base == 'G' || base == 'T' || base == 'C' || base == 'U');

    deoxy = _deoxy;

    pucker = _pucker;
    assert(pucker == "C2'-endo" || pucker == "C3'-endo");

    puckerNext = _puckerNext;
    assert(puckerNext == "C2'-endo" || puckerNext == "C3'-endo");
}

// depending on this and next pucker, sample delta0(O3*), epsilon(P), zeta(O5*), alpha(C5*), beta(C4*), gamma(C3*)
// depending on pucker, build sugar and use chi(base)
// if rna, build O2*
// ips expected : C5*, C4*, C3*
// ops : O3*; next's (P, O5*, C5*, C4*, C3*, O1P, O2P); O4*, C1*, C2*, O2* if rna, followed by base atoms
// A : 0 - N1, 1 - C2, 2 - N3, 3 - C4, 4 - C5, 5 - C6, 6 - N7, 7 - C8, 8 - N9, 9 - N6
// G : 0 - N1, 1 - C2, 2 - N3, 3 - C4, 4 - C5, 5 - C6, 6 - N7, 7 - C8, 8 - N9, 9 - O6, 10 - N2
// T : 0 - N1, 1 - C2, 2 - O2, 3 - N3, 4 - C4, 5 - C5, 6 - C6, 7 - O4, 8 - C5M
// C : 0 - N1, 1 - C2, 2 - O2, 3 - N3, 4 - C4, 5 - C5, 6 - C6, 7 - N4
// U : 0 - N1, 1 - C2, 2 - O2, 3 - N3, 4 - C4, 5 - C5, 6 - C6, 7 - O4
bool NAsuiteBuilder::build(VVF & pts) {
    SuiteVal sv;
    int p0 = 2; if(pucker == "C3'-endo") p0 = 3;
    int p1 = 2; if(puckerNext == "C3'-endo") p1 = 3;
    naSampler->sample(p0, p1, sv);
    int delta = sv.d0;
    int epsilon = sv.e;
    int zeta = sv.z;
    int alpha = sv.a;
    int beta = sv.b;
    int gamma = sv.g;
    //    if(verbose(7)) { cout << "sampled "; sv.describe(cout); }
    if(session && !checkAndAddIfNew(sv)) return false;

    buildBackbone(consts, pts[ip[0]], pts[ip[1]], pts[ip[2]],
        pts[op[0]], pts[op[1]], pts[op[2]], pts[op[3]], pts[op[4]], pts[op[5]],
        delta, epsilon, zeta, alpha, beta, gamma);

    //    if(verbose(7)) cout << "built backbone" << endl;

    NAbuilder::buildPhosphate(consts, pts[op[0]], pts[op[1]], pts[op[2]], pts[op[6]], pts[op[7]]);

    vector<float> o2(3);
    if(deoxy)
        NAbuilder::buildSugar(consts, pucker, pts[ip[0]], pts[ip[1]], pts[ip[2]], pts[op[0]],
                pts[op[8]], pts[op[9]], pts[op[10]], o2);
    else
        NAbuilder::buildSugar(consts, pucker, pts[ip[0]], pts[ip[1]], pts[ip[2]], pts[op[0]],
                pts[op[8]], pts[op[9]], pts[op[10]], pts[op[11]]);

    //    if(verbose(7)) cout << "built sugar" << endl;

    int bs = 12;
    if(deoxy) bs = 11;
    if(base == 'T' || base == 'C' || base == 'U') {
        NAbuilder::build_TCU_scaffold( consts, pucker, base,
            pts[ip[1]], pts[op[8]], pts[op[9]],
            pts[op[bs+0]], pts[op[bs+1]], pts[op[bs+2]], pts[op[bs+3]], pts[op[bs+4]], pts[op[bs+5]], pts[op[bs+6]]
                // 0 - N1, 1 - C2, 2 - O2, 3 - N3, 4 - C4, 5 - C5, 6 - C6
        );
        if(base == 'T') NAbuilder::build_T( consts, pts[ op[bs+1] ], pts[ op[bs+3] ], pts[ op[bs+4] ], pts[ op[bs+5] ],
                        pts[ op[bs+7] ], pts[ op[bs+8] ] ); // ip: C2, N3, C4, C5 op: O4, C5M
        else if(base == 'C') NAbuilder::build_C( consts, pts[ op[bs+1] ], pts[ op[bs+3] ], pts[ op[bs+4] ],
                        pts[ op[bs+7] ]); // ip: C2, N3, C4 op: N4
        else NAbuilder::build_U( consts, pts[ op[bs+1] ], pts[ op[bs+3] ], pts[ op[bs+4] ],
                        pts[ op[bs+7] ] ); // ip: C2, N3, C4 op: O4
    }
    else if(base == 'A' || base == 'G') {
        NAbuilder::build_AG_scaffold( consts, pucker, base,
                pts[ ip[1] ], pts[ op[8] ], pts[ op[9] ], // c4s, o4s, c1s
                pts[ op[bs+0] ], pts[ op[bs+1] ], pts[ op[bs+2] ], pts[ op[bs+3] ], pts[ op[bs+4] ], pts[ op[bs+5] ],
                pts[ op[bs+6] ], pts[ op[bs+7] ], pts[ op[bs+8] ] // 0 - N1, 1 - C2, 2 - N3, 3 - C4, 4 - C5, 5 - C6, 6 - N7, 7 - C8, 8 - N9
        );
        if(base == 'A') NAbuilder::build_A( consts, pts[op[bs+1]], pts[op[bs+0]], pts[op[bs+5]], pts[op[bs+9]] ); // c2, n1, c6, n6
        else NAbuilder::build_G( consts, pts[op[bs+1]], pts[op[bs+0]], pts[op[bs+5]], pts[op[bs+9]], pts[op[bs+10]] ); // c2, n1, c6, o6, n2
    }
    else { cout << "Cannot build nucleotide base " << base << endl; assert(0); }

    return true;
}

void NAsuiteBuilder::buildBackbone( Constants *consts, vector<float>& c5s, vector<float>& c4s, vector<float>& c3s,
        vector<float>& o3s, vector<float>& pn, vector<float>& o5sn, vector<float>& c5sn, vector<float>& c4sn, vector<float>& c3sn,
        float delta, float epsilon, float zeta, float alpha, float beta, float gamma) {
    find4thPoint(  o3s,  c5s,  c4s,  c3s, cget("NU_C3*_O3*"), cget("NU_C4*_C3*_O3*"), delta );
    find4thPoint(  pn,  c4s,  c3s,  o3s, cget("NU_P_O3*"), cget("NU_P_O3*_C3*"), epsilon );
    find4thPoint( o5sn,  c3s,  o3s,  pn, cget("NU_P_O5*"), cget("NU_O3*_P_O5*"), zeta );
    find4thPoint(  c5sn, o3s,   pn,  o5sn, cget("NU_O5*_C5*"), cget("NU_P_O5*_C5*"), alpha );
    find4thPoint(  c4sn,   pn,  o5sn,  c5sn, cget("NU_C5*_C4*"), cget("NU_O5*_C5*_C4*"), beta );
    find4thPoint(  c3sn,  o5sn,  c5sn,  c4sn, cget("NU_C4*_C3*"), cget("NU_C5*_C4*_C3*"), gamma );
}



void NAsuiteBuilder::deleteSessionInfo() { sessionSVs.clear(); }

bool NAsuiteBuilder::checkAndAddIfNew(SuiteVal & sv) {
    for(int i=0; i < sessionSVs.size(); i++)
        if(
            abs(sessionSVs[i].d0  - sv.d0) < 2 &&
            abs(sessionSVs[i].a  - sv.a) < 2 &&
            abs(sessionSVs[i].b  - sv.b) < 2 &&
            abs(sessionSVs[i].g  - sv.g) < 2 &&
            abs(sessionSVs[i].e  - sv.e) < 2 &&
            abs(sessionSVs[i].z  - sv.z) < 2 &&
            abs(sessionSVs[i].d1  - sv.d1) < 2
        ) return false;

    sessionSVs.push_back(sv);
    return true;
}

NAsuiteBuilder NAsuiteBuilder::makeCopy() { return NAsuiteBuilder(*this); }
