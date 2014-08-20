#include "NAbuilder.h"
#include "Constants.h"
#include "samplers/NAsampler.h"
#include "geometry/functions.h"

#include <assert.h>

#include <vector>
#include <iostream>
using namespace std;

NAbuilder::NAbuilder(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc,
        NAsampler* nas, char _base, bool _deoxy, char* _pucker)
    : Builder(ipInds, opInds, con, desc) {

    pucker = _pucker;
    assert(pucker == "C2'-endo" || pucker == "C3'-endo");

    base = _base;
    assert(base == 'A' || base == 'G' || base == 'T' || base == 'C' || base == 'U');

    naSampler = nas;
    deoxy = _deoxy;
}


// ips expected : O3*, P, O5*
// ops : 0 - O1P, 1 - O2P, 2 - C5*, 3 - C4*, 4 - C3*, 5 - O3*, 6 - next P, 7 - next O5*, 8 - O4*, 9 - C1*, 10 - C2*, 11 - O2* (if RNA), followed by base atoms
// A : 0 - N1, 1 - C2, 2 - N3, 3 - C4, 4 - C5, 5 - C6, 6 - N7, 7 - C8, 8 - N9, 9 - N6
// G : 0 - N1, 1 - C2, 2 - N3, 3 - C4, 4 - C5, 5 - C6, 6 - N7, 7 - C8, 8 - N9, 9 - O6, 10 - N2
// T : 0 - N1, 1 - C2, 2 - O2, 3 - N3, 4 - C4, 5 - C5, 6 - C6, 7 - O4, 8 - C5M
// C : 0 - N1, 1 - C2, 2 - O2, 3 - N3, 4 - C4, 5 - C5, 6 - C6, 7 - N4
// U : 0 - N1, 1 - C2, 2 - O2, 3 - N3, 4 - C4, 5 - C5, 6 - C6, 7 - O4
bool NAbuilder::build(VVF & pts) {
    naSampler->sample(base, pucker, alpha, beta, gamma, epsilon, zeta);
    
    NAbuilder::buildPhosphate( consts, pts[ip[0]], pts[ip[1]], pts[ip[2]], pts[op[0]], pts[op[1]] );

    float delta = cget("NU_C3E_C5*_C4*_C3*_O3*");
    if(pucker == "C2'-endo") delta = cget("NU_C2E_C5*_C4*_C3*_O3*");
    buildBackbone( consts, pts[ip[0]], pts[ip[1]], pts[ip[2]], // o3*, p, o5*
                    pts[op[2]], pts[op[3]], pts[op[4]], pts[op[5]], pts[op[6]], pts[op[7]],
                    alpha, beta, gamma, delta, epsilon, zeta );

    vector<float> o2(3);
    if(deoxy)
        buildSugar( consts, pucker, pts[op[2]], pts[op[3]], pts[op[4]], pts[op[5]], // c5*, c4*, c3*, o3*
                    pts[op[8]], pts[op[9]], pts[op[10]], o2 ); // o4*, c1*, c2*, o2*-dummy
    else
        buildSugar( consts, pucker, pts[op[2]], pts[op[3]], pts[op[4]], pts[op[5]], // c5*, c4*, c3*, o3*
                    pts[op[8]], pts[op[9]], pts[op[10]], pts[op[11]] ); // o4*, c1*, c2*, o2*

    int bs = 11; // base point indices start at bs in ops array, ie pts[ ops[bs] ] is in base
    if(!deoxy) bs = 12;

    if(base == 'T' || base == 'C' || base == 'U') {
        build_TCU_scaffold( consts, pucker, base,
            pts[op[3]], pts[op[8]], pts[op[9]],
            pts[op[bs+0]], pts[op[bs+1]], pts[op[bs+2]], pts[op[bs+3]], pts[op[bs+4]], pts[op[bs+5]], pts[op[bs+6]]
                // 0 - N1, 1 - C2, 2 - O2, 3 - N3, 4 - C4, 5 - C5, 6 - C6
        );
        if(base == 'T') build_T( consts, pts[ op[bs+1] ], pts[ op[bs+3] ], pts[ op[bs+4] ], pts[ op[bs+5] ],
                        pts[ op[bs+7] ], pts[ op[bs+8] ] ); // ip: C2, N3, C4, C5 op: O4, C5M
        else if(base == 'C') build_C( consts, pts[ op[bs+1] ], pts[ op[bs+3] ], pts[ op[bs+4] ],
                        pts[ op[bs+7] ]); // ip: C2, N3, C4 op: N4
        else build_U( consts, pts[ op[bs+1] ], pts[ op[bs+3] ], pts[ op[bs+4] ],
                        pts[ op[bs+7] ] ); // ip: C2, N3, C4 op: O4
    }
    else if(base == 'A' || base == 'G') {
        build_AG_scaffold( consts, pucker, base, pts[ op[3] ], pts[ op[8] ], pts[ op[9] ], // c4s, o4s, c1s
                pts[ op[bs+0] ], pts[ op[bs+1] ], pts[ op[bs+2] ], pts[ op[bs+3] ], pts[ op[bs+4] ], pts[ op[bs+5] ],
                pts[ op[bs+6] ], pts[ op[bs+7] ], pts[ op[bs+8] ] // 0 - N1, 1 - C2, 2 - N3, 3 - C4, 4 - C5, 5 - C6, 6 - N7, 7 - C8, 8 - N9
        );
        if(base == 'A') build_A( consts, pts[op[bs+1]], pts[op[bs+0]], pts[op[bs+5]], pts[op[bs+9]] ); // c2, n1, c6, n6
        else build_G( consts, pts[op[bs+1]], pts[op[bs+0]], pts[op[bs+5]], pts[op[bs+9]], pts[op[bs+10]] ); // c2, n1, c6, o6, n2
    }
    else { cout << "Cannot build nucleotide base " << base << endl; assert(0); }

    return true;
}

void NAbuilder::buildPhosphate(Constants *consts, vector<float>& o3, vector<float>& p, vector<float>& o5, vector<float>& o1p, vector<float>& o2p) {
    find4thPoint( o1p, o3, o5, p, cget("NU_P_O1P"), cget("NU_O5*_P_O1P"), cget("NU_O3*_O5*_P_O1P") );
    find4thPoint( o2p, o3, o5, p, cget("NU_P_O2P"), cget("NU_O5*_P_O2P"), cget("NU_O3*_O5*_P_O2P") );
}

// o3p = o3 of prev 'residue', pn o5n = p, o5 of next 'residue'
void NAbuilder::buildBackbone( Constants *consts, vector<float>& o3p, vector<float>& p, vector<float>& o5,
            vector<float>& c5, vector<float>& c4, vector<float>& c3, vector<float>& o3, vector<float>& pn, vector<float>& o5n,
            float alpha, float beta, float gamma, float delta, float epsilon, float zeta ) {

    find4thPoint(  c5, o3p,   p,  o5, cget("NU_O5*_C5*"), cget("NU_P_O5*_C5*"), alpha );
    find4thPoint(  c4,   p,  o5,  c5, cget("NU_C5*_C4*"), cget("NU_O5*_C5*_C4*"), beta );
    find4thPoint(  c3,  o5,  c5,  c4, cget("NU_C4*_C3*"), cget("NU_C5*_C4*_C3*"), gamma );
    find4thPoint(  o3,  c5,  c4,  c3, cget("NU_C3*_O3*"), cget("NU_C4*_C3*_O3*"), delta );
    find4thPoint(  pn,  c4,  c3,  o3, cget("NU_P_O3*"), cget("NU_P_O3*_C3*"), epsilon );
    find4thPoint( o5n,  c3,  o3,  pn, cget("NU_P_O5*"), cget("NU_O3*_P_O5*"), zeta );
}

void NAbuilder::buildSugar( Constants *consts, string& pucker, vector<float>& c5, vector<float>& c4, vector<float>& c3, vector<float>& o3,
            vector<float>& o4, vector<float>& c1, vector<float>& c2, vector<float>& o2 ) {

    assert(pucker == "C2'-endo" || pucker == "C3'-endo");
    float dihed0 = cget("NU_C2E_O4*_C4*_C3*_O3*");
    float dihed1 = cget("NU_C2E_C3*_C4*_O4*_C1*");
    float dihed2 = cget("NU_C2E_C4*_O4*_C1*_C2*");
    float dihed3 = cget("NU_C2E_O3*_C3*_C2*_O2*");
    if(pucker == "C3'-endo") {
        dihed0 = cget("NU_C3E_O4*_C4*_C3*_O3*");
        dihed1 = cget("NU_C3E_C3*_C4*_O4*_C1*");
        dihed2 = cget("NU_C3E_C4*_O4*_C1*_C2*");
        dihed3 = cget("NU_C3E_O3*_C3*_C2*_O2*");
    }

    find4thPoint(o4, o3, c3, c4, cget("NU_C4*_C3*"), cget("NU_O4*_C4*_C3*"), dihed0);
    find4thPoint(c1, c3, c4, o4, cget("NU_C1*_O4*"), cget("NU_C1*_O4*_C4*"), dihed1);
    find4thPoint(c2, c4, o4, c1, cget("NU_C2*_C1*"), cget("NU_O4*_C1*_C2*"), dihed2);
    find4thPoint(o2, o3, c3, c2, cget("NU_C2*_O2*"), cget("NU_C3*_C2*_O2*"), dihed3);
}

void NAbuilder::build_AG_scaffold( Constants *consts, string & pucker, char base, vector<float>& c4s, vector<float>& o4s, vector<float>& c1s,
        vector<float>& n1, vector<float>& c2, vector<float>& n3, vector<float>& c4, vector<float>& c5, vector<float>& c6,
        vector<float>& n7, vector<float>& c8, vector<float>& n9 ) {
//cout << "in build_AG_scaffold" << endl;
    string pstr = "NU_C2E";
    if(pucker == "C3'-endo") pstr = "NU_C3E";

    float dist, angle, dihed;
    string lstr, astr, dstr;
    string bstr = "NU_"; bstr += base;

    dstr = pstr + "_" + base + "_C4*_O4*_C1*_N9"; dihed = cget(dstr);
    lstr = bstr + "_N9_C1*"; dist = cget(lstr);
    astr = bstr + "_O4*_C1*_N9"; angle = cget(astr);
    find4thPoint(n9, c4s, o4s, c1s, dist, angle, dihed);

    dstr = pstr + "_" + base + "_O4*_C1*_N9_C4"; dihed = cget(dstr);
    lstr = bstr + "_N9_C4"; dist = cget(lstr);
    astr = bstr + "_C1*_N9_C4"; angle = cget(astr);
    find4thPoint(c4, o4s, c1s, n9, dist, angle, dihed);

    lstr = bstr + "_C4_C5"; dist = cget(lstr);
    astr = bstr + "_N9_C4_C5"; angle = cget(astr);
    find4thPoint(c5, c1s, n9, c4, dist, angle, 180);

    lstr = bstr + "_C5_N7"; dist = cget(lstr);
    astr = bstr + "_C4_C5_N7"; angle = cget(astr);
    find4thPoint(n7, n9, c4, c5, dist, angle, 0);

    lstr = bstr + "_N7_C8"; dist = cget(lstr);
    astr = bstr + "_C5_N7_C8"; angle = cget(astr);
    find4thPoint(c8, c4, c5, n7, dist, angle, 0);

    lstr = bstr + "_C5_C6"; dist = cget(lstr);
    astr = bstr + "_C4_C5_C6"; angle = cget(astr);
    find4thPoint(c6, n9, c4, c5, dist, angle, 180);

    lstr = bstr + "_C6_N1"; dist = cget(lstr);
    astr = bstr + "_C5_C6_N1"; angle = cget(astr);
    find4thPoint(n1, c4, c5, c6, dist, angle, 0);

    lstr = bstr + "_N1_C2"; dist = cget(lstr);
    astr = bstr + "_C6_N1_C2"; angle = cget(astr);
    find4thPoint(c2, c5, c6, n1, dist, angle, 0);

    lstr = bstr + "_C2_N3"; dist = cget(lstr);
    astr = bstr + "_N1_C2_N3"; angle = cget(astr);
    find4thPoint(n3, c6, n1, c2, dist, angle, 0);
}

void NAbuilder::build_A( Constants *consts, vector<float>& c2, vector<float>& n1, vector<float>& c6, vector<float>& n6 ) {
    find4thPoint(n6, c2, n1, c6, cget("NU_A_C6_N6"), cget("NU_A_N1_C6_N6"), 180);
}

void NAbuilder::build_G( Constants *consts, vector<float>& c2, vector<float>& n1, vector<float>& c6, vector<float>& o6, vector<float>& n2 ) {
    find4thPoint(o6, c2, n1, c6, cget("NU_G_C6_O6"), cget("NU_G_N1_C6_O6"), 180);
    find4thPoint(n2, c6, n1, c2, cget("NU_G_C2_N2"), cget("NU_G_N1_C2_N2"), 180);
}


void NAbuilder::build_TCU_scaffold( Constants *consts, string & pucker, char base,
        vector<float>& c4s, vector<float>& o4s, vector<float>& c1s, // ips, s for sugar
        vector<float>& n1, vector<float>& c2, vector<float>& o2, vector<float>& n3, vector<float>& c4, vector<float>& c5, vector<float>& c6 ) { // ops

    string pstr = "NU_C2E";
    if(pucker == "C3'-endo") pstr = "NU_C3E";

    float dist, angle, dihed;
    string lstr, astr, dstr;
    string bstr = "NU_"; bstr += base;

    dstr = pstr + "_" + base + "_C4*_O4*_C1*_N1"; dihed = cget(dstr);
    lstr = bstr + "_N1_C1*"; dist = cget(lstr);
    astr = bstr + "_O4*_C1*_N1"; angle = cget(astr);
    find4thPoint(n1, c4s, o4s, c1s, dist, angle, dihed);

    dstr = pstr + "_" + base + "_O4*_C1*_N1_C2"; dihed = cget(dstr);
    lstr = bstr + "_N1_C2"; dist = cget(lstr);
    astr = bstr + "_C1*_N1_C2"; angle = cget(astr);
    find4thPoint(c2, o4s, c1s, n1, dist, angle, dihed);

    lstr = bstr + "_C2_O2"; dist = cget(lstr);
    astr = bstr + "_N1_C2_O2"; angle = cget(astr);
    find4thPoint(o2, c1s, n1, c2, dist, angle, 0); // anti is favored in pyrimidines

    lstr = bstr + "_C2_N3"; dist = cget(lstr);
    astr = bstr + "_N1_C2_N3"; angle = cget(astr);
    find4thPoint(n3, c1s, n1, c2, dist, angle, 180);

    lstr = bstr + "_N3_C4"; dist = cget(lstr);
    astr = bstr + "_C2_N3_C4"; angle = cget(astr);
    find4thPoint(c4, n1, c2, n3, dist, angle, 0);

    lstr = bstr + "_C4_C5"; dist = cget(lstr);
    astr = bstr + "_N3_C4_C5"; angle = cget(astr);
    find4thPoint(c5, c2, n3, c4, dist, angle, 0);

    lstr = bstr + "_C5_C6"; dist = cget(lstr);
    astr = bstr + "_C4_C5_C6"; angle = cget(astr);
    find4thPoint(c6, n3, c4, c5, dist, angle, 0);
}

void NAbuilder::build_T(Constants *consts, vector<float>& c2, vector<float>& n3, vector<float>& c4, vector<float>& c5, vector<float>& o4, vector<float>& c5m) {
    find4thPoint(o4, c2, n3, c4, cget("NU_T_C4_O4"), cget("NU_T_N3_C4_O4"), 180);
    find4thPoint(c5m, n3, c4, c5, cget("NU_T_C5_C5M"), cget("NU_T_C4_C5_C5M"), 180);
}
void NAbuilder::build_C(Constants *consts, vector<float>& c2, vector<float>& n3, vector<float>& c4, vector<float>& n4) {
    find4thPoint(n4, c2, n3, c4, cget("NU_C_C4_N4"), cget("NU_C_N3_C4_N4"), 180);
}
void NAbuilder::build_U(Constants *consts, vector<float>& c2, vector<float>& n3, vector<float>& c4, vector<float>& o4) {
    find4thPoint(o4, c2, n3, c4, cget("NU_U_C4_O4"), cget("NU_U_N3_C4_O4"), 180);
}
