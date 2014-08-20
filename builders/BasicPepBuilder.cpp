#include <assert.h>
#include <math.h>
#include "BasicPepBuilder.h"
#include "Constants.h"
#include "samplers/PhipsiSampler.h"
#include "samplers/OmegaSampler.h"
#include "geometry/functions.h"
#include "misc/verbosity.h"
#include "misc/RanGen.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
BasicPepBuilder::BasicPepBuilder(vector<int>& ipInds, vector<int>& opInds, Constants* con, const char *resn, const char *desc,
        bool fwdbkwd) : Builder(ipInds, opInds, con, desc) {

    assert(ipInds.size() == 3 && (opInds.size() == 4 || opInds.size() == 5));

    fwd = fwdbkwd;

    CA_C = consts->get("CA_C");
    N_CA_C = consts->get("N_CA_C");
    C_N = consts->get("C_N");
    CA_C_N = consts->get("CA_C_N");
    N_CA = consts->get("N_CA");
    C_N_CA = consts->get("C_N_CA");
    C_O = consts->get("C_O");
    CA_C_O = consts->get("CA_C_O");
    N_C_O = consts->get("N_C_O");
    CA_CB = consts->get("CA_CB");
    if(strcmp(resn,"PRO")==0) N_CA_CB = consts->get("PRO_N_CA_CB");
    else N_CA_CB = consts->get("N_CA_CB");
    C_N_CA_CB = consts->get("C_N_CA_CB");
}

void BasicPepBuilder::deleteSessionInfo() { sessionPSO.clear(); }

bool BasicPepBuilder::checkAndAddIfNew(float phi, float psi, float omega) {
    if(sessionSize >= maxSessionSize) return false;
    //assuming phi,psi to be within -180,180. and omega to be -180 or 0. XXX
    assert(-180 <= phi && phi < 180);
    assert(-180 <= psi && psi < 180);
    assert(omega == 0 or omega == -180);
    long long int key = round(180+phi) + round(180+psi) * 72 + 72*72*(180+omega);
    sessionPSO.insert(key);
    if(sessionPSO.size() == sessionSize) return false;
    sessionSize ++; return true;
}

void BasicPepBuilder::sample(VVF & pts, float& phi, float& psi, float& omega) { exit(1); }

// fwd : ip[0,1,2] = C,N,CA op[0,1,2,3] = C,O,N,CA
// bwd : ip[0,1,2] = N,C,CA op[0,1,2,3] = N,C,O,CA
bool BasicPepBuilder::build(VVF & pts) {
    float phi=-999, psi=-999, omega=-999;
    sample(pts, phi, psi, omega); // this function can be overridden for sampling differently
    //cout << "phipsi sample " << phi <<" "<< psi <<" "<< omega << endl;
    if(phi < -200) { // revert to backup sampler
      //if(verbose(7)) cout << "sampling failed" << endl;
        return false;
    }

    //    if(verbose(7)) cout << "phipsi sample " << phi <<" "<< psi <<" "<< omega << endl;

    if(session) {
        if( ! checkAndAddIfNew(phi,psi,omega) ) {
	  //            if(verbose(7)) cout << "phipsi already sampled or session size exceeded" << endl;
            return false; // sampled in this session already
        }
    }
    return build1(pts, phi, psi, omega); // this can be used for using build functionality w/o sampling functionality
}

bool BasicPepBuilder::build1(VVF & pts, float phi, float psi, float omega) {
    if(fwd) return buildFwd(pts, phi, psi, omega);
    return buildBkwd(pts, phi, psi, omega);
}


bool BasicPepBuilder::buildBkwd(VVF & pts, float phi, float psi, float omega) {
    int iN = ip[0], iC = ip[1], iCA = ip[2], oN = op[0], oC = op[1], oO = op[2], oCA = op[3];
    //cout << "IO pts " << iN <<" "<< iC <<" "<< iCA <<" "<< oN <<" "<< oC <<" "<< oO <<" "<< oCA <<" "<< consts << endl;

    find4thPoint( pts[ oN ], // N
        pts[ iN ], pts[ iC ], pts[ iCA ],
        N_CA, N_CA_C, psi );
    find4thPoint( pts[ oC ], // C
        pts[ iC ], pts[ iCA ], pts[ oN ],
        C_N, C_N_CA, phi );
    find4thPoint( pts[ oO ], // O
        pts[ iCA ], pts[ oN ], pts[ oC ],
        C_O, N_C_O, withinPlusMinus180(omega+180) );

    if(op.size() > 3)
        find4thPoint( pts[ oCA ], // CA
            pts[ iCA ], pts[ oN ], pts[ oC ],
            CA_C, CA_C_N, omega );

    if(op.size() > 4) // build CB also
        find4thPoint( pts[ op[4] ],
            pts[ iC ], pts[ oN ], pts[ iCA ],
            CA_CB, N_CA_CB, C_N_CA_CB );
    //cout << pts[oC][0] <<" "<< pts[oC][1] <<" "<< pts[oC][2] << endl;
    //cout << pts[oO][0] <<" "<< pts[oO][1] <<" "<< pts[oO][2] << endl;
    //cout << pts[oN][0] <<" "<< pts[oN][1] <<" "<< pts[oN][2] << endl;
    //cout << pts[oCA][0] <<" "<< pts[oCA][1] <<" "<< pts[oCA][2] << endl;
    return true;
}

bool BasicPepBuilder::buildFwd(VVF & pts, float phi, float psi, float omega) {
    //cout << CA_C <<" "<< N_CA_C <<" "<< C_N <<" "<< CA_C_N <<" "<< N_CA <<" "<< C_N_CA <<" "<< C_O <<" "<< CA_C_O << endl;
    int iC = ip[0], iN = ip[1], iCA = ip[2], oC = op[0], oO = op[1], oN = op[2], oCA = op[3];
    //cout << "IO pts " << iC <<" "<< iN <<" "<< iCA <<" "<< oC <<" "<< oO <<" "<< oN <<" "<< oCA <<" "<< consts << endl;
    //cout << "fwd building " << pts[oCA][0] <<" "<< pts[oCA][1] <<" "<< pts[oCA][2] << endl;

    find4thPoint( pts[ oC ], // C
        pts[ iC ], pts[ iN ], pts[ iCA ],
        CA_C, N_CA_C, phi );
    find4thPoint( pts[ oN ], // N
        pts[ iN ], pts[ iCA ], pts[ oC ],
        C_N, CA_C_N, psi );
    find4thPoint( pts[ oO ], // O
        pts[ iN ], pts[ iCA ], pts[ oC ],
        C_O, CA_C_O, withinPlusMinus180(psi+180) );

    if(op.size() > 3)
        find4thPoint( pts[ oCA ], // CA
            pts[ iCA ], pts[ oC ], pts[ oN ],
            N_CA, C_N_CA, omega );

    if(op.size() > 4) // build CB also
        find4thPoint( pts[ op[4] ],
            pts[ oC ], pts[ iN ], pts[ iCA ],
            CA_CB, N_CA_CB, C_N_CA_CB );
    //cout << pts[oC][0] <<" "<< pts[oC][1] <<" "<< pts[oC][2] << endl;
    //cout << pts[oO][0] <<" "<< pts[oO][1] <<" "<< pts[oO][2] << endl;
    //cout << pts[oN][0] <<" "<< pts[oN][1] <<" "<< pts[oN][2] << endl;
    //cout << pts[oCA][0] <<" "<< pts[oCA][1] <<" "<< pts[oCA][2] << endl;
    // cout << "fwd done " << pts[oCA][0] <<" "<< pts[oCA][1] <<" "<< pts[oCA][2] << endl;
    return true;
}
