#include <assert.h>
#include "CBbuilder.h"
#include "geometry/functions.h"
#include "misc/verbosity.h"
#include "Constants.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>

CBbuilder::CBbuilder(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc, const char *resn)
	: Builder(ipInds, opInds, con, desc) {
    assert(ipInds.size() == 3);
    assert(opInds.size() == 1);
    CA_CB = consts->get("CA_CB");
    if(strcmp(resn,"PRO") == 0) N_CA_CB = consts->get("PRO_N_CA_CB");
    else                        N_CA_CB = consts->get("N_CA_CB");
    C_CA_CB = consts->get("C_CA_CB");
    N_C_CA_CB = consts->get("N_C_CA_CB");
    C_N_CA_CB = consts->get("C_N_CA_CB");
}

bool CBbuilder::buildSample(VVF & pts, int sampleIndex) {
    if(sampleIndex != 0) return false;
    build(pts);
    return true;
}

bool CBbuilder::build(VVF & pts) {
    int iN = ip[0], iC = ip[1], iCA = ip[2], oCB = op[0];
    //    if(verbose(7)) cout <<"IO pts "<< iN <<" "<< iCA <<" "<< iC <<" "<< oCB <<" "<< consts << endl;
    //if(verbose(7)) cout << "VALS " << CA_CB <<" "<< C_CA_CB <<" "<< N_C_CA_CB << endl;
    //find4thPoint( pts[ oCB ],
	//pts[ iN ], pts[ iC ], pts[ iCA ],
	//CA_CB, C_CA_CB, N_C_CA_CB );
    find4thPoint( pts[ oCB ],
        pts[ iC ], pts[ iN ], pts[ iCA ],
        CA_CB, N_CA_CB, C_N_CA_CB );
    //cout << "CONST VALS " << consts <<" "<< CA_CB << " " << C_CA_CB << " " << N_C_CA_CB << endl;
    //    if(verbose(7)) cout << pts[oCB][0] <<" "<< pts[oCB][1] <<" "<< pts[oCB][2] << endl;
    return true;
}

CBbuilder CBbuilder::makeCopy() { return CBbuilder(*this); }
