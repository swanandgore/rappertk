#include "DihedBuilder.h"

#include <assert.h>
#include "geometry/functions.h"
#include "misc/RanGen.h"

DihedBuilder::DihedBuilder(vector<int>& ipInds, vector<int>& opInds, const char *desc,
	float length, float angle, float dihed) : Builder(ipInds, opInds, NULL, desc) {
    if(ipInds.size() != 3 || opInds.size() != 1) { cout <<"point sizes incorrect for DihedBuilder"<< endl; exit(1); }
    len = length;
    ang = angle;
    dih = dihed;
}

bool DihedBuilder::build(VVF & pts) {
    //dih = 360. * ran01() - 180.;
    //cout << "DIH " << dih << endl;
    find4thPoint( pts[ op[0] ],
	    pts[ ip[0] ], pts[ ip[1] ], pts[ ip[2] ],
	    len, ang, dih );
    return true;
}

DihedBuilder DihedBuilder::makeCopy() {
    return DihedBuilder(*this);
}
