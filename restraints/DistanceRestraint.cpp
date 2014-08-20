#include <iostream>
using namespace std;

#include <math.h>
#include "DistanceRestraint.h"
#include "misc/verbosity.h"

DistanceRestraint::DistanceRestraint(vector<int>& pis, const char *desc, float mindist, float maxdist)
        : Restraint(pis, desc) {
    minsq = mindist * mindist;
    maxsq = maxdist * maxdist;
}

bool DistanceRestraint::satisfied(VVF & pts) {
  //    if(verbose(7)) cout << "checking DistanceRestraint" << endl;
    vector<float>& p = pts [ ptInds[0] ];
    vector<float>& q = pts [ ptInds[1] ];
    
    float dx = fabs(p[0] - q[0]);
    if(dx > maxsq) return false;
    float dy = fabs(p[1] - q[1]);
    if(dy > maxsq) return false;
    float dz = fabs(p[2] - q[2]);
    if(dz > maxsq) return false;

    float totdist = dx*dx + dy*dy + dz*dz;
    //    if(verbose(7)) cout << minsq <<" <? "<< totdist <<" <? "<< maxsq << endl;
    if(totdist < minsq || totdist > maxsq) return false;

    //    if(verbose(7)) cout << "ok" << endl;
    return true;
}

void DistanceRestraint::describe() {
    cout << "DistanceRestraint between " << sqrt(minsq) <<" "<< sqrt(maxsq) <<" on "<< ptInds[0] <<" "<< ptInds[1];
}
