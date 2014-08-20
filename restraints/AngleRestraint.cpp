#include <iostream>
using namespace std;

#include <math.h>
#include "AngleRestraint.h"
#include "misc/verbosity.h"
#include "geometry/functions.h"

AngleRestraint::AngleRestraint(vector<int>& pis, const char *desc, float min, float max)
        : Restraint(pis, desc), minAng(min), maxAng(max) {
}

bool AngleRestraint::satisfied(VVF & pts) {
    float a = calcAngle( pts [ ptInds[0] ], pts [ ptInds[1] ], pts [ ptInds[2] ] );
    if(a < minAng) return false;
    if(a > maxAng) return false;
    return true;
}

void AngleRestraint::describe() {
    cout << "AngleRestraint between " << minAng <<" "<< maxAng <<" on "<< ptInds[0] <<" "<< ptInds[1] <<" "<< ptInds[2];
}
