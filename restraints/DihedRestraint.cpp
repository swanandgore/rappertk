#include <iostream>
using namespace std;

#include <math.h>
#include "DihedRestraint.h"
#include "misc/verbosity.h"
#include "geometry/functions.h"

DihedRestraint::DihedRestraint(vector<int>& pis, const char *desc, vector<float> exp, vector<float> err)
        : Restraint(pis, desc), expected(exp), error(err) {
}

bool DihedRestraint::satisfied(VVF & pts) {
  //    if(verbose(7)) cout << "checking DihedRestraint" << endl;
    float a = calcDihed( pts [ ptInds[0] ], pts [ ptInds[1] ], pts [ ptInds[2] ] , pts [ ptInds[3] ] );
    for(int i=0; i < expected.size(); i++)
        if(dihedDiff(expected[i], a) < error[i]) return true;
    return false;
}

void DihedRestraint::describe() {
    cout << "DihedRestraint on "<< ptInds[0] <<" "<< ptInds[1] <<" "<< ptInds[2] <<" "<< ptInds[3];
}
