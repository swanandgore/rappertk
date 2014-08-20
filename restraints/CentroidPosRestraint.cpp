#include <math.h>
#include <iostream>
using namespace std;

#include "misc/verbosity.h"
#include "CentroidPosRestraint.h"

CentroidPosRestraint::CentroidPosRestraint(vector<int>& pis, const char *desc, vector<float>& c, float r)
        : Restraint(pis, desc) {
    cen = c;
    rad = r;
}

bool CentroidPosRestraint::satisfied(VVF & pts) {
    vector<float> centroid(3,0);
    int index = 0;
    for(int i=0; i < ptInds.size(); i++) {
        index = ptInds[i];
        for(int k=0; k < 3; k++) centroid[k] += pts[index][k];
    }
    for(int k=0; k < 3; k++) centroid[k] /= ptInds.size();

    //    if(verbose(7)) cout << "checking CentroidPosRestraint (" << cen[0] <<","<< cen[1] <<" "<< cen[2] <<") "
    //    << rad << " ("<< centroid[0] <<","<< centroid[1] <<","<< centroid[2] <<")" << endl;
    float dx = fabs(cen[0] - centroid[0]);
    if(dx > rad) return false;
    float dy = fabs(cen[1] - centroid[1]);
    if(dy > rad) return false;
    float dz = fabs(cen[2] - centroid[2]);
    if(dz > rad) return false;
    if(dx*dx + dy*dy + dz*dz > rad*rad) return false;

    //    if(verbose(7)) cout << "ok" << endl;
    return true;
}

void CentroidPosRestraint::describe() {
    cout << "CentroidPosRestraint cen " << cen[0] <<" "<< cen[1] <<" "<< cen[2] <<" rad "<< rad <<" on";
    for(int i=0; i < ptInds.size(); i++) cout << " " << ptInds[i];
}