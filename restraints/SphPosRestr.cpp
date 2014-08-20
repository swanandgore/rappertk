#include <iostream>
using namespace std;

#include <math.h>
#include "SphPosRestr.h"
#include "misc/verbosity.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
SphPosRestr::SphPosRestr(vector<int>& pis, const char *desc, vector<float>& p, float rmin, float rmax)
        : Restraint(pis, desc) {
    pos = p;
    radMin = rmin;
    radMax = rmax;
    if(rmin > rmax || rmin < 0 || rmax < 0) { cout << "Fatal error, SphPosRestr " << desc << " has incorrect radii " << rmin <<" "<< rmax << endl; exit(0); }
    radMinsq = radMin * radMin; //precalculate to avoid squaring later
    radMaxsq = radMax * radMax;
}

bool SphPosRestr::satisfied(VVF & pts) {
    vector<float>& p = pts[ptInds[0]];

    //    if(verbose(7))
    //    cout << "checking SphPosRestr (" << p[0] <<","<< p[1] <<" "<< p[2] <<") "
    //       << radMin<<":"<<radMax << " ("<< pos[0] <<","<< pos[1] <<","<< pos[2] <<") "
    //       << sqrt( (p[0]-pos[0])*(p[0]-pos[0]) + (p[1]-pos[1])*(p[1]-pos[1]) + (p[2]-pos[2])*(p[2]-pos[2]) ) << endl;

    float dx = fabs(p[0] - pos[0]);
    if(dx > radMax) return false;
    float dy = fabs(p[1] - pos[1]);
    if(dy > radMax) return false;
    float dz = fabs(p[2] - pos[2]);
    if(dz > radMax) return false;

    float radsq = dx*dx + dy*dy + dz*dz;
    if(radsq > radMaxsq) return false;
    if(radsq < radMinsq) return false;

    //    if(verbose(7)) cout << "ok" << endl;

    return true;
}

void SphPosRestr::describe() {
    cout << "SphPosRestr cen " << pos[0] <<" "<< pos[1] <<" "<< pos[2] <<" rad-min,max "<< radMin<<":"<<radMax <<" on "<< ptInds[0];
}
