#include "misc/verbosity.h"
#include "CAtraceGH.h"
#include <math.h>

#include <iostream>
#include<cstring>
#include<cstdlib>

using namespace std;


CAtraceGH::CAtraceGH(vector<int> & scList, vector<int> & positions, float scRed, float ssRed) {
    sc = scList; pos = positions;
    scReduction = scRed;
    ssReduction = ssRed;
    spacegroup = "";
}

// sc[i] is 0,1,2,3 for N,CA,C,O, >=4 for sidechains. =5 for sulphurs.
// sc[i] < 0 for dummy atoms, clashcheck is _not_ performed if any of 2 atoms is a dummy
bool CAtraceGH::clash(int pi1, int pi2, float r1, float r2, vector<float> & p1, vector<float> & p2, bool image) {
  //    if(verbose(7)) cout << "checking overlap-" << image <<" "<< pi1 <<" "<< pi2 << endl;

    if(sc[pi1] < 0 || sc[pi2] < 0) return false; // no dummies please

    int posi1 = pos[pi1], posi2 = pos[pi2];

    if(!image) {
        // same residue atoms dont clash
        if(posi1 == posi2) return false;
    
        // some clashchecks for atoms in adjacent residues can be ignored
        if( abs(posi1-posi2) == 1 ) {
            if(posi1 > posi2 && sc[pi1] == 6 && sc[pi2] == 2) return false; // to take care of PRO CD and prev C clash
            if(posi1 < posi2 && sc[pi1] == 2 && sc[pi2] == 6) return false; // to take care of PRO CD and prev C clash
            int a, b; // a ahead, b behind
            if( posi1-posi2 == 1 ) { a = pi1; b = pi2; }
            else /*( pos[pi2]-pos[pi1] == 1 )*/ { a = pi2; b = pi1; }
            if(sc[a] == 0 || sc[a] == 1) // if ahead is N or CA
                if(sc[b] == 1 || sc[b] == 2 || sc[b] == 3) // if behind is CA,C,O
                    return false;
        }
    }

    float radsum = r1 + r2;

    if(sc[pi1] >= 4 || sc[pi2] >= 4) // sidechain reduction
        radsum = radsum * scReduction;
    if(!image && sc[pi1] == 5 && sc[pi2] == 5) // disulphide reduction
        radsum = radsum - ssReduction * scReduction;

    float dx = fabs( p1[0] - p2[0] );
    if(dx > radsum) return false;
    float dy = fabs( p1[1] - p2[1] );
    if(dy > radsum) return false;
    float dz = fabs( p1[2] - p2[2] );
    if(dz > radsum) return false;

    if(dx*dx + dy*dy + dz*dz > radsum*radsum) return false;

    /*    if(verbose(7)) {
        cout << posi1<<":"<<pi1 <<" "<< posi2<<":"<<pi2 << " Overlap-" <<image <<" "<< r1 <<" "<< r2
        <<" ("<< p1[0]<<" "<<p1[1]<<" "<<p1[2]
        <<") ("<< p2[0]<<" "<<p2[1]<<" "<<p2[2] <<") "
        << dx <<" "<< dy <<" "<< dz <<" ("<< sqrt(dx*dx+dy*dy+dz*dz) <<","<< radsum << ") : clashing" << endl;
    }
    */
    return true; // really clash
}



// returned spacegroup is "" if no info has been set for cell,symmetry
string CAtraceGH::cellsym(float& _a, float& _b, float& _c, float& _alpha, float& _beta, float& _gamma) {
    _a = a; _b = b; _c = c;
    _alpha = alpha; _beta = beta; _gamma = gamma;
    return spacegroup;
}

void CAtraceGH::setCellsym(float _a, float _b, float _c, float _alpha, float _beta, float _gamma, const char* sg) {
    a = _a; b = _b; c = _c;
    alpha = _alpha; beta = _beta; gamma = _gamma;
    spacegroup = string(sg);
}
