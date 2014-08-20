#include "CNCaBuilder.h"
#include "Constants.h"
#include "geometry/functions.h"
#include "misc/RanGen.h"

CNCaBuilder::CNCaBuilder(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc,
        vector<float> & capos, float caposTol)
        : Builder(ipInds, opInds, con, desc) {
    CApos = capos;
    CAposTol = caposTol;
}

// given CA, build N and C before it
// choose point N at N-CA distance and then choose C using C-N-CA and C-N
bool CNCaBuilder::build(VVF & pts) {
    vector<float> v(3), a(3), b(3);
    randomNormalVector(v);
    findYZgivenX(v, a, b);

    vector<float> & c = pts[ op[0] ];
    vector<float> & n = pts[ op[1] ];

    vector<float> & ca = pts[ op[2] ]; // add noise to ca position around CApos
    float noise = ran01() * CAposTol;
    vector<float> noi(3); randomNormalVector(noi);
    //cout << "CNCA " << CApos[0] <<" "<< CApos[1] <<" "<< CApos[2] << endl;
    for(int i=0; i < 3; i++) ca[i] = CApos[i] + noi[i]*noise;
    //cout << "CNCA " << ca[0] <<" "<< ca[1] <<" "<< ca[2] << endl;

    float d = cget("N_CA");
    for(int i=0; i < 3; i++) n[i] = ca[i] + v[i] * d;

    float ang = DEGREES_TO_RADIANS(180. - cget("C_N_CA")); d = cget("C_N");
    float d1 = d * cos(ang), d2 = d * sin(ang);
    for(int i=0; i < 3; i++) c[i] = n[i] + d1 * v[i] + d2 * a[i];

    return true;
}
