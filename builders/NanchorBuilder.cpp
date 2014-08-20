#include <assert.h>
#include <math.h>

#include "misc/RanGen.h"
#include "misc/verbosity.h"
#include "geometry/functions.h"
#include "NanchorBuilder.h"
#include "Constants.h"

NanchorBuilder::NanchorBuilder(vector<int> & ipInds, vector<int> & opInds, Constants *con, const char *desc,
    vector<float>& CA0, float R0, vector<float>& CA1, float R1) 
        : Builder(ipInds, opInds, con, desc)
{
    ca0 = CA0; ca1 = CA1;
    r0 = R0; r1 = R1;
}

// inputs : none, outputs : CA,C,O,N,CA
// build N Ca C O N Ca 'ab-initio'.
// first choose cis or trans, and accordingly ca0, ca1.
// then choose a plane which will contain ca0, c0, n1, ca1 and place c0, n1
// randomly choose a psi and and build n0
bool NanchorBuilder::build(VVF & pts) {
    // is trans possible for ca0-ca1?
    bool trans = false;
    if( calcDist(ca0,ca1) + r0 + r1 > 3.85 ) trans = true;

    float dist = 2.81;
    if(trans) dist = 3.81;

    int CA0 = op[0], C0 = op[1], O0 = op[2], N1 = op[3], CA1 = op[4];
    vector<float> v(3,0), dir(3,0); float noi;
    while(1) { // find ca0, ca1
        randomNormalVector(v); //---
        noi = r0 * ran01();
        linear_combination(pts[CA0], 1, ca0, noi, v);
        randomNormalVector(v); //---
        noi = r1 * ran01();
        linear_combination(pts[CA1], 1, ca1, noi, v);
        linear_combination(dir, 1, pts[CA1], -1, pts[CA0]);
        normalize_point(dir,dir);
        linear_combination(pts[CA1], 1, pts[CA0], dist, dir);
        if( calcDist(pts[CA1],ca1) > r1 ) continue;
        break;
    }

    //cout << "trans " << trans <<" "<< dist << endl;
    //cout << "CA0 "<< pts[CA0][0]<<" "<<pts[CA0][1]<<" "<<pts[CA0][2] << endl;
    //cout << "CA1 "<< pts[CA1][0]<<" "<<pts[CA1][1]<<" "<<pts[CA1][2] << endl;

    vector<float> caAxis(3,0), p1(3,0), p2(3,0);
    linear_combination(caAxis, 1, pts[CA1], -1, pts[CA0]);
    normalize_point(caAxis, caAxis);
    findYZgivenX(caAxis, p1, p2);

    float c_n = cget("C_N");
    float n_ca = cget("N_CA");
    float theta = DEGREES_TO_RADIANS(cget("C_N_CA"));
    float c_ca = sqrt( c_n*c_n + n_ca*n_ca - 2*c_n*n_ca*cos(theta) ); // cos rule
    //cout << "C_CA " << c_ca << endl;

    float ca_c = cget("CA_C");
    theta = DEGREES_TO_RADIANS(cget("CA_C_N"));
    float ca_n = sqrt( ca_c*ca_c + c_n*c_n - 2*ca_c*c_n*cos(theta) ); // cos rule
    //cout << "CA_N " << ca_n << endl;

    float base, perp;
    findTriangleBaseHeight(c_ca, ca_c, dist, base, perp);
    //cout << "base-perp " << base <<" "<< perp << endl;
    linear_combination(pts[C0], 1, pts[CA0], base, caAxis);
    linear_combination(pts[C0], 1, pts[C0], perp, p1); // C

    findTriangleBaseHeight(n_ca, ca_n, dist, base, perp);
    //cout << "base-perp " << base <<" "<< perp << endl;
    linear_combination(pts[N1], 1, pts[CA0], base, caAxis);
    if(trans) linear_combination(pts[N1], 1, pts[N1], -1*perp, p1);
    else linear_combination(pts[N1], 1, pts[N1], perp, p1); // N

    find4thPoint(pts[O0], pts[N1], pts[CA0], pts[C0], cget("C_O"), cget("CA_C_O"), -180); // O

    assert( calcDist(ca0,pts[CA0]) < r0 );
    //    if(verbose(7)) cout << "ca0pts1 " << ca0[0]<<" "<<ca0[1]<<" "<<ca0[2]<<" - "<<pts[CA0][0]<<" "<<pts[CA0][1]<<" "<<pts[CA0][2] << endl;
    assert( calcDist(ca1,pts[CA1]) < r1 );

    return true;
}

NanchorBuilder NanchorBuilder::makeCopy() { return NanchorBuilder(*this); }
