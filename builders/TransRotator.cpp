#include "TransRotator.h"
#include "geometry/functions.h"

TransRotator::TransRotator(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc, vector<float> c1, float r1, vector<float> c2, float r2)
    : Builder(ipInds, opInds, con, desc), cen1(c1), cen2(c2), rad1(r1), rad2(r2) {}

bool TransRotator::build(VVF & pts) {
    // sample 2 points in the spheres, w/o changing distance AB
    // assume AB to be 1st 2 pts in op points
    bool done = false;
    vector<float> a(3), b(3), ab(3);
    float diff, D = calcDist(pts[op[0]],pts[op[1]]);
    for(int natt=0; natt < 50; natt++) {
        sphereVolSample(cen1, rad1, a);
        sphereVolSample(cen2, rad2, b);
        diff = (D - calcDist(a,b)) / 2.;
        linear_combination(ab, 1, b, -1, a);
        normalize_point(ab, ab);
        linear_combination(a, 1, a, -1*diff, ab);
        linear_combination(b, 1, b, +1*diff, ab);
        if(calcDist(cen1,a) < rad1 && calcDist(cen2,b) < rad2) { done=true; break; }
    }
    if(!done) return false;
    // now translate points s.t. A,A' coincide
    vector<float> aa(3);
    linear_combination(aa, 1, a, -1, pts[op[0]]);
    for(int oi = 0; oi < op.size(); oi++) linear_combination(pts[op[oi]], 1, aa, 1, pts[op[oi]]);
    // now find operator that rotates AB to AB" if A were origin
    vector<float> oldAB(3), rotAxis(3);
    vector<vector<float> > rotOp;
    linear_combination(oldAB, 1, pts[op[1]], -1, pts[op[0]]);
    normalize_point(oldAB, oldAB);
    //cross_product(rotAxis, ab, oldAB); XXX
    cross_product(rotAxis, oldAB, ab);
    findRotnOperator(rotOp, rotAxis, calcAngle(ab,oldAB));
    // translate each point by -A, apply rotation operator, translate back by A
    for(int oi = 1; oi < op.size(); oi++) {
        linear_combination(pts[op[oi]], -1, pts[op[0]], 1, pts[op[oi]]);
        rotate(rotOp, pts[op[oi]]);
        linear_combination(pts[op[oi]], +1, pts[op[0]], 1, pts[op[oi]]);
    }
    return true;
}
