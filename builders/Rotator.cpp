#include "Rotator.h"
#include "geometry/functions.h"
#include "misc/RanGen.h"

#include <vector>
using namespace std;

Rotator::Rotator(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc, float rangeMin, float rangeMax, float stepSize)
    : Builder(ipInds, opInds, con, desc), min(rangeMin), max(rangeMax), step(stepSize) {}

bool Rotator::build(vector<vector<float> > & pts) {
    // sample angle to rotate thru
    float angle = floor((max-min+0.)/step * ran01()) * step + min;
    // rotate op pts around axis formed by ip pts
    vector<float> axis(3);
    linear_combination(axis, 1, pts[ip[0]], -1, pts[ip[1]]);
    normalize_point(axis, axis);
    vector<vector<float> > rotOp;
    findRotnOperator(rotOp, axis, angle);
    for(int oi=0; oi < op.size(); oi++) {
        linear_combination(pts[op[oi]], 1, pts[op[oi]], -1, pts[ip[0]]);
        rotate(rotOp, pts[op[oi]]);
        linear_combination(pts[op[oi]], 1, pts[op[oi]], +1, pts[ip[0]]);
    }
    return true;
}
