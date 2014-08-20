#ifndef __TRANSROTATOR_H__
#define __TRANSROTATOR_H__

#include "Builder.h"

// inputs are none, outputs are points A, B, P1,P2, ...
// New positions for A,B are sampled from their given spheres, so that d(A,B) dsnt change
// A,B are transformed to A'B', and same transform is appld to rest of the points
//  note that this builder dsnt change any lengths, angles or dihedrals
class TransRotator : public Builder {
public :
    TransRotator(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc, vector<float> c1, float r1, vector<float> c2, float r2);
    bool build(vector<vector<float> > & pts);
private :
    vector<float> cen1, cen2;
    float rad1, rad2;
};

#endif//__TRANSROTATOR_H__
