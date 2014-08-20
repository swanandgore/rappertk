#ifndef __ROTATOR_H__
#define __ROTATOR_H__

#include "Builder.h"

class Rotator : public Builder {
public :
    Rotator(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc, float rangeMin, float rangeMax, float stepSize);
    bool build(vector<vector<float> > & pts);
private :
    float min, max, step;
};

#endif//__ROTATOR_H__
