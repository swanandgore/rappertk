#ifndef __BUILDERGROUP_H__
#define __BUILDERGROUP_H__

#include "Builder.h"
#include <vector>
using namespace std;

class BuilderGroup : public Builder {
public:
    BuilderGroup(vector<Builder*>& builders, const char* desc);
    bool build(vector<vector<float> > & pts);
    BuilderGroup makeCopy();
private:
    vector<Builder*> builders;
};

#endif // __BUILDERGROUP_H__
