#ifndef __OPTFRAGPLACER_H__
#define __OPTFRAGPLACER_H__

#include "Builder.h"

#include <vector>
using namespace std;

class FragSampler;

class OptFragPlacer : public Builder {
public :
    OptFragPlacer(vector<int>& ipInds, vector<int>& opInds, FragSampler *fsamp, vector<int> & fpis, vector<int> & fpos, const char* desc);
    bool build(VVF & pts);
    OptFragPlacer makeCopy();
protected:
    FragSampler *fs;
    vector<int> fragips, fragops;
};

#endif//__OPTFRAGPLACER_H__
