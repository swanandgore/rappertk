#ifndef __RESTRAINT_H__
#define __RESTRAINT_H__

#include <vector>
#include <string>
using namespace std;

typedef vector<vector<float> > VVF;

class Restraint {
public :
    Restraint(vector<int>& pis, const char *desc);
    virtual bool satisfied(VVF & pts);
    virtual ~Restraint();
    vector<int>& getInds();
    virtual void describe();
    const char * name();
protected :
    vector<int> ptInds;
    string description;
};

#endif //__RESTRAINT_H__
