#ifndef __BUILDER_H__
#define __BUILDER_H__

/**
 * all builders will be derived from this baseclass.
 */

#include <vector>
#include <string>
using namespace std;

typedef vector<vector<float> > VVF;

#define cget(X) (consts->get((X)))

class Constants;

class Builder {
public:
    Builder(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc);
    Builder();
    virtual void describe();
    virtual ~Builder();
    virtual bool build(VVF & pts);
    virtual bool buildSample(VVF & pts, int sampleIndex);
    vector<int>& getIP();
    vector<int>& getOP();
    const char* name();

    void startSession(int ssize);
    void endSession();
protected:
    vector<int> ip, op;
    Constants *consts;
    string description;

    bool session; // is session on? if yes, remember and verify samples, else not
    int maxSessionSize, sessionSize; // dont sample more than this number of unique samples
    virtual void deleteSessionInfo();
};

#endif // __BUILDER_H__
