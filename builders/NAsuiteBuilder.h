#ifndef __NASUITEBUILDER_H__
#define __NASUITEBUILDER_H__

#include "Builder.h"
#include "samplers/NAsuiteSampler.h"

class Constants;

class NAsuiteBuilder : public Builder {
public:
    NAsuiteBuilder(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc,
        NAsuiteSampler* nas, char _base, bool _deoxy, char* _pucker, char* _puckerNext);
    bool build(VVF & pts);
    NAsuiteBuilder makeCopy();
private:
    NAsuiteSampler *naSampler;
    string pucker, puckerNext;
    char base;
    bool deoxy;
    static void buildBackbone( Constants *consts, vector<float>& c5s, vector<float>& c4s, vector<float>& c3s,
        vector<float>& o3s, vector<float>& pn, vector<float>& o5sn, vector<float>& c5sn, vector<float>& c4sn, vector<float>& c3sn,
        float delta, float epsilon, float zeta, float alpha, float beta, float gamma );

    vector<SuiteVal> sessionSVs;
    bool checkAndAddIfNew(SuiteVal & sv);
    void deleteSessionInfo();
};

#endif // __NASUITEBUILDER_H__
