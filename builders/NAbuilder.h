#ifndef __NABUILDER_H__
#define __NABUILDER_H__

#include "Builder.h"

class NAsampler;
class Constants;

class NAbuilder : public Builder {
public:
    NAbuilder(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc,
        NAsampler* nas, char _base, bool deoxy, char* _pucker);
    bool build(VVF & pts);
    static void buildPhosphate(Constants *consts, vector<float>& o3, vector<float>& p, vector<float>& o5, vector<float>& o1p, vector<float>& o2p);
    static void buildBackbone( Constants *consts, vector<float>& o31, vector<float>& p, vector<float>& o5,
            vector<float>& c5, vector<float>& c4, vector<float>& c3, vector<float>& o3, vector<float>& pn, vector<float>& o5n,
            float alpha, float beta, float gamma, float delta, float epsilon, float zeta );
    static void buildSugar( Constants *consts, string& pucker, vector<float>& c5, vector<float>& c4, vector<float>& c3, vector<float>& o3,
            vector<float>& o4, vector<float>& c1, vector<float>& c2, vector<float>& o2 );
    static void build_AG_scaffold( Constants *consts, string & pucker, char base, vector<float>& c4s, vector<float>& o4s, vector<float>& c1s,
        vector<float>& n1, vector<float>& c2, vector<float>& n3, vector<float>& c4, vector<float>& c5, vector<float>& c6,
        vector<float>& n7, vector<float>& c8, vector<float>& n9 );
    static void build_A( Constants *consts, vector<float>& c2, vector<float>& n1, vector<float>& c6, vector<float>& n6 );
    static void build_G( Constants *consts, vector<float>& c2, vector<float>& n1, vector<float>& c6, vector<float>& o6, vector<float>& n2 );
    static void build_TCU_scaffold( Constants *consts, string & pucker, char base,
        vector<float>& c4s, vector<float>& o4s, vector<float>& c1s, // ips, s for sugar
        vector<float>& n1, vector<float>& c2, vector<float>& o2, vector<float>& n3, vector<float>& c4, vector<float>& c5, vector<float>& c6 );
    static void build_T(Constants *consts, vector<float>& c2, vector<float>& n3, vector<float>& c4, vector<float>& c5, vector<float>& o4, vector<float>& c5m);
    static void build_C(Constants *consts, vector<float>& c2, vector<float>& n3, vector<float>& c4, vector<float>& n4);
    static void build_U(Constants *consts, vector<float>& c2, vector<float>& n3, vector<float>& c4, vector<float>& o4);

private:
    float alpha, beta, gamma, epsilon, zeta; // delta and chi are taken from constants, others are sampled
    char base;
    string pucker;
    NAsampler *naSampler;
    bool deoxy;
};

#endif //__NABUILDER_H__
