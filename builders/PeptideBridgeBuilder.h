#ifndef __PEPTIDEBRIDGEBUILDER_H__
#define __PEPTIDEBRIDGEBUILDER_H__

#include "Builder.h"

#include <set>
#include <map>
#include <string>
using namespace std;

class PeptideBridgeBuilder : public Builder {
public :
    PeptideBridgeBuilder(vector<int>& ipInds, vector<int>& opInds, Constants* con, const char* desc, const char *r0, const char *r1, const char *r2, const char *rdpath);
    bool build(VVF & pts);
    PeptideBridgeBuilder makeCopy();
private :
    string RATdataPath, resn0, resn1, resn2;
    bool checkAndAddIfNew(float phi1, float psi1, float phi2, float psi2, float d1, float d2);
    bool checkAndAddIfNew(float d1, float d2);
    void deleteSessionInfo();
    float CA_C, N_CA_C, C_N, CA_C_N, C_O, CA_C_O, TAU_QUALITY;
    set <long long int> keys; map<int,int> ctkeys;

    void buildCBs(VVF & pts);
public :
    static void setCTtrials(int ct);
    static void setThetastep(int ts);
    static int cttrials, thetastep;
    float transProp0, transProp1;
};

#endif // __PEPTIDEBRIDGEBUILDER_H__
