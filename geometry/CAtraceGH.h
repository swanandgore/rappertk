#ifndef __CATRACEGH_H__
#define __CATRACEGH_H__

#include "GridHelper.h"
#include <string>
#include <set>
#include <map>
using namespace std;

class CAtraceGH : public GridHelper {
public :
    // positions are residue positions, adjacent positions are considered connected
    // so for chain breaks, put arbitrary gap (say 10)
    CAtraceGH(vector<int> & scList, vector<int> & positions, float scRed, float ssRed);
    bool clash(int pi1, int pi2, float r1, float r2, vector<float> & p1, vector<float> & p2, bool image);
    string cellsym(float& a, float& b, float& c, float& alpha, float& beta, float& gamma);
    void setCellsym(float a, float b, float c, float alpha, float beta, float gamma, const char* sg);
    virtual ~CAtraceGH() {}
private :
    // residue-position of each atom in the chain. for atoms are in adjacent positions, some clashcheck are not done!!
    vector<int> pos;
    // which atoms are in sidechain? vdw rad will be lowered for this
    // sc[i] is 0,1,2,3 for N,CA,C,O. 4 or more for sidechain atoms. 5 for sulphurs that may form disulphides
    vector<int> sc;
    float scReduction, ssReduction;

    float a,b,c,alpha,beta,gamma; string spacegroup;
};

#endif // __CATRACEGH_H__
