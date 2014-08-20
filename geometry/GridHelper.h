#ifndef __GRIDHELPER_H__
#define __GRIDHELPER_H__

#include<vector>
#include<string>
using namespace std;

class GridHelper {
public :
    virtual bool clash(int pi1, int pi2, float r1, float r2, vector<float> & p1, vector<float> & p2, bool image) = 0;
    virtual string cellsym(float& a, float& b, float& c, float& alpha, float& beta, float& gamma) = 0;
    virtual void setCellsym(float a, float b, float c, float alpha, float beta, float gamma, const char * sg) = 0;
};

#endif // __GRIDHELPER_H__
