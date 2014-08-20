#ifndef __RATDATA_H__
#define __RATDATA_H__

#include <vector>
#include <string>
#include <map>
using namespace std;

class RAT {
public:
    RAT(int i, int j, int k) : r(i), a(j), t(k) {}
    int r, a, t;
};

struct compareRATs {
    bool operator() (const RAT & r1, const RAT & r2) { // return true if <
        if(r1.r < r2.r) return true;
        else if(r1.r > r2.r) return false;
        if(r1.a < r2.a) return true;
        else if(r1.a > r2.a) return false;
        if(r1.t < r2.t) return true;
        else if(r1.t > r2.t) return false;
        return false;
    }
};


typedef RAT PSO;

typedef map< RAT , vector<PSO>, compareRATs > RATmap;

class RATdata {
public:
    static RATdata & bwdinstance(const char *fn);
    static RATdata & fwdinstance(const char *fn);
    void contents();
    void range(string & resn, float r, float a, float t,
        vector<PSO>::iterator & bi, vector<PSO>::iterator & ei);
private :
    void fillmap(string &resn, RAT &rat);
    RATdata();
    RATdata(const char *filename);

    map< string, RATmap > ratmap;
};

#endif // __RATDATA_H__
