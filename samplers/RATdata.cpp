#include <math.h>

#include <iostream>
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
using namespace std;

#include "misc/verbosity.h"
#include "RATdata.h"

void RATdata::fillmap(string &resn, RAT &rat) {
    if(ratmap.find(resn) == ratmap.end())
        ratmap[resn] = RATmap();
    if(ratmap[resn].find(rat) == ratmap[resn].end())
        ratmap[resn][rat] = vector<PSO>();
}

RATdata::RATdata(const char *filename) {
    char line[500], resn[10];
    int phi, psi, omega, r, a, t;
    FILE *fp = fopen(filename, "r");
    while(fgets(line, 500, fp) != NULL) {
        sscanf(line, "%s %d %d %d %d %d %d", resn, &phi, &psi, &omega, &r, &a, &t);
        string sresn(resn);
        //cout << "read " << resn <<" "<< phi <<" "<< psi <<" "<< omega <<" "<< r <<" "<< a <<" "<< t << endl;
        RAT rat(r,a,t);
        fillmap(sresn, rat);
        ratmap[sresn][rat].push_back( PSO(phi,psi,omega) );
    }
}

void RATdata::contents() {
    cout << "reporting " << ratmap.size() << endl;
    for(map<string,RATmap>::iterator it = ratmap.begin(); it != ratmap.end(); ++it) {
        cout << it->first << " : " << it->second.size() << endl;
        RATmap & rmap = it->second;
        for(RATmap::iterator rit = rmap.begin(); rit != rmap.end(); ++rit) {
            const RAT & r = rit->first;
            vector<PSO> & vp = rit->second;
            cout << "\t" << r.r <<" "<< r.a <<" "<< r.t <<" : "<< vp.size() << endl;
            for(int i=0; i < vp.size(); i++)
                cout << "\t\t\t" << vp[i].r <<" "<< vp[i].a <<" "<< vp[i].t << endl;
        }
    }
}

// assume that dihedrals are within -180 and 180
int dihdiff(int a, int b) {
    int d = abs(a-b);
    if(d > 180) d = 360 - d;
    return d;
}


void RATdata::range(string & resn, float r, float a, float t,
        vector<PSO>::iterator & bi, vector<PSO>::iterator & ei) {

    int ri = (int) round(r*10);
    int ai = (int) round(a/5.) * 5;
    int ti = (int) round(t/5.) * 5;
    RAT rat(ri, ai, ti);
    fillmap(resn, rat);
    bi = ratmap[resn][rat].begin();
    ei = ratmap[resn][rat].end();
    //    if(verbose(7)) cout << "RAT " << resn <<" "<< ri <<" "<< ai <<" "<< ti <<" = "<< ei-bi << " phipsiomega" << endl;
}

RATdata & RATdata::bwdinstance(const char *fn) {
    static RATdata ratdata1(fn);
    return ratdata1;
}
RATdata & RATdata::fwdinstance(const char *fn) {
    static RATdata ratdata2(fn);
    return ratdata2;
}
