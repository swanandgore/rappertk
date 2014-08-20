#ifndef __SCLSAMPLE_H__
#define __SCLSAMPLE_H__

#include <vector>
#include <map>
using namespace std;

class SCLsampler {
public :
    SCLsampler(const char *crdfilename, const char *propfilename, char aa);
    int sample(int phi, int psi, vector<vector<float> > ** pts); // sample and copy the into ptinds in pts
    int sample(int phi, int psi, vector<vector<float> > ** pts, int sampleIndex); // copy the specified sample into ptinds in pts
    SCLsampler makeCopy();
    void addSample(float prop, vector<vector<float> > newrot);

private :
    void copyRotamer(int roti, vector<vector<float> > & pts, vector<int> & ptinds);
    void findNbrOccuPhipsi(int phi, int psi, int & newphi, int & newpsi);
    void properPhipsi(int & phi, int & psi);

    vector<vector<vector<float> > > rotamers; // rotamer coordinates
    map<int,map<int,vector<int> > > jump; // jumptable from phi-psi-random to index in ps2rot[phi][psi]. ps2rot[phi][psi][index] is right rotamer-index in rotamers
    map<int,map<int,vector<int> > > ps2rot; // 40 degree bins, 40 degree bins, rotamer indices
    static int MAGFAC, BINSIZE;
};

#endif // __SCLSAMPLE_H__
