#include "SCLsampler.h"
#include "misc/RanGen.h"
#include "geometry/functions.h"
#include "misc/verbosity.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
using namespace std;

int SCLsampler::MAGFAC = 1000;
int SCLsampler::BINSIZE = 40;

SCLsampler::SCLsampler(const char *crdfilename, const char *propfilename, char aa) {
  //    if(verbose(6)) cout << "SCLsampler " << crdfilename <<" "<< propfilename <<" "<< aa << endl;
    char line[500];
    map<int,map<int,vector<float> > > ps2prop; // BINSIZE degree bins, BINSIZE degree bins, propensity for rotamer indices

    for(int phi=-180; phi < 180; phi+=BINSIZE) {
        ps2rot[phi] = map<int, vector<int> >();
        ps2prop[phi] = map<int, vector<float> >();
        for(int psi=-180; psi < 180; psi+=BINSIZE) {
            ps2rot[phi][psi] = vector<int>();
            ps2prop[phi][psi] = vector<float>();
        }
    }

    bool modelOn = false;
    FILE *fp = fopen(crdfilename, "r");
    while(fgets(line, 500, fp) != NULL) {
        if(strlen(line) > 21 && line[21] != aa) continue;
        if(strlen(line) > 76 && line[76] == 'H') continue;
        if(strlen(line) > 13 && line[13] == 'H') continue;
        if(strstr(line,"ATOM") != line) continue;
        if(strstr(line, " N  ") == &line[12]) rotamers.push_back(vector<vector<float> >());
        //if(!modelOn && strstr(line,"MODEL") == line) { modelOn = true; rotamers.push_back(vector<vector<float> >()) ; continue; }
        //if(strstr(line,"ENDMDL") == line) { modelOn = false; continue; }
        int nrot = rotamers.size()-1;
        rotamers[nrot].push_back(vector<float>(3));
        int ncrd = rotamers[nrot].size()-1;
        sscanf(&line[30], "%8f", &rotamers[nrot][ncrd][0]);
        sscanf(&line[38], "%8f", &rotamers[nrot][ncrd][1]);
        sscanf(&line[46], "%8f", &rotamers[nrot][ncrd][2]);
    }
    fclose(fp);
    for(int i=0; i < rotamers.size(); i++) {
      //        if(verbose(6)) cout << "rotsize "<< aa <<" "<< i <<" "<< rotamers[i].size() << endl;
        if(rotamers[i].size() != rotamers[0].size()) { cout << "rotamer[i] and [0] sizes dont match!" << endl; exit(0); }
    }

    fp = fopen(propfilename, "r");
    int rotind, phi, psi; float prop, fphi, fpsi;
    while(fgets(line, 500, fp) != NULL) {
        if(line[0] != aa || line[1] != ' ') continue;
        sscanf(&line[2], "%d %f %f %f", &rotind, &fphi, &fpsi, &prop);
        phi = fphi ; psi = fpsi;
        if(prop <= 0.0001) continue; //prop = 0.0001;
        ps2rot[phi][psi].push_back(rotind);
        ps2prop[phi][psi].push_back(prop);
    }
    fclose(fp);

    //populate jumptables
    for(int phi=-180; phi < 180; phi+=BINSIZE) {
        jump[phi] = map<int, vector<int> >();
        for(int psi=-180; psi < 180; psi+=BINSIZE)
            jump[phi][psi] = vector<int>(MAGFAC, -1);
    }

    for(int phi=-180; phi < 180; phi+=BINSIZE) {
        for(int psi=-180; psi < 180; psi+=BINSIZE) {
            if(ps2prop[phi][psi].size() == 0) continue;
            for(int i=1; i < ps2prop[phi][psi].size(); i++)
                ps2prop[phi][psi][i] += ps2prop[phi][psi][i-1];
            for(int i=0; i < ps2prop[phi][psi].size(); i++) {
                ps2prop[phi][psi][i] /= ps2prop[phi][psi][ ps2prop[phi][psi].size()-1 ];
                ps2prop[phi][psi][i] = round(ps2prop[phi][psi][i] * MAGFAC);
            }
            //cout << phi <<" "<< psi <<" "<< ps2rot[phi][psi].size() << endl;
            //cout << ps2prop[phi][psi][0] <<" "<< ps2prop[phi][psi][ ps2prop[phi][psi].size()-1 ] << endl;
            if( ps2prop[phi][psi][0] <= 0 || ps2prop[phi][psi][ ps2prop[phi][psi].size()-1 ] != MAGFAC )
            { cout <<"count error "<< ps2prop[phi][psi][0] <<" "<< ps2prop[phi][psi][ ps2prop[phi][psi].size()-1 ] << endl; exit(0); }
            for(int i=0; i < ps2prop[phi][psi].size(); i++) {
                int start = 0, end = ps2prop[phi][psi][i];
                if(i > 0) start = ps2prop[phi][psi][i-1];
                for(int k=start; k < end; k++) jump[phi][psi][k] = i; //ps2rot[phi][psi][i];
            }
        }
    }
    //cout << "done sampler ctor" << endl;
}

// for phi-psi bin with no rotamers, find closest occupied one, assuming BINSIZE degree binning
void SCLsampler::findNbrOccuPhipsi(int phi, int psi, int & newphi, int & newpsi) {
    if(ps2rot[phi][psi].size() > 0) { newphi = phi; newpsi = psi; return;}

    for(int ni = 1; ni < 180/BINSIZE; ni ++) { // levels of nbrhoods; shd be enough i guess?
        int low = -1*ni*BINSIZE, high = ni*BINSIZE, step=ni*BINSIZE;
        for(int dphi=low; dphi <= high; dphi+=step) {
            newphi = withinPlusMinus180 (phi + dphi);
            for(int dpsi=low; dpsi <= high; dpsi+=step) {
                newpsi = withinPlusMinus180 (psi + dpsi);
                if(ps2rot[newphi][newpsi].size() > 0) return;
            }
        }
    }
    newphi = phi; newpsi = psi;
}

void SCLsampler::copyRotamer(int roti, vector<vector<float> > & pts, vector<int> & ptinds) {
    for(int pi=0; pi < rotamers[roti].size(); pi++) {
        int di = ptinds[pi];
        pts[di] = rotamers[roti][pi];
    }
}

void SCLsampler::properPhipsi(int & phi, int & psi) {
    //cout << phi <<" given "<< psi << endl;
    int tmpphi = floor( (phi+180.)/BINSIZE ) * BINSIZE - 180;
    int tmppsi = floor( (psi+180.)/BINSIZE ) * BINSIZE - 180;
    findNbrOccuPhipsi(tmpphi, tmppsi, phi, psi);
    //cout << phi <<" new   "<< psi << endl;
}

int SCLsampler::sample(int phi, int psi, vector<vector<float> > ** rotpts) {
    properPhipsi(phi, psi);

    int ri = int(floor( jump[phi][psi].size() * ran01() ));
    int rrii = jump[phi][psi][ri];
    int roti = ps2rot[phi][psi][rrii];
    //cout << "here " << ri <<" "<< rrii <<" "<< ps2rot[phi][psi].size() <<" "<< roti <<" "<< rotamers.size() << endl;

    (*rotpts) = & rotamers[roti];
    return rrii;
}

int SCLsampler::sample(int phi, int psi, vector<vector<float> > ** rotpts, int sampleIndex) {
    properPhipsi(phi, psi);

    if(sampleIndex >= ps2rot[phi][psi].size()) return -1;
    int roti = ps2rot[phi][psi][sampleIndex];
    if(rotpts != NULL) (*rotpts) = & rotamers[roti]; // essential for access from findMinNumSamples in prepareChain
    return sampleIndex;
}

SCLsampler SCLsampler::makeCopy() { return SCLsampler(*this); }

void SCLsampler::addSample(float prop, vector<vector<float> > newrot) {
    rotamers.push_back(newrot);
    if(rotamers.size() > 1 && rotamers[0].size() != newrot.size()) {
        cout << "fatal SCL::addSample error " << rotamers[0].size() <<" "<< newrot.size() << endl; exit(1);
    }
    //for(int i=0; i < newrot.size(); i++) {
        //for(int k=0; k < newrot[i].size(); k++)
            //cout << newrot[i][k] <<" ";
        //cout << endl;
    //}
    int ni = rotamers.size()-1;
    int phi, psi;
    for(map<int,map<int,vector<int> > >::iterator phit = ps2rot.begin(); phit != ps2rot.end(); ++phit) {
        phi = phit->first;
        for(map<int,vector<int> >::iterator psit = ps2rot[phi].begin(); psit != ps2rot[phi].end(); ++psit) {
            psi = psit->first;
            if(ps2rot[phi][psi].size() == 0) continue;
            ps2rot[phi][psi].push_back(ni);
            int newsteps = int(floor(jump[phi][psi].size() * prop));
            for(int i=0; i < newsteps; i++) jump[phi][psi].push_back(ps2rot[phi][psi].size()-1);
        }
    }
}


////{{{
//int SCLsampler::sample(int phi, int psi, vector<vector<float> > & pts, vector<int> & ptinds) {
//    properPhipsi(phi, psi);
//
//    int ri = int(floor( MAGFAC * ran01() ));
//    int rrii = jump[phi][psi][ri];
//    int roti = ps2rot[phi][psi][rrii];
//    copyRotamer(roti, pts, ptinds);
//    return rrii;
//}
//int SCLsampler::sample(int phi, int psi, vector<vector<float> > & pts, vector<int> & ptinds, int sampleIndex) {
//    properPhipsi(phi, psi);
//
//    if(sampleIndex >= ps2rot[phi][psi].size()) return -1;
//    int roti = ps2rot[phi][psi][sampleIndex];
//    copyRotamer(roti, pts, ptinds);
//    return sampleIndex;
//}
////}}}
