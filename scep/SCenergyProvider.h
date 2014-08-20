#ifndef __SCENERGYPROVIDER_H__
#define __SCENERGYPROVIDER_H__

#include <string>
#include <map>
#include <vector>
using namespace std;

class Builder;
class Grid;

class SCenergyProvider {
public :
    SCenergyProvider(vector<Builder*> & chibuilders, vector<vector<float> > & pts, Grid *gr);
    void addMTZinfo(const char *mtzFN, const char *FOlabel, const char *FClabel, const char *PHIlabel, int MAPtype, float ESmin, float ESmax);
    void addMAPinfo(const char *mtzFN, float ESmin, float ESmax);
    float scoreGivenRotamers();
    float pairEn(int sci1, int sci2, int ri1=-1, int ri2=-1);
    bool interact(int sci1, int sci2);
    float selfEn(int sci, int ri=-1);
    int numRot(int sci); // number of valid rotamers ie which dont clash against mainchain. ri has to be less than this in functions above
    //float getMaxCBDist(int sci); float interCBdist(int sci1, int sci2);
    void collapse(vector<int> SCIs, int A);
    void join(vector<int> & finalAssig);
    void buildAssig(vector<int> & assig);
    void energyRandAssig();
    void resetEnergy();
    float findEsofar(vector<int> & SCIs, int J, vector<int> & assign, bool verbose=0);
    void deeGoldstein(int sci);

    void exhaustiveEnergy();
    void exEn(vector<int> & assig, int pos);

private :
    void calcIntxns(int sci1, int sci2);
    float findEbound(vector<int> & SCIs, int J);
    void assign(vector<int> & SCIs, int J, vector<int> & currAssig);
    void calcIEpairmin(int sci1, int sci2);
    void calcEselfmin(int sci);
    void orderOnEself(int sci);

    void printEpair();

    vector<vector<int> > validRotamers;
    vector<vector<float> > eself; vector<float> eselfmin;
    map<int, map<int, map<int, map<int, float> > > > epair; // i<j
    map<int, map<int, float> > epairmin; // i<j
    map<int, map<int, bool> > ipair; // i<j

    // store the related assignments for artipt's valid rotameric states
    // eg artiAssig[artiSci][artiRI][sci] = ri
    // artiRI varies from 0 till validRotamers[artiSci].size()
    // ri varies from 0 till validRotamers[sci].size()
    vector<int> artiOrder;
    vector< vector<map<int,int> > > artiAssig;
    vector< vector<float> > artiEn;

    Grid *grid;
    vector<Builder*> chibs;
    vector<vector<float> > * pts;

    float Ebest; vector<int> bestAssig;

    string mtzfn, folabel, fclabel, philabel;
    int maptype;
    float esmin, esmax, esmean;
    bool mtzmap;
};

#endif//__SCENERGYPROVIDER_H__
