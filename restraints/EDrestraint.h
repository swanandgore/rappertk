#ifndef __EDRESTRAINT_H__
#define __EDRESTRAINT_H__

#include "Restraint.h"

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <string>
using std::string;


class EDrestraint : public Restraint {
public :
    static EDrestraint makeEDrestraintFromMTZ(vector<int>& pis, const char *desc, const char *mtzfn, const char *folabel, const char *fclabel, const char *philabel,
            int maptype, float zero, float max, float cutoff);
    static EDrestraint makeEDrestraintFromMap(vector<int>& pis, const char *desc, const char *mapfn, float zero, float max, float cutoff);
    static float findBadfit(const char *mtzfn, char *folabel, char *fclabel, char *philabel, vector<vector<float> > & pts, vector<int> & ptinds);
    static float findBadfit1(const char *mapfn, const char *mapfn2, vector<vector<float> > & pts, vector<int> & ptinds);
    bool satisfied(VVF & pts);
    float scoreSigma(VVF & pts);
    float scoreAround(VVF & pts);
    float evalDensity(float radius);
    float score(VVF & pts);
    void describe();
    EDrestraint(vector<int> & pis, const char *desc);
private :
    void *xmap; // clipper's xmap properly initialized
    float densityCutoff; // minimum density above which the density is considered 'good'
    float zero, max, cutoff; // maximum density below which the density is considered bad
    float mean, stdev; // of density over asymmetric unit

public:
    static float getPenalty(); // for penalty constant to scale up negative density occupation
    static void setPenalty(float np);
    static float negpen;
    static float scatterRad; // for scattere radius. used for scoreAround() and calculating correlation coefficients
    static float getScatRad();
    static void setScatRad(float sr);
};


class SingleMap {
public :
    static SingleMap* instance();
    void * map(const char *mtzfn, const char *folabel, const char *fclabel, const char *philabel, int maptype, float& mean, float& stdev);
    void * map(const char *mapfilename, float & mean, float & stdev);
    void addMap(string mtz, void *xmap);
private :
    SingleMap() {}
    std::map< string, void* > xmaps;
    std::map< string, float > mapMeans;
    std::map< string, float > mapStdevs;
};
#endif//__EDRESTRAINT_H__
