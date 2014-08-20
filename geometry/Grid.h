#ifndef __GRID_H__
#define __GRID_H__

#include <iostream>
#include <set>
#include <vector>
#include <map>
using namespace std;

class VoxelIndex {
public :
    VoxelIndex() : i(-999), j(-999), k(-999) {}
    VoxelIndex(int a, int b, int c) : i(a), j(b), k(c) {}
    int i, j, k;
};

struct compareVoxelIndices { // return true if strictly less
    bool operator()(const VoxelIndex &vi1, const VoxelIndex &vi2) {
	bool retval;
	if(vi1.i < vi2.i) retval = true;
	else if(vi1.i > vi2.i) retval = false;
	else {
	    if(vi1.j < vi2.j) retval = true;
	    else if(vi1.j > vi2.j) retval = false;
	    else {
		if(vi1.k < vi2.k) retval = true;
		else if(vi1.k > vi2.k) retval = false;
		else retval = false;
	    }
	}
	//cout << "comparing "
	//    << vi1.i <<" "<< vi1.j <<" "<< vi1.k <<" : "
	//    << vi2.i <<" "<< vi2.j <<" "<< vi2.k <<" : "
	//    << retval << endl;
	return retval;
    }
};
typedef map< VoxelIndex, set<int>, compareVoxelIndices > VoxelMap;

class MapIntIntFloat { // written as substitute to map<int, map<int,float> >
public:
    MapIntIntFloat();
    void put(int a, int b, float f);
    float get(int a, int b);
    bool contains(int a, int b);
private:
    map<int, map<int,float> > miif;
};
typedef MapIntIntFloat OverlapMap;

class GridHelper;

class Grid {
public :
    Grid(int resol, vector<vector<float> >* pts, vector<float>* rads, GridHelper *gh);
    ~Grid();
    bool add(vector<int>& ptindices);
    void remove(int ptindex);
    void removes(vector<int> & ptindices);
    static void xyz2ijk(float x, float y, float z, int& vi, int& vj, int& vk, int resol);
    bool justAdd(vector<int>& ptinds);
    void printIJK(vector<float>& pt);
    bool withinEnvelope(vector<float> & pt, float rad);

    void popCopyGridPoints(Grid & g); // for use in population-based conf.sampling
    void describe();
private :
    bool voxelClashcheck(int pi, VoxelIndex& vi);
    void printVox(VoxelIndex& vi);
    void justAdd(VoxelIndex& vi, int pi);

    vector<vector<float> >* points;
    vector<float>* radii;
    GridHelper *gh;
    map<int, VoxelIndex> pi2vi; //point index -> voxel index
    VoxelMap vi2pi; // voxel index -> point-indices it contains
    float maxrad; //maximum radius observed so far
    int resolution;
    char *sg; char *cell;
    float specialCutoff, farCutoff;
    map<int, vector<VoxelIndex> > friends; // pi -> voxels with pi's friends
    map<VoxelIndex, map<int,vector<float> >, compareVoxelIndices > imageMap; // voxel,pi -> pi's friends's coordinate in that voxel
    vector<float> llb, rtf;
};

#endif // __GRID_H__

#ifdef DONT_DEFINE_ME
bool Grid::overlappingSpheres(int pi1, int pi2) {
    float r1 = (*radii)[pi1], r2 = (*radii)[pi2];

    float totdist = r1 + r2; // find min dist betn the points with reductions if any
    if(overlapReductions->contains(pi1,pi2))
        totdist -= overlapReductions->get(pi1,pi2);
    if(totdist < 0) totdist = 0;

    bool retval = true;
    float dx, dy, dz;
    vector<float>& p1 = (*points)[pi1]; vector<float>& p2 = (*points)[pi2];
    dx = fabs( p1[0] - p2[0] );
    if(dx > totdist) retval = false;
    else {
        dy = fabs( p1[1] - p2[1] );
        if(dy > totdist) retval = false;
        else {
            dz = fabs( p1[2] - p2[2] );
            if(dz > totdist) retval = false;
            else if( dx*dx+dy*dy+dz*dz >= totdist*totdist ) retval = false;
        }
    }
    if(retval && verbose(7)) {
        cout << pi1 <<" "<< pi2 << " Overlap " << r1 <<" "<< r2
        <<" ("<< p1[0]<<" "<<p1[1]<<" "<<p1[2]
        <<") ("<< p2[0]<<" "<<p2[1]<<" "<<p2[2] <<") "
        << dx <<" "<< dy <<" "<< dz <<" ("<< dx*dx+dy*dy+dz*dz <<","<< totdist*totdist << ") : "<< retval << endl;
    }
    return retval;
}
#endif // DONT_DEFINE_ME
