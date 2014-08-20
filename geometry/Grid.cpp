#include "Grid.h"
#include "GridHelper.h"
#include "misc/verbosity.h"

#include <assert.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

#ifdef CLIPPERGEOM
#include "clipper/core/coords.h"
using clipper::Coord_orth;
using clipper::Coord_frac;
#include "clipper/core/cell.h"
#include "clipper/core/spacegroup.h"

void generateSymmpos(clipper::Cell * cell, clipper::Spacegroup * spgr,
        vector<float> & point, vector<vector<float> > & images,
        vector<float> & llb, vector<float> & rtf, float specialCutoff, float farCutoff) {
    Coord_orth target(point[0], point[1], point[2]);
    Coord_frac tcf = target.coord_frac( *cell );

    if(spgr->num_symops() < 1) { cout << "check spacegroup again! number of operators is " << spgr->num_symops() << endl; exit(1); }
    vector<float> im1(3); Coord_orth trn;

    for ( int isym = 1; isym < spgr->num_symops(); isym++ ) {
        //cout << spgr->symop(isym).format() << endl;
        //cout << "first ("<< scf[0] <<" "<< scf[1] <<" "<< scf[2] <<") ("<< image[0] <<" "<< image[1] <<" "<< image[2] <<")"<< endl;
        Coord_frac scf = spgr->symop(isym) * tcf;
        //for(int i=0; i < 1; i++) for(int j=0; j < 1; j++) for(int k=0; k < 1; k++) {
        for(int i=-2; i < 3; i++) for(int j=-2; j < 3; j++) for(int k=-2; k < 3; k++) {
            trn = (scf + Coord_frac(i,j,k)).coord_orth(*cell);
            im1[0] = trn[0]; im1[1] = trn[1]; im1[2] = trn[2];
            //if( fabs(im1[0]-target[0]) < specialCutoff && fabs(im1[1]-target[1]) < specialCutoff && fabs(im1[2]-target[2]) < specialCutoff ) continue; // special position
            if( fabs(im1[0] - target[0]) > farCutoff || fabs(im1[1] - target[1]) > farCutoff || fabs(im1[2] - target[2]) > farCutoff ) continue; // too far away
            if( rtf[0] < im1[0] || llb[0] > im1[0] ) continue; // not in box
            if( rtf[1] < im1[1] || llb[1] > im1[1] ) continue; // not in box
            if( rtf[2] < im1[2] || llb[2] > im1[2] ) continue; // not in box
            images.push_back(im1);
        }
    }
    //cout << "IMAGE " << point[0] <<" "<< point[1] <<" "<< point[2] << endl;
    //for(int im=0; im < images.size(); im++)
    //    cout << "IMAGE " << images[im][0] <<" "<< images[im][1] <<" "<< images[im][2] << endl;
}
#endif//CLIPPERGEOM

void findBoundingBox(vector<vector<float> > & pts, vector<float> & radii, vector<float> & llb, vector<float> & rtf, float margin) {
    llb.clear(); rtf.clear();
    for(int i=0; i < 3; i++) { llb.push_back(1e10); rtf.push_back(-1e10); }
    float maxrad = -100;
    for(int i=0; i < pts.size(); i++) {
        for(int k=0; k < 3; k++) {
            if(rtf[k] < pts[i][k]) rtf[k] = pts[i][k];
            if(llb[k] > pts[i][k]) llb[k] = pts[i][k];
        }
        if(maxrad < radii[i]) maxrad = radii[i];
    }
    for(int i=0; i < 3; i++) { llb[i] -= (margin+maxrad); rtf[i] += (margin+maxrad); }
}


void Grid::describe() {
    cout << "NUMPTS in GRID " << pi2vi.size() << endl;
}


Grid::Grid(int resol, vector<vector<float> >* pts, vector<float>* rads, GridHelper *GH) {
    resolution = resol;
    points = pts;
    radii = rads;
    maxrad = 1;
    gh = GH;
    sg = NULL; cell = NULL;
#ifdef CLIPPERGEOM
    llb.push_back(-1e10); llb.push_back(-1e10); llb.push_back(-1e10);
    rtf.push_back(1e10); rtf.push_back(1e10); rtf.push_back(1e10);
    findBoundingBox(*pts, *rads, llb, rtf, 20.); // approximately figure out the box containing points, bloated by 10A + max rad
    float a, b, c, A, B, C; string spacegroup;
    spacegroup = gh->cellsym(a, b, c, A, B, C);
    if(spacegroup == "") return;
    cell = reinterpret_cast<char*> (new clipper::Cell( clipper::Cell_descr(a,b,c,A,B,C) ));
    sg = reinterpret_cast<char*> (new clipper::Spacegroup( clipper::Spgr_descr(spacegroup.c_str()) ));
    specialCutoff = 0.2;
    farCutoff = ((clipper::Cell*)cell)->a();
    if(farCutoff < ((clipper::Cell*)cell)->b()) farCutoff = ((clipper::Cell*)cell)->b();
    if(farCutoff < ((clipper::Cell*)cell)->c()) farCutoff = ((clipper::Cell*)cell)->c();
    farCutoff = 500.;
    //cout << cell << endl;
    //cout << (reinterpret_cast<clipper::Spacegroup*>(sg))->num_symops() << endl;
#endif//CLIPPERGEOM
}

Grid::~Grid() {
#ifdef CLIPPERGEOM
    //cout << "here " << int(sg) <<" "<< int(cell) << endl;
    //if(sg != NULL) delete ((clipper::Spacegroup*)sg);
    //if(cell != NULL) delete ((clipper::Cell*)cell);
    sg = NULL; cell = NULL;
    //cout << "now here" << endl;
#endif//CLIPPERGEOM
}

void Grid::removes(vector<int> & pis) {
    for(int i=0; i < pis.size(); i++) remove(pis[i]);
}

void Grid::remove(int pi) {
    VoxelIndex& vi = pi2vi[pi];
    //cout << "removing " << pi << " from " << vi.i <<" "<< vi.j <<" "<< vi.k << endl;
    vi2pi[vi].erase(pi); pi2vi.erase(pi);
#ifdef CLIPPERGEOM
    if(friends.find(pi) != friends.end()) {
        for(int fi=0; fi < friends[pi].size(); fi++) {
            vi = friends[pi][fi];
            imageMap[vi].erase(pi);
        }
        friends.erase(pi);
    }
#endif//CLIPPERGEOM
}

bool Grid::voxelClashcheck(int pi, VoxelIndex& vi) {
    //cout << "VoxelClashCheck " << pi <<" ("<< vi.i <<" "<< vi.j <<" "<< vi.k <<") (" <<
    //    (*points)[pi][0]<<" "<< (*points)[pi][1]<<" "<< (*points)[pi][2]<<")"<< endl;
    if(vi2pi.find(vi) != vi2pi.end()) {
        set<int>& ptis = vi2pi[vi];
        for(set<int>::iterator it = ptis.begin(); it != ptis.end(); it++)
            if(true == gh->clash(*it, pi, (*radii)[*it], (*radii)[pi], (*points)[*it], (*points)[pi], 0))
                return true;
    }
#ifdef CLIPPERGEOM
    if(imageMap.find(vi) != imageMap.end()) {
        for(map<int, vector<float> >::iterator it = imageMap[vi].begin(); it != imageMap[vi].end(); ++it)
            if(true == gh->clash(it->first, pi, (*radii)[it->first], (*radii)[pi], it->second, (*points)[pi], 1))
                return true;
    }
#endif//CLIPPERGEOM
    return false;
}

void Grid::justAdd(VoxelIndex& vi, int pi) {
    if( vi2pi.find(vi) == vi2pi.end() ) vi2pi[vi] = set<int>();
    vi2pi[vi].insert(pi); pi2vi[pi] = vi;
    float radius = (*radii)[pi];
    if(maxrad < radius) maxrad = radius;
}

void Grid::printVox(VoxelIndex& vi) {
    cout << "VoxelInd " << vi.i<<" "<<vi.j<<" "<<vi.k<<" : ";
    if(vi2pi.find(vi) == vi2pi.end()) {
        cout << "empty" << endl; return;
    }
    set<int>& ptis = vi2pi[vi];
    for(set<int>::iterator it = ptis.begin(); it != ptis.end(); it++) cout << *it << " ";
    cout << endl;
}

// returns false only when a point clashes with one of its images
bool Grid::justAdd(vector<int>& ptinds) {
    VoxelIndex vi; bool retval = true;
    for(unsigned int pi=0; pi < ptinds.size(); pi++) {
        int pti = ptinds[pi];
        vector<float>& pt = (*points)[ pti ];
        xyz2ijk(pt[0],pt[1],pt[2], vi.i,vi.j,vi.k, resolution);
        //if(verbose(7)) cout << "justAdd-ing point " << pti <<" : "<< pt[0]<<" "<<pt[1]<<" "<<pt[2]<<" : "<<vi.i<<" "<<vi.j<<" "<<vi.k<< endl;
        justAdd(vi, pti);
#ifdef CLIPPERGEOM
        if(sg!=NULL) {
            vector< vector<float> > images;
            generateSymmpos((clipper::Cell*)cell, (clipper::Spacegroup*)sg, pt, images, llb, rtf, specialCutoff, farCutoff);
            if(friends.find(pti) == friends.end()) friends[pti] = vector<VoxelIndex>();
            for(int im=0; im < images.size(); im++) {
                xyz2ijk(images[im][0],images[im][1],images[im][2], vi.i,vi.j,vi.k, resolution);
                if(imageMap.find(vi) == imageMap.end()) imageMap[vi] = map<int, vector<float> >();
                imageMap[vi][pti] = images[im];
                friends[pti].push_back(vi);
//cout << "SEEE " << (*points)[pti][0] <<" "<< (*points)[pti][1] <<" "<< (*points)[pti][2] << endl;
//cout << "SEEE " << images[im][0] <<" "<< images[im][1] <<" "<< images[im][2] << endl;
                if(true == gh->clash(pti, pti, (*radii)[pti], (*radii)[pti], images[im], (*points)[pti], 1)) retval = false;
            }
        }
#endif//CLIPPERGEOM
    }
    return retval;
}

// envelope is the the union of spheres centered on points in the grid with given radius
// find if pt lies in the envelope
bool Grid::withinEnvelope(vector<float> & pt, float rad) {
    VoxelIndex vi;
    float dx, dy, dz, radsq = rad*rad;
    xyz2ijk(pt[0],pt[1],pt[2], vi.i,vi.j,vi.k, resolution);
    int margin = (int) ( 1 + ceil(rad / resolution) );
    for(int i = vi.i-margin; i <=  vi.i + margin; i++)
        for(int j = vi.j-margin; j <=  vi.j + margin; j++)
            for(int k = vi.k-margin; k <=  vi.k + margin; k++) {
                VoxelIndex voxind(i,j,k);
                if(vi2pi.find(voxind) == vi2pi.end()) continue;
                set<int>::iterator it = vi2pi[voxind].begin(), end = vi2pi[voxind].end();
                for(; it != end; ++it) {
                    dx = fabs((*points)[*it][0] - pt[0]);
                    if(dx > rad) continue;
                    dy = fabs((*points)[*it][1] - pt[1]);
                    if(dy > rad) continue;
                    dz = fabs((*points)[*it][2] - pt[2]);
                    if(dz > rad) continue;
                    if(dx*dx + dy*dy + dz*dz < radsq) return true;
                }
            }
    return false;
}

bool Grid::add(vector<int>& ptinds) {
    //for(int i=0; i < ptinds.size(); i++) cout << "ADD " << ptinds[i] << endl;
    VoxelIndex vi;
    map<int,VoxelIndex> pi2viMap;
    for(unsigned int pi=0; pi < ptinds.size(); pi++) {
        int pti = ptinds[pi];
        vector<float>& pt = (*points)[pti];
        xyz2ijk(pt[0],pt[1],pt[2], vi.i,vi.j,vi.k, resolution);
        pi2viMap[pti] = vi;
        int margin = (int) ( ceil( (maxrad + (*radii)[pti]) / resolution ) );
        //cout << "margin " << margin <<" "<< maxrad <<" "<< (*radii)[pti] << endl;
        bool clashing = false;
        for(int i = vi.i-margin; !clashing && i <=  vi.i + margin; i++)
            for(int j = vi.j-margin; !clashing && j <=  vi.j + margin; j++)
                for(int k = vi.k-margin; !clashing && k <=  vi.k + margin; k++) {
                    VoxelIndex voxind(i,j,k);
                    //cout << "VoxelClashCheck " << pti <<" ("<< voxind.i <<" "<< voxind.j <<" "<< voxind.k <<") (" <<
                    //    (*points)[pti][0]<<" "<< (*points)[pti][1]<<" "<< (*points)[pti][2]<<")";
                    //if(vi2pi.find(voxind) == vi2pi.end()) cout << " absent"; cout << endl;
                    //if(vi2pi.find(voxind) == vi2pi.end()) continue;
                    if( voxelClashcheck(pti, voxind) ) clashing = true;
                }
        if(clashing) return false;
    }
    if(false == justAdd(ptinds)) { removes(ptinds); return false; }
    return true;
}

void Grid::xyz2ijk(float x, float y, float z, int& vi, int& vj, int& vk, int resol) {
    vi = (int) round(x / resol);
    vj = (int) round(y / resol);
    vk = (int) round(z / resol);
}

void Grid::printIJK(vector<float> &pt) {
    int i, j, k;
    xyz2ijk(pt[0], pt[1], pt[2], i, j, k, resolution);
    cout << "(" << pt[0] <<","<< pt[1] << "," << pt[2] << ") -> ";
    cout << "(" << i <<","<< j << "," << k << ")" << endl;
}

// for population-search strategy
// copy arg's contents into oneself, ie pi2vi, vi2pi, maxrad
// other members need not be copied bcz points and grids are one-to-one and points are copied
// resolution, overlap-reduction, radii are same for all
void Grid::popCopyGridPoints(Grid & g) {
    pi2vi = g.pi2vi;
    vi2pi = g.vi2pi;
    maxrad = g.maxrad;
    vector<vector<float> > & pts = *points;
    vector<vector<float> > & gpts = *(g.points);
    pts = gpts;

#ifdef CLIPPERGEOM
    sg = g.sg; cell = g.cell;
    specialCutoff = g.specialCutoff; farCutoff = g.farCutoff;
    friends = g.friends; imageMap = g.imageMap;
#endif//CLIPPERGEOM
}

//{{{ OverlapMap functions
MapIntIntFloat::MapIntIntFloat() {}

bool MapIntIntFloat::contains(int a, int b) {
    if(miif.find(a) == miif.end()) return false;
    if(miif[a].find(b) == miif[a].end()) return false;
    return true;
}

float MapIntIntFloat::get(int a, int b) {
    return miif[a][b];
}

void MapIntIntFloat::put(int a, int b, float f) {
    if(miif.find(a) == miif.end()) miif[a] = map<int,float>();
    if(miif.find(b) == miif.end()) miif[b] = map<int,float>();
    miif[a][b] = f;
    miif[b][a] = f;
}
//}}}


