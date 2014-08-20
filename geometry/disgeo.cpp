#include <assert.h>
#include <vector>
#include <set>
#include <map>
#include <iostream>
using namespace std;

#define lower(mat,i,j) (i>j?mat[i][j]:mat[j][i])
#define upper(mat,i,j) (i<j?mat[i][j]:mat[j][i])

typedef pair<int,int> PII;
typedef set<PII> SPII;

template<class floatmat>
void checkModifyBounds(int i, int j, int k, SPII & changedBounds, floatmat & ul) {
    float uij = upper(ul,i,j), ujk = upper(ul,j,k), uki = upper(ul,k,i);
    //cout << "IJK " << i <<" "<< j <<" "<< k <<" "<< uij <<" "<< ujk <<" " << uki << endl;
    if     (uij > ujk + uki) {upper(ul,i,j) = ujk + uki; changedBounds.insert(PII(i,j));}
    else if(ujk > uki + uij) {upper(ul,j,k) = uki + uij; changedBounds.insert(PII(j,k));}
    else if(uki > uij + ujk) {upper(ul,k,i) = uki + uij; changedBounds.insert(PII(i,k));}
    //else cout << "smoothing not done" << endl;
}

//assume that lower triangular (i>j) is lower bounds, upper is upper bounds (i<j) and dont use diagonal
template<class floatmat>
void boundSmoothing(floatmat & ul, int numAtoms) {
    for(int i=0; i < numAtoms; i++) {
	for(int j=0; j < numAtoms; j++) {
	    cout << ul[i][j] << " ";
	    assert(upper(ul,i,j) >= lower(ul,i,j));
	}
	cout << endl;
    }
    
    // find initial changed limits for edges i<j
    SPII changedBounds;
    for(int i=0; i < numAtoms; i++) {
	for(int j=i+1; j < numAtoms; j++)
	    for(int k=j+1; k < numAtoms; k++) { // all triangles
		checkModifyBounds(i,j,k, changedBounds, ul);
	    }
    }
    // for each changed bound, check all triangles
    while(changedBounds.size() > 0) {
	//cout << "changedBounds--------------------" << endl;
	//for(SPII::iterator it = changedBounds.begin(); it != changedBounds.end(); it++)
	//    cout << it->first <<","<< it->second << endl;
	PII bound = *(changedBounds.begin());
	changedBounds.erase(bound);
	for(int k=0; k < numAtoms; k++) {
	    if(k == bound.first || k == bound.second) continue;
	    checkModifyBounds(bound.first, bound.second, k, changedBounds, ul);
	}
    }

    cout << "---------------------------------" << endl;
    for(int i=0; i < numAtoms; i++) {
	for(int j=0; j < numAtoms; j++) {
	    cout << ul[i][j] << " ";
	    assert(upper(ul,i,j) >= lower(ul,i,j));
	}
	cout << endl;
    }
    cout << "---------------------------------" << endl;
}

main() {
    int mlen = 20;
    vector<vector<float> > flmat(mlen);
    for(int i=0; i < mlen; i++) flmat[i].resize(mlen);
    for(int i=0; i < mlen; i++)
	for(int k=i+1; k < mlen; k++) {
	    lower(flmat, i, k) = 1;
	    upper(flmat, i, k) = 999;
	}
    for(int i=0; i < mlen-1; i++) upper(flmat, i, i+1) = 1;
    upper(flmat, 0, mlen-1) = 1;

    boundSmoothing(flmat, mlen);

    for(int i=0; i < mlen; i++)
	cout << "final " << i << " " << mlen-1 << " " << upper(flmat, i, mlen-1) << endl;
}
