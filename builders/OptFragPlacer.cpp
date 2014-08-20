#include "OptFragPlacer.h"
#include "samplers/FragSampler.h"
#include "geometry/functions.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
OptFragPlacer::OptFragPlacer(vector<int>& ipInds, vector<int>& opInds, FragSampler *fsamp, vector<int>& fis, vector<int>& fos, const char* desc)
        : Builder(ipInds, opInds, NULL, desc), fs(fsamp), fragips(fis), fragops(fos) {
    if(fragips.size() < 3 || fragips.size() != ip.size() || fragops.size() != op.size()) {
        cout <<"Sizes of corresponding pointsets not same or too small "<< ip.size() <<" "<< fragips.size() <<" "<< op.size() <<" "<< fragops.size() << endl;
        exit(1);
    }
}

// this builder will sample a fragment, superpose the sample's fraginds onto ipInds in pts to get the optimal transformation, and apply it to whole fragment
bool OptFragPlacer::build(VVF & pts) {
    vector<vector<float > > sample;
    fs->sample(sample);
    if(sample.size() != fragips.size() + fragops.size()) { cout << "fragsampler sample has unexpected size" << endl; exit(1); }
    if(ip.size() != fragips.size() or op.size() != fragops.size()) { cout << "fragsampler sample inpu sizes incorrect" << endl; exit(1); }

    vector<vector<float> > from, onto;
    for(int i=0; i < fragips.size(); i++) {
        from.push_back(sample[fragips[i]]);
        onto.push_back(pts[ip[i]]);
        //cout << sample[fragips[i]][0]<<" "<< sample[fragips[i]][1]<<" "<< sample[fragips[i]][2]<<" "<<
        //  pts[ip[i]][0]<<" "<< pts[ip[i]][1]<<" "<< pts[ip[i]][2]<<" "<<endl;
    }

    //cout << calcDist(from[0], from[1]) <<" "<< calcDist(onto[0], onto[1]) << endl;
    //cout << calcDist(from[1], from[2]) <<" "<< calcDist(onto[1], onto[2]) << endl;
    //cout << calcAngle(from[0], from[1], from[2]) <<" "<< calcAngle(onto[0], onto[1], onto[2]) << endl;
    //exit(0);

    vector<float> t1(3,0), t2(3,0); vector<vector<float> > rot(3, vector<float>(3,0));
    float rmsd = findSuperpositionTransform(from, onto, t1, rot, t2);
    //cout << "superpose " << rmsd << endl;

    for(int i=0; i < fragops.size(); i++) {
        //cout << i <<" "<< op[i] <<" "<< fragops[i] <<" "<< fragops.size() << endl;
        pts[op[i]] = sample[fragops[i]];
        TRTtransform1(pts[op[i]], t1, rot, t2);
    }
    return true;
}

OptFragPlacer OptFragPlacer::makeCopy() {
    return OptFragPlacer(*this);
}
