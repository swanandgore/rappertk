#include <math.h>
#include "FragSampler.h"
#include "misc/RanGen.h"

FragSampler::FragSampler(int numfrags, vector<vector<float> > & pts) {
    if(pts.size() % numfrags != 0) { cout << "frag input size " << pts.size() << " is not divisible by numsamples " << numfrags << endl; exit(1); }
    for(int i=0; i < numfrags; i++) {
        int start = pts.size() / numfrags * i;
        int end = pts.size() / numfrags * (i+1);
        samples.push_back(vector<vector<float> >());
        for(int k=start; k < end; k++) samples[i].push_back(vector<float>(pts[k]));
    }
}

int FragSampler::sample(vector<vector<float> > & pts) {
    int myrand = (int) floor ( pts.size() * ran01() );
    pts = samples[myrand];
}

int FragSampler::maxUniqueSamples() { return samples.size(); }
