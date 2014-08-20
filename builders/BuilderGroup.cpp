#include "misc/verbosity.h"
#include "BuilderGroup.h"
#include "Constants.h"

#include <iostream>
#include <set>
using namespace std;

BuilderGroup::BuilderGroup(vector<Builder*>& bldrs, const char* desc) {
    builders = bldrs;
    description = desc;

    //cout <<"HERE" << endl;
    set<int> tip, top;
    for(int bi=0; bi < builders.size(); bi++) {
        vector<int> & biIP = builders[bi]->getIP();
        vector<int> & biOP = builders[bi]->getOP();
        //for(int i=0; i < biIP.size(); i++) cout << biIP[i] <<" "; cout << " : ";
        //for(int i=0; i < biOP.size(); i++) cout << biOP[i] <<" "; cout << endl;
        for(int i=0; i < biIP.size(); i++) tip.insert(biIP[i]);
        for(int i=0; i < biOP.size(); i++) top.insert(biOP[i]);
    }
    ip.clear(); op.clear();
    for(set<int>::iterator it = tip.begin(); it != tip.end(); ++it) // put tip[i] in ip if it is not in top
        if(top.find(*it) == top.end()) ip.push_back(*it);
    for(set<int>::iterator it = top.begin(); it != top.end(); ++it) // put top[i] in op unconditionally
        op.push_back(*it);

    //for(int i=0; i < ip.size(); i++) cout << ip[i] << " "; cout << " : ";
    //for(int i=0; i < op.size(); i++) cout << op[i] << " "; cout << endl;
    //cout <<"HERE------------------------------" << endl;
}

bool BuilderGroup::build(vector<vector<float> > & pts) {
    for(int bi=0; bi < builders.size(); bi++) {
        bool bstate = builders[bi]->build(pts);
        //if(verbose(7)) cout << "trying subbuilder " << builders[bi]->name() <<" : "<< bstate << endl;
        if(! bstate ) return false;
    }
    return true;
}

// assuming that it is safe to use constituent builders as they are in stateless mode
BuilderGroup BuilderGroup::makeCopy() {
    return BuilderGroup(*this);
}
