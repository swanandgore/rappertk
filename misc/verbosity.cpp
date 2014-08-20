#include "verbosity.h"
#include <iostream>
using namespace std;

Verbosity::Verbosity() { verbo = 0; }

void Verbosity::setVerbo(int level) {
    verbo = level;
    //cout << this << " VERBOSITY set at " << verbo << endl;
}

bool Verbosity::isLevelLEq(int level) {
    //cout << this << " VERBOSITY " << verbo << endl;
    return level <= verbo;
}

Verbosity & Verbosity::instance() {
    static Verbosity verbinst;
    return verbinst;
}

void setVerbosity(int v) {
    Verbosity::instance().setVerbo(v);
}

bool verbose(int level) {
    //static Verbosity & myinst = Verbosity::instance();
    //return myinst.isLevelGEq(level);
    return Verbosity::instance().isLevelLEq(level);
}
