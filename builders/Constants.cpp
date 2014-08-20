#include "Constants.h"

#include <assert.h>
#include <iostream>
using namespace std;

Constants::Constants() {}

void Constants::set(string& k, float v) {
    assert(floatConsts.find(k) == floatConsts.end());
    floatConsts[k] = v;
}

void Constants::set(const char* k, float v) {
    string s(k);
    this->set(s, v);
}

float Constants::get(string& k) {
    map<string,float>::iterator it = floatConsts.find(k);
    assert(it != floatConsts.end());
    return it->second; //floatConsts[k];
}
float Constants::get(const char* k) {
    string str(k);
    return get(str);
}
