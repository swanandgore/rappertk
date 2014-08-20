#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <string>
#include <map>
using namespace std;

class Constants {
public:
    Constants();
    void set(const char*, float);
    void set(string&, float);
    float get(const char *);
    float get(string&);
private:
    map<string,float> floatConsts;
};

#endif //__CONSTANTS_H__
