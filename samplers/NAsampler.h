#ifndef __NASAMPLER_H__
#define __NASAMPLER_H__

#include <string>
using namespace std;

class NAsampler {
public :
    NAsampler();
    void sample(char base, string & pucker, float& alpha, float& beta, float& gamma, float& epsilon, float& zeta);
private :
};

#endif // __NASAMPLER_H__
