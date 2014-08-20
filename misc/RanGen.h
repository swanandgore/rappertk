#ifndef __RANGEN_H__
#define __RANGEN_H__

class RanGen {
public :
    static RanGen & instance();
    double ran_uni_0_1();
    void seedme(int r);
private :
    RanGen();
    double ran2(long *idum);
    static long seed;
};

#include "stdlib.h"
#define ran01() (( RanGen::instance().ran_uni_0_1() ))
//#define ran01() (( (rand()+0.)/(RAND_MAX+0.) ))

#endif // __RANGEN_H__
