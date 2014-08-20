#ifndef __ENVELOPERESTRAINT_H__
#define __ENVELOPERESTRAINT_H__

#include "Restraint.h"

class Grid;

class EnvelopeRestraint : public Restraint {
public :
    EnvelopeRestraint(vector<int> & pis, const char *desc, Grid *grid, float radius);
    bool satisfied(VVF & pts);
    void describe();
private :
    Grid *grid;
    float envRad;
};

#endif//__ENVELOPERESTRAINT_H__
