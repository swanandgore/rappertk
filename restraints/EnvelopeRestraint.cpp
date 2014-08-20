#include "EnvelopeRestraint.h"
#include "geometry/Grid.h"

EnvelopeRestraint::EnvelopeRestraint(vector<int> & pis, const char *desc, Grid *gr, float radius)
    : Restraint(pis,desc), grid(gr), envRad(radius)
{}

bool EnvelopeRestraint::satisfied(VVF & pts) {
    for(int i=0; i < ptInds.size(); i++)
        if( ! grid->withinEnvelope(pts[ptInds[i]], envRad) ) return false;
    return true;
}

void EnvelopeRestraint::describe() { cout << "EnvelopeRestraint with radius " << envRad; }
