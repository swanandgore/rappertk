#include <math.h>
#include "PeptideBuilder.h"
#include "Constants.h"
#include "samplers/PhipsiSampler.h"
#include "samplers/OmegaSampler.h"
#include "geometry/functions.h"
#include "misc/verbosity.h"

PeptideBuilder::PeptideBuilder(vector<int>& ipInds, vector<int>& opInds,
        Constants* con, const char *resn, const char *desc, PhipsiSampler* p, OmegaSampler *o, bool fwdbkwd)
        : BasicPepBuilder(ipInds, opInds, con, resn, desc, fwdbkwd)
{
    pps = p; os = o;
}

void PeptideBuilder::sample(VVF & pts, float & phi, float & psi, float & omega) {
    int p, s, o;
    pps->sample(p,s); os->sample(o);
    phi = p; psi = s; omega = o;
    //    if(verbose(7)) cout << "sampled phi-psi-omega " << phi <<" "<< psi <<" "<< omega << endl;
}

PeptideBuilder PeptideBuilder::makeCopy() {
    return PeptideBuilder(*this);
}