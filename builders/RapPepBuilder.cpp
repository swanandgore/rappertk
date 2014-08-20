#include <math.h>
#include "RapPepBuilder.h"
#include "Constants.h"
#include "samplers/RapSampler.h"
#include "samplers/OmegaSampler.h"
#include "geometry/functions.h"
#include "misc/verbosity.h"

RapPepBuilder::RapPepBuilder(vector<int>& ipInds, vector<int>& opInds,
        Constants* con, const char *resn, const char *desc, RapSampler *rs, OmegaSampler *o, bool fwdbkwd)
        : BasicPepBuilder(ipInds, opInds, con, resn, desc, fwdbkwd)
{
    raps = rs; os = o;
}

void RapPepBuilder::sample(VVF & pts, float & phi, float & psi, float & omega) {
    raps->sample(phi, psi);
    int o;
    os->sample(o); omega = o;
    //if(verbose(7)) cout << "sampled phi-psi-omega " << phi <<" "<< psi <<" "<< omega << endl;
}
