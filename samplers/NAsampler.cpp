#include <assert.h>
#include "NAsampler.h"
#include "misc/RanGen.h"

#include <stdlib.h>

NAsampler::NAsampler() {}

// delta and chi are read from constants according to base and pucker, other dihedrals need sampling
void NAsampler::sample(char base, string & pucker, float& alpha, float& beta, float& gamma, float& epsilon, float& zeta) {
    alpha = ran01() * 360;
    beta = ran01() * 360;
    gamma = ran01() * 360;
    epsilon = ran01() * 360;
    zeta = ran01() * 360;
}
