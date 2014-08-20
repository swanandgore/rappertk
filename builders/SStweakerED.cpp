#include "SStweakerED.h"

SStweakerED::SStweakerED(vector<int>& opInds, const char* desc) : op(opInds), description(desc) {}

// translate and rotate randomly within the given limits
void SStweakerED::build(VVF & pts) {}

vector<int>& SStweakerED::getOP() { return op; }

// count how many points are in > 1 sigma density
float SStweakerED::energy() {
    
}

const char* SStweakerED::name() { return description.c_str(); }
