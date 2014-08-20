#include <stdlib.h>

#include "misc/RanGen.h"
#include "OmegaSampler.h"


#include <map>
using std::map;
#include <string>
using std::string;

#include<iostream>; 
using namespace std ; 

OmegaSampler::OmegaSampler(float trprop) : transprop(trprop) {}

int OmegaSampler::sample(int& omega) {
  //  cout << "trans" << transprop << " " <<  ran01() <<  endl ; 
  if(ran01() > transprop) { omega = 0; return 0; }
  else { omega = -180; return 1; }
}

int OmegaSampler::maxUniqueSamples() { return 2; }
