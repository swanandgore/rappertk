%module samplers

%{
#include "PhipsiSampler.h"
#include "InformedPhipsiSampler.h"
#include "BBdepChiSampler.h"
#include "PRLsampler.h"
#include "SCLsampler.h"
#include "RapSampler.h"
#include "OmegaSampler.h"
#include "NAsampler.h"
#include "NAsuiteSampler.h"
#include "RATdata.h"
#include "FragSampler.h"
%}

%include "std_vector.i"
namespace std {
}

%include "PhipsiSampler.h"
%include "InformedPhipsiSampler.h"
%include "OmegaSampler.h"
%include "BBdepChiSampler.h"
%include "PRLsampler.h"
%include "NAsampler.h"
%include "NAsuiteSampler.h"
%include "RapSampler.h"
%include "RATdata.h"
%include "FragSampler.h"
%include "SCLsampler.h"
