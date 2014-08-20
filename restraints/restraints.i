%module restraints

%{
#include "Restraint.h"
#include "AngleRestraint.h"
#include "DihedRestraint.h"
#include "DistanceRestraint.h"
#include "SphPosRestr.h"
#include "CentroidPosRestraint.h"
#include "RATrestraint.h"
#include "EDrestraint.h"
#include "EnvelopeRestraint.h"
%}

%include "std_vector.i"
namespace std {
}

%include "Restraint.h"
%include "DistanceRestraint.h"
%include "AngleRestraint.h"
%include "DihedRestraint.h"
%include "SphPosRestr.h"
%include "CentroidPosRestraint.h"
%include "RATrestraint.h"
%include "EDrestraint.h"
%include "EnvelopeRestraint.h"
