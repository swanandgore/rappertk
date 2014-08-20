%module builders

%{
#include "Builder.h"
#include "DihedBuilder.h"
#include "PeptideBuilder.h"
#include "PeptideBridgeBuilder.h"
#include "RapPepBuilder.h"
#include "BasicPepBuilder.h"
#include "InformedPeptideBuilder.h"
#include "CNCaBuilder.h"
#include "CBbuilder.h"
#include "ChiBuilder.h"
#include "Constants.h"
#include "NAbuilder.h"
#include "NAsuiteBuilder.h"
#include "NanchorBuilder.h"
#include "BuilderGroup.h"
#include "ProtResBuilder.h"
#include "TransRotator.h"
#include "Rotator.h"
#include "OptFragPlacer.h"
#include "SCLbuilder.h"
%}

%include "std_vector.i"
namespace std {
    %template(VecInt) vector<int>;
    %template(VecVecInt) vector< vector<int> >;
    %template(VecFloat) vector<float>;
    %template(VecVecFloat) vector<vector <float> >;
    %template(VecBuilder) vector<Builder*>;
}

%include "Builder.h"
%include "DihedBuilder.h"
%include "BasicPepBuilder.h"
%include "PeptideBuilder.h"
%include "PeptideBridgeBuilder.h"
%include "CBbuilder.h"
%include "ChiBuilder.h"
%include "Constants.h"
%include "NAbuilder.h"
%include "NAsuiteBuilder.h"
%include "CNCaBuilder.h"
%include "InformedPeptideBuilder.h"
%include "NanchorBuilder.h"
%include "BuilderGroup.h"
%include "RapPepBuilder.h"
%include "ProtResBuilder.h"
%include "TransRotator.h"
%include "Rotator.h"
%include "OptFragPlacer.h"
%include "SCLbuilder.h"
