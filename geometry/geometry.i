%module geometry

%{
#include "Grid.h"
#include "GridHelper.h"
#include "CAtraceGH.h"
#include "functions.h"
%}

%include "std_vector.i"
namespace std {
    %template(VecFloat) vector<float>;
    %template(VecLong) vector<long>;
    %template(VecInt) vector<int>;
    %template(VecVecInt) vector<vector<int> >;
}

%include "Grid.h"
%include "GridHelper.h"
%include "CAtraceGH.h"
%include "functions.h"
