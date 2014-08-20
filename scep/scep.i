%module scep

%{
#include "SCenergyProvider.h"
#include "biconn.h"
%}

%include "std_vector.i"
namespace std {
    %template(VecInt) vector<int>;
    %template(VecVecInt) vector< vector<int> >;
    %template(VecFloat) vector<float>;
    %template(VecVecFloat) vector<vector <float> >;
    %template(VecBuilder) vector<Builder*>;
}

%include "SCenergyProvider.h"
%include "biconn.h"
