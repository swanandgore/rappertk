#include <assert.h>

#include "Builder.h"
#include "Constants.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
#include <iostream>
using namespace std;

Builder::Builder(vector<int>& ipInds, vector<int>& opInds, Constants *con, const char* desc)
    : ip(ipInds), op(opInds), consts(con), description(desc), session(false), maxSessionSize(10000000) {}

Builder::Builder() {}

Builder::~Builder() {}

void Builder::describe() {
    cout << "(";
    for(int i=0; i < ip.size(); i++) cout << ip[i] << " ";
    cout << ") (";
    for(int i=0; i < op.size(); i++) cout << op[i] << " ";
    cout << ") ";
    cout << description;
}

const char* Builder::name() {
    return description.c_str();
}

bool Builder::build(VVF & pts) { exit(0); return false; }

bool Builder::buildSample(VVF & pts, int sampleIndex) { cout << "sample-index building not implemeted" << endl; exit(0); }

vector<int>& Builder::getIP() { return ip; }

vector<int>& Builder::getOP() { return op; }

void Builder::startSession(int ssize) {
    if(session == true) { cout<< "Cant start another session when already in session, fatal error" << endl; exit(0); }
    maxSessionSize = ssize; sessionSize = 0;
    deleteSessionInfo();
    session = true;
}

void Builder::endSession() {
    deleteSessionInfo();
    session = false;
}

void Builder::deleteSessionInfo() {}
