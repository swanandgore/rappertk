
#include "Restraint.h"
#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>
#include <iostream>
using namespace std;

Restraint::Restraint(vector<int>& pis, const char *desc)
    : ptInds(pis), description(desc) {}

bool Restraint::satisfied(VVF & pts) { return true; }

Restraint::~Restraint() {}

vector<int>& Restraint::getInds() { return ptInds; }

void Restraint::describe() {
    cout << "this restraint description is incomplete." << endl;
    exit(1);
}

const char * Restraint::name() { return description.c_str(); }
