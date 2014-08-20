#include <vector>
#include "functions.h"
using namespace std;

int main() {
    vector<float> X(3), Y(3), Z(3);
    for(int i=0; i < 1000; i++) {
        randomNormalVector(X);
        cout << X[0] <<" "<< X[1] <<" "<< X[2] << endl;
        continue;
        X[0] = 1; X[1] = 0; X[2] = 0;
        findYZgivenX(X, Y, Z);
        cout << "X " << X[0] <<" "<< X[1] <<" "<< X[2] << endl;
        cout << "Y " << Y[0] <<" "<< Y[1] <<" "<< Y[2] << endl;
        cout << "Z " << Z[0] <<" "<< Z[1] <<" "<< Z[2] << endl;
    }
    return 0;

//    vector<float> c1(3), c2(3), q(3), p1(3), p2(3);
//    float r1=1, r2=1, angle=135;
//
//    c1[0] = 0; c1[1] = 0; c1[2] = 0;
//    c2[0] = 1; c2[1] = 1; c2[2] = 0;
//    q[0] = 0; q[1] = -1; q[2] = 0;
//
//    cout << findSphereSphereAngleIntx(r1, c1, r2, c2, q, angle, p1, p2) << endl;
}
