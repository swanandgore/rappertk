#include <assert.h>
#include "functions.h"
#include "misc/RanGen.h"

// find the optimal superposition transform (translation, rotation, translation)
// return rmsd resulting from it
// both pointsets are unchanged
float findSuperpositionTransform( vector<vector<float> > & from, vector<vector<float> > & onto, vector<float> & t1, vector<vector<float> > & rot, vector<float> & t2 ) {
    if(from.size() != onto.size()) { cout << "superpose called with unequal pointsets" << endl; exit(1); }
    // find center of mass and translate the pointsets so that they are at origin
    vector<float> CMfrom(3), CMonto(3);
    for(int i=0; i < 3; i++) { CMfrom[i] = 0; CMonto[i] = 0.; }
    for(int pi=0; pi < from.size(); pi++)
        for(int i=0; i < 3; i++) {
            CMfrom[i] += from[pi][i];
            CMonto[i] += onto[pi][i];
        }
    for(int i=0; i < 3; i++) { CMfrom[i] /= from.size(); CMonto[i] /= from.size(); }
    // find distance matrix
    vector<float> q(4,0), m(3,0), p(3,0); // minus and plus
    vector<vector<float> > Q(4,q); // quaternion
    //for(int qi=0; qi < 4; qi++) cout << "Q " << Q[qi][0] <<" "<< Q[qi][1] <<" "<< Q[qi][2] <<" "<< Q[qi][3] << endl;
    float xm, ym, zm, xp, yp, zp;
    for(int pi=0; pi < from.size(); pi++) {
        xm = from[pi][0]-CMfrom[0] - onto[pi][0]+CMonto[0]; xp = onto[pi][0]-CMonto[0] + from[pi][0]-CMfrom[0];
        ym = from[pi][1]-CMfrom[1] - onto[pi][1]+CMonto[1]; yp = onto[pi][1]-CMonto[1] + from[pi][1]-CMfrom[1];
        zm = from[pi][2]-CMfrom[2] - onto[pi][2]+CMonto[2]; zp = onto[pi][2]-CMonto[2] + from[pi][2]-CMfrom[2];
        Q[0][0] += xm*xm + ym*ym + zm*zm;
        Q[1][1] += yp*yp + zp*zp + xm*xm;
        Q[2][2] += xp*xp + zp*zp + ym*ym;
        Q[3][3] += xp*xp + yp*yp + zm*zm;
        Q[0][1] += yp*zm - ym*zp;
        Q[0][2] += xm*zp - xp*zm;
        Q[0][3] += xp*ym - xm*yp;
        Q[1][2] += xm*ym - xp*yp;
        Q[1][3] += xm*zm - xp*zp;
        Q[2][3] += ym*zm - yp*zp;
    }
    Q[1][0]=Q[0][1];
    Q[2][0]=Q[0][2];
    Q[2][1]=Q[1][2];
    Q[3][0]=Q[0][3];
    Q[3][1]=Q[1][3];
    Q[3][2]=Q[2][3];

    //for(int qi=0; qi < 4; qi++) cout << "Q " << Q[qi][0] <<" "<< Q[qi][1] <<" "<< Q[qi][2] <<" "<< Q[qi][3] << endl;

    vector<float> eigvals(4,0);
    vector<vector<float> > eigvecs(4,eigvals);
    if(Jacobi(Q, eigvals, eigvecs)) { cout << "Jacobi did not converge" << endl; exit(1); }
    //for(int ei=0; ei < 4; ei++) cout <<"EIG "<< eigvals[ei] <<" : "<< eigvecs[0][ei] <<" "<< eigvecs[1][ei] <<" "<< eigvecs[2][ei] <<" "<< eigvecs[3][ei] << endl;

    float rmsd;
    if(eigvals[3] < 0) rmsd = 0;
    else rmsd = sqrt(eigvals[3]/from.size());

    // construct rotation matrix from eigenvector with lowest eigenvalue
    float q1 = eigvecs[0][3], q2 = eigvecs[1][3], q3 = eigvecs[2][3], q4 = eigvecs[3][3];
    float q1sqr=q1*q1, q2sqr=q2*q2, q3sqr=q3*q3, q4sqr=q4*q4;
    float q1q2=q1*q2, q1q3=q1*q3, q1q4=q1*q4;
    float q2q3=q2*q3, q2q4=q2*q4, q3q4=q3*q4;

    rot[0][0] = q1sqr + q2sqr - q3sqr -q4sqr;
    rot[0][1] = 2.0 * (q2q3 - q1q4);
    rot[0][2] = 2.0 * (q2q4 + q1q3);

    rot[1][0] = 2.0 * (q2q3 + q1q4);
    rot[1][1] = q1sqr + q3sqr - q2sqr - q4sqr;
    rot[1][2] = 2.0 * (q3q4 - q1q2);

    rot[2][0] = 2.0 * (q2q4 - q1q3);
    rot[2][1] = 2.0 * (q3q4 + q1q2);
    rot[2][2] = q1sqr + q4sqr -q2sqr - q3sqr;

    for(int i=0; i < 3; i++) t1[i] = -1 * CMfrom[i];
    for(int i=0; i < 3; i++) t2[i] = CMonto[i];

    return rmsd;
}

// translate, rotate, translate
void TRTtransform(vector<vector<float> > & pts, vector<float> & t1, vector<vector<float> > & rot, vector<float> & t2) {
    for(int pi=0; pi < pts.size(); pi++) TRTtransform1(pts[pi], t1, rot, t2);
}

void TRTtransform1(vector<float> & pt, vector<float> & t1, vector<vector<float> > & rot, vector<float> & t2) {
    for(int i=0; i < 3; i++) pt[i] += t1[i];
    float temp[3];
    for(int i=0; i < 3; i++) temp[i] = rot[i][0]*pt[0] + rot[i][1]*pt[1] + rot[i][2]*pt[2];
    for(int i=0; i < 3; i++) pt[i] = temp[i];
    for(int i=0; i < 3; i++) pt[i] += t2[i];
}


// find base as projection of b on c, and perp as height on base c
void findTriangleBaseHeight(float a, float b, float c, float &base, float &perp) {
    base = (b*b + c*c - a*a) / (2*c);
    perp = sqrt(b*b - base*base); 
}   

// generates a random normalized vector
// is uniformly distributed over sphere surface
// taken from sphere point picking in mathworld
void randomNormalVector(vector<float> & v) {
//    while(1) {
//        v[0] = -1. + 2.* ran01();
//        v[1] = -1. + 2.* ran01();
//        v[2] = 1 - v[0]*v[0] + v[1]*v[1];
//        if(v[2] < 0) continue;
//        v[2] = sqrt(v[2]);
//        return;
//    }
    float theta = 2. * M_PI * ran01();
    float phi = acos(2. * ran01() - 1);
    v[0] = sin(phi) * cos(theta);
    v[1] = sin(phi) * sin(theta);
    v[2] = cos(phi);
}

void sphereVolSample(vector<float>& cen, float rad, vector<float>& sample) {
    float dx, dy, dz;
    while(1) {
        dx = (ran01() * 2 - 1) * rad;
        dy = (ran01() * 2 - 1) * rad;
        dz = rad*rad - dx*dx - dy*dy;
        if(dz < 0) continue;
        dz = sqrt(dz);
        if(ran01() > 0.5) dz *= -1;
        sample[0] = cen[0] + dx;
        sample[1] = cen[1] + dy;
        sample[2] = cen[2] + dz;
        return;
    }
}


float dihedDiff(float d0, float d1) {
    d0 = withinPlusMinus180(d0);
    d1 = withinPlusMinus180(d1);
    float diff = d0 - d1;
    if(diff > 180) diff = 360 - diff;
    else if(diff < -180) diff = -360 - diff;
    return diff;
}

float withinPlusMinus180(float num) {
    while(num < -180) num += 360;
    while(num >= 180) num -= 360;
    return num;
}

void linear_combination(vector<float>& pt, float a, vector<float>& A, float b, vector<float>& B) {
    pt[0] = a*A[0]+b*B[0];
    pt[1] = a*A[1]+b*B[1];
    pt[2] = a*A[2]+b*B[2];
}

float dot_product(const vector<float>& a, const vector<float>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

void point_sub(vector<float>& p, const vector<float>& a, const vector<float>& b) {
    p[0] = a[0] - b[0]; p[1] = a[1] - b[1]; p[2] = a[2] - b[2];
}

float magnitude(vector<float>& p) { return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]); }

void normalize_point(vector<float>& n, vector<float>& p) {
    float mag = magnitude(p);
    n[0]=p[0]/mag; n[1]=p[1]/mag; n[2]=p[2]/mag;
}

void cross_product(vector<float>& p, const vector<float>& a, const vector<float>& b) {
    p[0] = a[1]*b[2] - a[2]*b[1];
    p[1] = a[2]*b[0] - a[0]*b[2];
    p[2] = a[0]*b[1] - a[1]*b[0];
}

void find4thPoint(vector<float>& p4,
	vector<float>& p1, vector<float>& p2, vector<float>& p3,
	float dist, float ang, float dihed) {
    vector<float> n1(3), n2(3), a(3), b(3);
    point_sub(a, p1,p2); point_sub(b, p3,p2);
    cross_product(n1, a,b); normalize_point(n1,n1);
    cross_product(n2, b,n1); normalize_point(n2,n2);

    double Sang = sin( DEGREES_TO_RADIANS(ang) );
    double Cang = cos( DEGREES_TO_RADIANS(ang) );
    double Sdihed = sin( DEGREES_TO_RADIANS(dihed) );
    double Cdihed = cos( DEGREES_TO_RADIANS(dihed) );
    normalize_point(b,b);
//cout << point_string(b) << point_string(n1) << point_string(n2) << endl;
    linear_combination(p4, Sang*Cdihed, n2, -1*Cang, b);
    linear_combination(p4, -1*Sang*Sdihed, n1, 1, p4);
//cout << point_string(p4) << endl;
    linear_combination(p4, dist, p4, 1, p3);
}

float calcDist(vector<float> & A, vector<float> & B) {
    vector<float> d(3);
    point_sub(d, A, B);
    return magnitude(d);
}

float calcAngle(vector<float> & A, vector<float> & B, vector<float> & C) {
    vector<float> BA(3), BC(3);
    linear_combination(BA, 1, A, -1, B);
    linear_combination(BC, 1, C, -1, B);
    normalize_point(BA, BA); normalize_point(BC, BC);
    float dp = dot_product(BA,BC);
    if(dp > 1) dp = 1;
    else if(dp < -1) dp = -1;
    return RADIANS_TO_DEGREES( acos(dp) );
}

float calcDihed(vector<float> &a, vector<float> &b, vector<float> &c, vector<float> &d) {
    //cout << "a " << a[0] <<" "<< a[1] <<" "<< a[2] << endl;
    //cout << "b " << b[0] <<" "<< b[1] <<" "<< b[2] << endl;
    //cout << "c " << c[0] <<" "<< c[1] <<" "<< c[2] << endl;
    //cout << "d " << d[0] <<" "<< d[1] <<" "<< d[2] << endl;

    vector<float> ba(3), bc(3), cb(3), cd(3);
    point_sub(ba, a, b); point_sub(bc, c, b);
    point_sub(cb, b, c); point_sub(cd, d, c);
    //cout << "ba " << ba[0] <<" "<< ba[1] <<" "<< ba[2] << endl;
    //cout << "bc " << bc[0] <<" "<< bc[1] <<" "<< bc[2] << endl;
    //cout << "cb " << cb[0] <<" "<< cb[1] <<" "<< cb[2] << endl;
    //cout << "cd " << cd[0] <<" "<< cd[1] <<" "<< cd[2] << endl;
    vector<float> ba_bc(3), cb_cd(3);
    cross_product(cb_cd, cb, cd); normalize_point(cb_cd, cb_cd);
    cross_product(ba_bc, ba, bc); normalize_point(ba_bc, ba_bc);
    //cout << "cb_cd " << cb_cd[0] <<" "<< cb_cd[1] <<" "<< cb_cd[2] << endl;
    //cout << "ba_bc " << ba_bc[0] <<" "<< ba_bc[1] <<" "<< ba_bc[2] << endl;
    float dp = dot_product(cb_cd, ba_bc);
    if(dp > 1) dp = 1;
    if(dp < -1) dp = -1;
    float angle = RADIANS_TO_DEGREES ( acos(dp) );
    //cout << angle <<" "<< dp << endl;
    vector<float> cp(3); cross_product(cp, ba_bc,cb_cd);
    if ( dot_product(cp,bc) < 0 ) angle = -1*angle;
    return angle;
}

// given a vector X, find normal vectors Y,Z such that X,Y,Z are mutually perpendicular
void findYZgivenX( vector<float> & givenX, vector<float> & Y, vector<float> & Z ) {
    vector<float> X(3);
    normalize_point(X, givenX);

    int a=0, b=1, c=2;
    if( fabs(X[a]) < 1e-2 ) { a=1; b=2; c=0; }
    if( fabs(X[a]) < 1e-2 ) { a=2; b=0; c=1; }
    assert ( fabs(X[a]) > 1e-2 );

    Y[b] = -1. + (2.*ran01()); Y[c] = -1. + (2.*ran01());
    //cout << "rands " << Y[b] << "\nrands " << Y[c] << endl;
    Y[a] = -1 * (X[b]*Y[b] + X[c]*Y[c]) / X[a];
    normalize_point(Y, Y);

    cross_product(Z, X, Y);
    normalize_point(Z, Z);
}

// on intersection of spheres c1,r1 and c2,r2 find p such that p-c1-q is of certain value, and lies on intersection of 2 spheres
// this is solved by sampling points on the circle of intersection and checking the angle
int findSphereSphereAngleIntx(float r1, vector<float> & c1, float r2, vector<float> & c2, vector<float> & q, float desiredAngle,
            vector<float>& p1, vector<float>& p2) {

    vector<float> c1c2(3), Y(3), Z(3), c(3);
    linear_combination(c1c2, 1, c2, -1, c1);
    double cc = magnitude(c1c2);
    cout << cc << endl;
    if(cc > r1+r2) return 0;

    normalize_point(c1c2, c1c2);
    findYZgivenX(c1c2, Y, Z);
    cout << "c1c2 " << c1c2[0] << " " << c1c2[1] << " " << c1c2[2] << endl;
    cout << "Y " << Y[0] << " " << Y[1] << " " << Y[2] << endl;
    cout << "Z " << Z[0] << " " << Z[1] << " " << Z[2] << endl;

    double c1c = (r1*r1 + cc*cc - r2*r2) / (2*cc);
    double r = sqrt(r1*r1 - c1c*c1c);
    linear_combination(c, 1, c1, c1c, c1c2);
    cout << "C " << c[0] << " " << c[1] << " " << c[2] << endl;

    float lastAngle = -999, angle, lastTheta = -999;
    float startTheta = 0, stopTheta = 360, step = 5; // clockwise

    int nsol = 0;
    vector<float> p(3);
    for(float theta = startTheta; theta <= stopTheta; theta += step) {
        double th = M_PI * theta / 180;
        linear_combination(p, 1, c, r*sin(th), Y);
        linear_combination(p, 1, p, r*cos(th), Z);
        angle = calcAngle(q, c1, p);
        cout << "P " << lastAngle << " " << desiredAngle << " " << angle << " " << p[0] << " " << p[1] << " " << p[2] << endl;
        if(lastAngle != -999 && (desiredAngle-angle)*(desiredAngle-lastAngle) <= 0) {
            double th = M_PI * (theta + lastTheta)/360.;
            if(nsol==0) { linear_combination(p1, 1, c, r*sin(th), Y); linear_combination(p1, 1, p1, r*cos(th), Z); }
            if(nsol==1) { linear_combination(p2, 1, c, r*sin(th), Y); linear_combination(p2, 1, p2, r*cos(th), Z); }
            assert(nsol!=2);
            nsol++;
            cout << "NSOL " << nsol << endl;
        }
        lastAngle = angle;
        lastTheta = theta;
    }
    return nsol;
}

float calcAngle(vector<float> & A, vector<float> & B) {
    vector<float> a = A, b = B;
    normalize_point(a, A); normalize_point(b, B);
    return RADIANS_TO_DEGREES(acos(dot_product(A,B)));
}

// as mentioned in AxisRotation.pdf in the same directory
void findRotnOperator(vector<vector<float> > & op, vector<float> & axis, float angle) {
    op.resize(3); op[0].resize(3); op[1].resize(3); op[2].resize(3);

    vector<float> n(3);
    normalize_point(n, axis);

    angle = DEGREES_TO_RADIANS(angle);
    float ct = cos(angle), st = sin(angle);

    op[0][0] = n[0]*n[0] + ct * (n[1]*n[1] + n[2]*n[2]);
    op[1][1] = n[1]*n[1] + ct * (n[2]*n[2] + n[0]*n[0]);
    op[2][2] = n[2]*n[2] + ct * (n[0]*n[0] + n[1]*n[1]);

    op[0][1] = n[0]*n[1]*(1-ct) - n[2]*st;
    op[1][0] = n[0]*n[1]*(1-ct) + n[2]*st;

    op[0][2] = n[0]*n[2]*(1-ct) + n[1]*st;
    op[2][0] = n[0]*n[2]*(1-ct) - n[1]*st;

    op[1][2] = n[1]*n[2]*(1-ct) - n[0]*st;
    op[2][1] = n[1]*n[2]*(1-ct) + n[0]*st;
}

void rotate(vector<vector<float> >& rotOp, vector<float>& pt) {
    vector<float> newpt(pt);
    newpt[0] = rotOp[0][0]*pt[0] + rotOp[0][1]*pt[1] + rotOp[0][2]*pt[2];
    newpt[1] = rotOp[1][0]*pt[0] + rotOp[1][1]*pt[1] + rotOp[1][2]*pt[2];
    newpt[2] = rotOp[2][0]*pt[0] + rotOp[2][1]*pt[1] + rotOp[2][2]*pt[2];
    pt = newpt;
}
