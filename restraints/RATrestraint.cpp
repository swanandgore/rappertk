#include "RATrestraint.h"
#include "../samplers/RATdata.h"
#include "../geometry/functions.h"
#include "../misc/RanGen.h"

bool notInSphere(vector<float> & p, vector<float> & cen, float rad) {
    float dx = fabs(p[0] - cen[0]);
    if(dx > rad) return true;
    float dy = fabs(p[1] - cen[1]);
    if(dy > rad) return true;
    float dz = fabs(p[2] - cen[2]);
    if(dz > rad) return true;

    if(dx*dx + dy*dy + dz*dz > rad*rad) return true;

    return false;
}

RATrestraint::RATrestraint(vector<int> & pis, const char *desc, vector<float> & p, float r, const char *resname, RATdata* rd)
    : Restraint(pis, desc), cen(p), rad(r), resn(resname), ratdata(rd)
{}

void sphereVolumeSample(vector<float> & cen, float rad, vector<float> & sample) {
    
}

#ifdef DONT_DEFINE_ME
bool RATrestraint::satisfied(VVF & pts) {
    vector<float> & C  = pts[ ptInds[0] ];
    vector<float> & N  = pts[ ptInds[1] ];
    vector<float> & CA = pts[ ptInds[2] ];

    // find the circle presented by the sphere to CA, and iterate square grid on it
    float cacaDist = calcDist(CA,cen);
    float alpha = asin(rad/cacaDist);
    float tangent = sqrt(cacaDist*cacaDist - rad*rad);
    float r = tangent * sin(alpha);
    float d = tangent * cos(alpha);

    vector<float> CAcen(3,0), c(3,0), X(3), Y(3), pt(3), v(3);
    linear_combination(CAcen, 1, cen, -1, CA);
    normalize_point(CAcen, CAcen);
    findYZgivenX(CAcen, X, Y);
    linear_combination(c, 1, CA, d, CAcen);

    vector<PSO>::iterator bi, ei;
    float rstep = r / 5., interCAdist, a, t;
    if(rstep > 0.2) rstep = 0.2; // rstep is max 0.2
    int numPSO = 0, minPSO = 10;
    for(float ri = -r; ri <= r; ri += rstep)
        for(float rj = -r; rj <= r; rj += rstep) {
            linear_combination(pt, 1,  c, ri, X);
            linear_combination(pt, 1, pt, rj, Y);
            a = calcAngle(N, CA, pt);
            t = calcDihed(C, N, CA, pt);
            interCAdist = 3.8;
            find4thPoint(v, C, N, CA, interCAdist, a, t);
            if(notInSphere(v, cen, rad)) {
                interCAdist = 2.8;
                find4thPoint(v, C, N, CA, interCAdist, a, t);
                if(notInSphere(v, cen, rad)) continue;
            }
            ratdata->range(resn, interCAdist, a, t, bi, ei);
            numPSO += (ei - bi);
            if(numPSO >= minPSO) return true;
        }
    return false;
}
#endif // DONT_DEFINE_ME

bool RATrestraint::satisfied(VVF & pts) {
    vector<float> & C  = pts[ ptInds[0] ];
    vector<float> & N  = pts[ ptInds[1] ];
    vector<float> & CA = pts[ ptInds[2] ];
    // take 50 samples in the sphere cen,rad.
    // satisfied if at least minPSO phi-psi-omega triples are found
    vector<float> noi(3,0), ca(3,0), CAca(3,0);
    vector<PSO>::iterator bi, ei;
    float interCAdist, a, t;
    int numPSO = 0, minPSO = 10;
    for(int si = 0; si < 500; si++) {
        randomNormalVector(noi);
        linear_combination(ca, 1, cen, ran01() * rad, noi);
        linear_combination(CAca, 1, ca, -1, CA);
        normalize_point(CAca, CAca);
        interCAdist = 3.8;
        linear_combination(ca, 1, CA, interCAdist, CAca);
        if(notInSphere(ca, cen, rad)) {
            linear_combination(ca, 1, CA, interCAdist, CAca);
            if(notInSphere(ca, cen, rad)) continue;
        }
        a = calcAngle(N, CA, ca);
        t = calcDihed(C, N, CA, ca);
//cout << "rat " << r <<" "<< a <<" "<< t << endl;
        ratdata->range(resn, interCAdist, a, t, bi, ei);
        numPSO += (ei - bi);
        if(numPSO >= minPSO) return true;
    }
    //cout << "found numPSO " << numPSO << endl;
    return false;
}
