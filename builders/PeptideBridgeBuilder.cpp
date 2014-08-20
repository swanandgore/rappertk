#include <math.h>
#include "PeptideBridgeBuilder.h"
#include "geometry/functions.h"
#include "builders/Constants.h"
#include "samplers/RATdata.h"
#include "misc/RanGen.h"

int PeptideBridgeBuilder::cttrials = 25;
int PeptideBridgeBuilder::thetastep = 5;
void PeptideBridgeBuilder::setCTtrials(int ct) { cttrials = ct; }
void PeptideBridgeBuilder::setThetastep(int ts) { thetastep = ts; }

PeptideBridgeBuilder::PeptideBridgeBuilder(vector<int>& ipInds, vector<int>& opInds, Constants* con, const char* desc,
        const char *r0, const char *r1, const char *r2, const char *rdpath)
    : Builder(ipInds, opInds, con, desc) , RATdataPath(rdpath), resn0(r0), resn1(r1), resn2(r2) {

    CA_C = consts->get("CA_C");
    N_CA_C = consts->get("N_CA_C");
    C_N = consts->get("C_N");
    CA_C_N = consts->get("CA_C_N");
    C_O = consts->get("C_O");
    CA_C_O = consts->get("CA_C_O");
    TAU_QUALITY = consts->get("TAU_QUALITY");

    ctkeys.clear();
    ctkeys[0] = 0;
    ctkeys[1] = 0;
    ctkeys[10] = 0;
    ctkeys[11] = 0;

    transProp0 = 0.99, transProp1 = 0.99;
    if(resn1 == "PRO") transProp0 = 0.95;
    if(resn2 == "PRO") transProp1 = 0.95;

}

// expected inputs are :
// IP : C-(-1), N-0, CA-0, CA-2, C-2
// OP : C-0, O-0, N-1, CA-1, C-1, O-1, N-2
bool PeptideBridgeBuilder::build(VVF & pts) {
    sessionSize ++;
    vector<float> & c = pts[ip[0]];
    vector<float> & n0 = pts[ip[1]];
    vector<float> & ca0 = pts[ip[2]];

    vector<float> & c0 = pts[op[0]];
    vector<float> & o0 = pts[op[1]];
    vector<float> & n1 = pts[op[2]];
    vector<float> & ca1 = pts[op[3]];
    vector<float> & c1 = pts[op[4]];
    vector<float> & o1 = pts[op[5]];
    vector<float> & n2 = pts[op[6]];

    vector<float> & ca2 = pts[ip[3]];
    vector<float> & c2 = pts[ip[4]];

    vector<float> caca(3), cen(3);
    linear_combination(caca, 1, ca2, -1, ca0);
    float cadist = magnitude(caca);
    float interCA = 3.81, interCA_1 = 3.81;

    bool interCAfound = false;
    for(int i=0; i < 1; i++) {
        interCA = 3.81; interCA_1 = 3.81;
        if(ran01() > transProp0) interCA = 2.8; if(ran01() > transProp1) interCA_1 = 2.8;
        if( (cadist < interCA + interCA_1) && (checkAndAddIfNew(interCA, interCA_1)) ) { interCAfound = true; break; }
    }
    if( ! interCAfound ) { ca1 = c; return true; }

    normalize_point(caca, caca);
    vector<float> Y(3), Z(3);
    findYZgivenX(caca, Y, Z);
    float cad = (interCA*interCA - interCA_1*interCA_1 + cadist*cadist) / (2*cadist) ;
    linear_combination(cen, 1, ca0, cad, caca);
    float r = sqrt(interCA*interCA - cad*cad);

    float a1, t1, a2, t2, phi, psi, tau, caAng;
    vector<PSO>::iterator bi1, ei1, bi2, ei2;
    bool done = false;
    for(float theta=0; theta < 360; theta += thetastep) {
        float th = M_PI/180. * theta;
        linear_combination(ca1, 1, cen, r*sin(th), Y); // find CA-1 on the middle circle
        linear_combination(ca1, 1, ca1, r*cos(th), Z);
        caAng = calcAngle(ca0, ca1, ca2);
        if(caAng < 70 || caAng > 150) continue; // this angle has to be good
        RATdata & ratdata = RATdata::fwdinstance(RATdataPath.c_str());
        a1 = calcAngle(n0, ca0, ca1);
        t1 = calcDihed(c, n0, ca0, ca1);
        ratdata.range(resn0, interCA, a1, t1, bi1, ei1);
        for(; bi1 != ei1; bi1++) {
            phi = bi1->r; psi = bi1->a;
            find4thPoint( c0, c, n0, ca0, CA_C, N_CA_C, phi );
            find4thPoint( n1, n0, ca0, c0, C_N, CA_C_N, psi );
            find4thPoint( o0, n0, ca0, c0, C_O, CA_C_O, withinPlusMinus180(psi+180) );
            a2 = calcAngle(n1, ca1, ca2);
            t2 = calcDihed(c0, n1, ca1, ca2);
            ratdata.range(resn1, interCA_1, a2, t2, bi2, ei2);
            for(; bi2 != ei2; bi2++) {
                if(!checkAndAddIfNew(phi,psi,bi2->r,bi2->a,interCA,interCA_1)) continue;
                phi = bi2->r; psi = bi2->a;
                find4thPoint( c1, c0, n1, ca1, CA_C, N_CA_C, phi );
                find4thPoint( n2, n1, ca1, c1, C_N, CA_C_N, psi );
                find4thPoint( o1, n1, ca1, c1, C_O, CA_C_O, withinPlusMinus180(psi+180) );
                tau = calcAngle(n2, ca2, c2); // TAU
                if(fabs(N_CA_C-tau) < TAU_QUALITY) done = true; // TODO check for phi-psi around CA-2 also
                if(done) break;
            }
            if(done) break;
        }
        if(done) break;
    }
    if(!done) { // mess up output coordinates so that they get rejected in clash-check
        ca1 = c;
    }
    //if (interCA_1 < 3.8) cout << "INTERCA " << interCA <<" "<< interCA_1 << endl;

    // now build the CBeta's also
    buildCBs(pts);

    return true;
}


void PeptideBridgeBuilder::buildCBs(VVF & pts) {
    if(op.size() == 7) return; // no CBs asked to build
    float N_CA_CB = 0,
        CA_CB = consts->get("CA_CB"),
        C_N_CA_CB = consts->get("C_N_CA_CB");
    int cbi = 7;
    if(resn0 != "GLY") {
        if(resn0 == "PRO") N_CA_CB = consts->get("PRO_N_CA_CB");
        else               N_CA_CB = consts->get("N_CA_CB");
        find4thPoint( pts[ op[cbi] ],
            pts[ op[0] ], pts[ ip[1] ], pts[ ip[2] ],
            CA_CB, N_CA_CB, C_N_CA_CB );
        cbi++;
    }
    if(resn1 != "GLY") {
        if(resn1 == "PRO") N_CA_CB = consts->get("PRO_N_CA_CB");
        else               N_CA_CB = consts->get("N_CA_CB");
        find4thPoint( pts[ op[cbi] ],
            pts[ op[4] ], pts[ op[2] ], pts[ op[3] ],
            CA_CB, N_CA_CB, C_N_CA_CB );
        cbi++;
    }
    if(resn2 != "GLY") {
        if(resn2 == "PRO") N_CA_CB = consts->get("PRO_N_CA_CB");
        else               N_CA_CB = consts->get("N_CA_CB");
        find4thPoint( pts[ op[cbi] ],
            pts[ ip[4] ], pts[ op[6] ], pts[ ip[3] ],
            CA_CB, N_CA_CB, C_N_CA_CB );
        cbi++;
    }
}


bool PeptideBridgeBuilder::checkAndAddIfNew(float d1, float d2) {
    int thekey = 0;
    if(d1 < 3. && d2 < 3.) thekey = 0;
    if(d1 < 3. && d2 > 3.) thekey = 1;
    if(d1 > 3. && d2 < 3.) thekey = 10;
    if(d1 > 3. && d2 > 3.) thekey = 11;
    ctkeys[thekey] += 1;
    if(ctkeys[thekey] > cttrials) return false;
    return true;
}

bool PeptideBridgeBuilder::checkAndAddIfNew(float phi1, float psi1, float phi2, float psi2, float d1, float d2) {
    int a = 36, b = 5;
    long long int key = round(a+phi1/b) + round(a+psi1/b) * 72 + round(a+phi2/b) *72*72 + round(a+psi2/b) * 72*72*72;
    if(d1 > 3.) key += 72*72*72*72;
    if(d2 > 3.) key += 72*72*72*72*2;
    if(keys.find(key) != keys.end()) return false;
    keys.insert(key);
    return true;
}

void PeptideBridgeBuilder::deleteSessionInfo() {
    keys.clear();
    ctkeys[0] = 0;
    ctkeys[1] = 0;
    ctkeys[10] = 0;
    ctkeys[11] = 0;
}

PeptideBridgeBuilder PeptideBridgeBuilder::makeCopy() {
    return PeptideBridgeBuilder(*this);
}
