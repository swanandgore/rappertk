#include <math.h>
#include <assert.h>

#include "ChiBuilder.h"

#include "samplers/BBdepChiSampler.h"
#include "geometry/functions.h"
#include "misc/verbosity.h"
#include "Constants.h"

#include<cstring>
#include<stdio.h>
#include<cstdlib>
#include<math.h>

ChiBuilder::ChiBuilder(vector<int>& ipInds, vector<int>& opInds,
    Constants *consts, const char* desc, const char *rn, BBdepChiSampler *cs)
        : Builder(ipInds, opInds, consts, desc) {
    resname = string(rn);
    chisampler = cs;
}

ChiBuilder ChiBuilder::makeCopy() { return ChiBuilder(*this); }

bool ChiBuilder::checkAndAddIfNew(int ci) {
    if(sessionSize >= maxSessionSize) return false;
    sessionChis.insert(ci);
    if(sessionSize == sessionChis.size()) return false;
    sessionSize ++; return true;
}

void ChiBuilder::deleteSessionInfo() { sessionChis.clear(); }

bool ChiBuilder::buildSample(VVF & pts, int sampleIndex) {
    int iCp = ip[0], iN = ip[1], iCA = ip[2], iC = ip[3], iNn = ip[4];
    float phi = calcDihed(pts[iCp], pts[iN], pts[iCA], pts[iC]);
    float psi = calcDihed(pts[iN], pts[iCA], pts[iC], pts[iNn]);

    float sample[4]; sample[0]=-9999; sample[1]=-9999; sample[2]=-9999; sample[3]=-9999;
    int ci = chisampler->sample(phi, psi, sample, sampleIndex);
    if(ci != sampleIndex) return false;

    //    if(verbose(7)) cout << "phipsi " << phi <<" "<< psi << " chi "
    //   << sample[0] <<" "<< sample[1] <<" "<< sample[2] <<" "<< sample[3] << endl;

    if(session && !checkAndAddIfNew(ci)) {
      //if(verbose(7)) cout << "chis already sampled" << endl;
        return false; // sampled in this session already
    }
    return build(pts, sample);
}

// ip : C, N, CA, C, N, CB
bool ChiBuilder::build(VVF & pts) {
    int iCp = ip[0], iN = ip[1], iCA = ip[2], iC = ip[3], iNn = ip[4];
    float phi = calcDihed(pts[iCp], pts[iN], pts[iCA], pts[iC]);
    float psi = calcDihed(pts[iN], pts[iCA], pts[iC], pts[iNn]);

    float sample[4]; sample[0]=-9999; sample[1]=-9999; sample[2]=-9999; sample[3]=-9999;
    int ci = chisampler->sample(phi, psi, sample);

    //    if(verbose(7)) cout << "phipsi " << phi <<" "<< psi << " chi "
    //   << sample[0] <<" "<< sample[1] <<" "<< sample[2] <<" "<< sample[3] << endl;

    if(session && !checkAndAddIfNew(ci)) {
      //if(verbose(7)) cout << "chis already sampled" << endl;
        return false; // sampled in this session already
    }
    return build(pts, sample);
}

bool ChiBuilder::build(VVF & pts, float *sample) {
    if(resname == "VAL") buildSC_VAL(pts, sample);
    else if(resname == "ARG") buildSC_ARG(pts, sample);
    else if(resname == "LYS") buildSC_LYS(pts, sample);
    else if(resname == "GLU") buildSC_GLU(pts, sample);
    else if(resname == "ILE") buildSC_ILE(pts, sample);
    else if(resname == "ASN") buildSC_ASN(pts, sample);
    else if(resname == "LEU") buildSC_LEU(pts, sample);
    else if(resname == "HIS") buildSC_HIS(pts, sample);
    else if(resname == "THR") buildSC_THR(pts, sample);
    else if(resname == "MET") buildSC_MET(pts, sample);
    else if(resname == "MSE") buildSC_MET(pts, sample);
    else if(resname == "GLN") buildSC_GLN(pts, sample);
    else if(resname == "ASP") buildSC_ASP(pts, sample);
    else if(resname=="PRO") buildSC_PRO(pts, sample);
    else if(resname=="SER") buildSC_SER(pts, sample);
    else if(resname=="PHE") buildSC_PHE(pts, sample);
    else if(resname=="TYR") buildSC_TYR(pts, sample);
    else if(resname=="CYS") buildSC_CYS(pts, sample);
    else if(resname=="TRP") buildSC_TRP(pts, sample);
    else {
        cerr << "No function to handle sidechain type " << resname << endl;
        exit(1);
    }
    return true;
}

// ip : C, N, CA, C, N, CB
// op : CG, CD, NE, CZ, NH1, NH2
void ChiBuilder::buildSC_ARG(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& cd = pts[ op[1] ];
    vector<float>& ne = pts[ op[2] ];
    vector<float>& cz = pts[ op[3] ];
    vector<float>& nh1 = pts[ op[4] ];
    vector<float>& nh2 = pts[ op[5] ];

    find4thPoint(cg, n, ca, cb, cget("ARG_CB_CG"), cget("ARG_CA_CB_CG"), chis[0]);
    find4thPoint(cd, ca, cb, cg, cget("ARG_CG_CD"), cget("ARG_CB_CG_CD"), chis[1]);
    find4thPoint(ne, cb, cg, cd, cget("ARG_CD_NE"), cget("ARG_CG_CD_NE"), chis[2]);
    find4thPoint(cz, cg, cd, ne, cget("ARG_NE_CZ"), cget("ARG_CD_NE_CZ"), chis[3]);
    find4thPoint(nh1, cd, ne, cz, cget("ARG_CZ_NH1"), cget("ARG_NE_CZ_NH1"), 0.);
    find4thPoint(nh2, cd, ne, cz, cget("ARG_CZ_NH2"), cget("ARG_NE_CZ_NH2"), -180.);
}

// ip : C, N, CA, C, N, CB
// op : CG1, CG2
void ChiBuilder::buildSC_VAL(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg1 = pts[ op[0] ];
    vector<float>& cg2 = pts[ op[1] ];

    float chi1 = chis[0];
    find4thPoint(cg1, n, ca, cb, cget("VAL_CB_CG1"), cget("VAL_CA_CB_CG1"), chi1);

    chi1 = withinPlusMinus180(chi1+120);
    find4thPoint(cg2, n, ca, cb, cget("VAL_CB_CG2"), cget("VAL_CA_CB_CG2"), chi1);
}

// ip : C, N, CA, C, N, CB
// op : CG, CD, CE, NZ
void ChiBuilder::buildSC_LYS(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& cd = pts[ op[1] ];
    vector<float>& ce = pts[ op[2] ];
    vector<float>& nz = pts[ op[3] ];

	find4thPoint(cg, n, ca, cb, cget("LYS_CB_CG"), cget("LYS_CA_CB_CG"), chis[0]);
	find4thPoint(cd, ca, cb, cg, cget("LYS_CG_CD"), cget("LYS_CB_CG_CD"), chis[1]);
	find4thPoint(ce, cb, cg, cd, cget("LYS_CD_CE"), cget("LYS_CG_CD_CE"), chis[2]);
	find4thPoint(nz, cg, cd, ce, cget("LYS_CE_NZ"), cget("LYS_CD_CE_NZ"), chis[3]);
}

// ip : C, N, CA, C, N, CB
// op : CG, CD, OE1, OE2
void ChiBuilder::buildSC_GLU(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& cd = pts[ op[1] ];
    vector<float>& oe1 = pts[ op[2] ];
    vector<float>& oe2 = pts[ op[3] ];

	find4thPoint(cg, n, ca, cb, cget("GLU_CB_CG"), cget("GLU_CA_CB_CG"), chis[0]);
	find4thPoint(cd, ca, cb, cg, cget("GLU_CG_CD"), cget("GLU_CB_CG_CD"), chis[1]);
	find4thPoint(oe1, cb, cg, cd, cget("GLU_CD_OE1"), cget("GLU_CG_CD_OE1"), chis[2]);
	float chi2 = withinPlusMinus180(chis[2]+180);
	find4thPoint(oe2, cb, cg, cd, cget("GLU_CD_OE2"), cget("GLU_CG_CD_OE2"), chi2);
}

// ip : C, N, CA, C, N, CB
// op : CG1, CG2, CD1
void ChiBuilder::buildSC_ILE(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg1 = pts[ op[0] ];
    vector<float>& cg2 = pts[ op[1] ];
    vector<float>& cd1 = pts[ op[2] ];

	find4thPoint(cg1, n, ca, cb, cget("ILE_CB_CG1"), cget("ILE_CA_CB_CG1"), chis[0]);
	float chi0 = withinPlusMinus180(chis[0]-120);
	find4thPoint(cg2, n, ca, cb, cget("ILE_CB_CG2"), cget("ILE_CA_CB_CG2"), chi0);
	find4thPoint(cd1, ca, cb, cg1, cget("ILE_CG1_CD1"), cget("ILE_CB_CG1_CD1"), chis[1]);
}

// ip : C, N, CA, C, N, CB
// op : CG, OD1, ND2
void ChiBuilder::buildSC_ASN(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& od1 = pts[ op[1] ];
    vector<float>& nd2 = pts[ op[2] ];

	find4thPoint(cg, n, ca, cb, cget("ASN_CB_CG"), cget("ASN_CA_CB_CG"), chis[0]);
	find4thPoint(od1, ca, cb, cg, cget("ASN_CG_OD1"), cget("ASN_CB_CG_OD1"), chis[1]);
	float chi1 = withinPlusMinus180(chis[1]+180);
	find4thPoint(nd2, ca, cb, cg, cget("ASN_CG_ND2"), cget("ASN_CB_CG_ND2"), chi1);
}

// ip : C, N, CA, C, N, CB
// op : CG, CD1, CD2
void ChiBuilder::buildSC_LEU(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& cd1 = pts[ op[1] ];
    vector<float>& cd2 = pts[ op[2] ];

	find4thPoint(cg, n, ca, cb, cget("LEU_CB_CG"), cget("LEU_CA_CB_CG"), chis[0]);
	find4thPoint(cd1, ca, cb, cg, cget("LEU_CG_CD1"), cget("LEU_CB_CG_CD1"), chis[1]);
	float chi1 = withinPlusMinus180(chis[1]+120);
	find4thPoint(cd2, ca, cb, cg, cget("LEU_CG_CD2"), cget("LEU_CB_CG_CD2"), chi1);
}

// ip : C, N, CA, C, N, CB
// op : CG, CD1, CD2
void ChiBuilder::buildSC_HIS(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& nd1 = pts[ op[1] ];
    vector<float>& cd2 = pts[ op[2] ];
    vector<float>& ce1 = pts[ op[3] ];
    vector<float>& ne2 = pts[ op[4] ];

	find4thPoint(cg, n, ca, cb, cget("HIS_CB_CG"), cget("HIS_CA_CB_CG"), chis[0]);
	find4thPoint(nd1, ca, cb, cg, cget("HIS_CG_ND1"), cget("HIS_CB_CG_ND1"), chis[1]);
	find4thPoint(ce1, cb, cg, nd1, cget("HIS_ND1_CE1"), cget("HIS_CG_ND1_CE1"), 180.);
	find4thPoint(ne2, cg, nd1, ce1, cget("HIS_CE1_NE2"), cget("HIS_ND1_CE1_NE2"), 0.);
	find4thPoint(cd2, nd1, ce1, ne2, cget("HIS_NE2_CD2"), cget("HIS_CE1_NE2_CD2"), 0.);
}


// ip : C, N, CA, C, N, CB
// op : OG1, CG2
void ChiBuilder::buildSC_THR(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& og1 = pts[ op[0] ];
    vector<float>& cg2 = pts[ op[1] ];

    find4thPoint(og1, n, ca, cb, cget("THR_CB_OG1"), cget("THR_CA_CB_OG1"), chis[0]);
    float chi1 = withinPlusMinus180(chis[0]-120);
    find4thPoint(cg2, n, ca, cb, cget("THR_CB_CG2"), cget("THR_CA_CB_CG2"), chi1);
}

// ip : C, N, CA, C, N, CB
// op : CG, SD, CE
void ChiBuilder::buildSC_MET(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& sd = pts[ op[1] ];
    vector<float>& ce = pts[ op[2] ];

	find4thPoint(cg, n, ca, cb, cget("MET_CB_CG"), cget("MET_CA_CB_CG"), chis[0]);
	find4thPoint(sd, ca, cb, cg, cget("MET_CG_SD"), cget("MET_CB_CG_SD"), chis[1]);
	find4thPoint(ce, cb, cg, sd, cget("MET_SD_CE"), cget("MET_CG_SD_CE"), chis[2]);
}

// ip : C, N, CA, C, N, CB
// op : CG, CD, OE1, NE2
void ChiBuilder::buildSC_GLN(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& cd = pts[ op[1] ];
    vector<float>& oe1 = pts[ op[2] ];
    vector<float>& ne2 = pts[ op[3] ];

	find4thPoint(cg, n, ca, cb, cget("GLN_CB_CG"), cget("GLN_CA_CB_CG"), chis[0]);
	find4thPoint(cd, ca, cb, cg, cget("GLN_CG_CD"), cget("GLN_CB_CG_CD"), chis[1]);
	find4thPoint(oe1, cb, cg, cd, cget("GLN_CD_OE1"), cget("GLN_CG_CD_OE1"), chis[2]);
	float chi3 = withinPlusMinus180(chis[2]+180);
	find4thPoint(ne2, cb, cg, cd, cget("GLN_CD_NE2"), cget("GLN_CG_CD_NE2"), chi3);
}

// ip : C, N, CA, C, N, CB
// op : CG, OD1, OD2
void ChiBuilder::buildSC_ASP(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& od1 = pts[ op[1] ];
    vector<float>& od2 = pts[ op[2] ];

	find4thPoint(cg, n, ca, cb, cget("ASP_CB_CG"), cget("ASP_CA_CB_CG"), chis[0]);
	find4thPoint(od1, ca, cb, cg, cget("ASP_CG_OD1"), cget("ASP_CB_CG_OD1"), chis[1]);
	float chi2 = withinPlusMinus180(chis[1]+180);
	find4thPoint(od2, ca, cb, cg, cget("ASP_CG_OD2"), cget("ASP_CB_CG_OD2"), chi2);
}

void ChiBuilder::buildSC_PRO(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& cd = pts[ op[1] ];

	find4thPoint(cg, n, ca, cb, cget("PRO_CB_CG"), cget("PRO_CA_CB_CG"), chis[0]);
	find4thPoint(cd, ca, cb, cg, cget("PRO_CG_CD"), cget("PRO_CB_CG_CD"), chis[1]);
}
void ChiBuilder::buildSC_SER(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& og = pts[ op[0] ];

	find4thPoint(og, n, ca, cb, cget("SER_CB_OG"), cget("SER_CA_CB_OG"), chis[0]);
}
void ChiBuilder::buildSC_PHE(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& cd1 = pts[ op[1] ];
    vector<float>& cd2 = pts[ op[2] ];
    vector<float>& ce1 = pts[ op[3] ];
    vector<float>& ce2 = pts[ op[4] ];
    vector<float>& cz = pts[ op[5] ];

	find4thPoint(cg, n, ca, cb, cget("PHE_CB_CG"), cget("PHE_CA_CB_CG"), chis[0]);
	find4thPoint(cd1, ca, cb, cg, cget("PHE_CG_CD1"), cget("PHE_CB_CG_CD1"), chis[1]);
	float chi2 = withinPlusMinus180(chis[1]+180);
	find4thPoint(cd2, ca, cb, cg, cget("PHE_CG_CD2"), cget("PHE_CB_CG_CD2"), chi2);
	find4thPoint(ce1, cb, cg, cd1, cget("PHE_CD1_CE1"), cget("PHE_CG_CD1_CE1"), -180.);
	find4thPoint(ce2, cb, cg, cd2, cget("PHE_CD2_CE2"), cget("PHE_CG_CD2_CE2"), -180.);
	find4thPoint(cz, cg, cd1, ce1, cget("PHE_CE1_CZ"), cget("PHE_CD1_CE1_CZ"), 0.);
}

void ChiBuilder::buildSC_TYR(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& cd1 = pts[ op[1] ];
    vector<float>& cd2 = pts[ op[2] ];
    vector<float>& ce1 = pts[ op[3] ];
    vector<float>& ce2 = pts[ op[4] ];
    vector<float>& cz = pts[ op[5] ];
    vector<float>& oh = pts[ op[6] ];

	find4thPoint(cg, n, ca, cb, cget("TYR_CB_CG"), cget("TYR_CA_CB_CG"), chis[0]);
	find4thPoint(cd1, ca, cb, cg, cget("TYR_CG_CD1"), cget("TYR_CB_CG_CD1"), chis[1]);
	float chi2 = withinPlusMinus180(chis[1]+180);
	find4thPoint(cd2, ca, cb, cg, cget("TYR_CG_CD2"), cget("TYR_CB_CG_CD2"), chi2);
	find4thPoint(ce1, cb, cg, cd1, cget("TYR_CD1_CE1"), cget("TYR_CG_CD1_CE1"), -180.);
	find4thPoint(ce2, cb, cg, cd2, cget("TYR_CD2_CE2"), cget("TYR_CG_CD2_CE2"), -180.);
	find4thPoint(cz, cg, cd1, ce1, cget("TYR_CE1_CZ"), cget("TYR_CD1_CE1_CZ"), 0.);
	find4thPoint(oh, ce1, ce2, cz, cget("TYR_CE1_CZ"), 120, -180.);
}

void ChiBuilder::buildSC_CYS(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& sg = pts[ op[0] ];

	find4thPoint(sg, n, ca, cb, cget("CYS_CB_SG"), cget("CYS_CA_CB_SG"), chis[0]);
}
void ChiBuilder::buildSC_TRP(VVF & pts, float *chis) {
    vector<float>& n  = pts[ ip[1] ];
    vector<float>& ca = pts[ ip[2] ];
    vector<float>& cb = pts[ ip[5] ];

    vector<float>& cg = pts[ op[0] ];
    vector<float>& cd1 = pts[ op[1] ];
    vector<float>& cd2 = pts[ op[2] ];
    vector<float>& ne1 = pts[ op[3] ];
    vector<float>& ce2 = pts[ op[4] ];
    vector<float>& ce3 = pts[ op[5] ];
    vector<float>& cz2 = pts[ op[6] ];
    vector<float>& cz3 = pts[ op[7] ];
    vector<float>& ch2 = pts[ op[8] ];

	find4thPoint(cg, n, ca, cb, cget("TRP_CB_CG"), cget("TRP_CA_CB_CG"), chis[0]);
	find4thPoint(cd1, ca, cb, cg, cget("TRP_CG_CD1"), cget("TRP_CB_CG_CD1"), chis[1]);
	float chi2 = withinPlusMinus180(chis[1]+180);
	find4thPoint(cd2, ca, cb, cg, cget("TRP_CG_CD2"), cget("TRP_CB_CG_CD2"), chi2);
	find4thPoint(ne1, cb, cg, cd1, cget("TRP_CD1_NE1"), cget("TRP_CG_CD1_NE1"), -180.);
	find4thPoint(ce2, cb, cg, cd2, cget("TRP_CD2_CE2"), cget("TRP_CG_CD2_CE2"), -180.);
	find4thPoint(cz2, cg, cd2, ce2, cget("TRP_CE2_CZ2"), cget("TRP_CD2_CE2_CZ2"), -180.);
	find4thPoint(ch2, cd2, ce2, cz2, cget("TRP_CZ2_CH2"), cget("TRP_CE2_CZ2_CH2"), 0.);
	find4thPoint(ce3, cz2, ce2, cd2, cget("TRP_CD2_CE3"), cget("TRP_CE2_CD2_CE3"), 0.);
	find4thPoint(cz3, ce2, cd2, ce3, cget("TRP_CE3_CZ3"), cget("TRP_CD2_CE3_CZ3"), 0.);
}
