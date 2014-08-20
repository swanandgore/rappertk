import string, math, sys

from geometry import calcDist, calcAngle, calcDihed
import data
from builders import PeptideBuilder, VecInt, VecFloat, VecVecFloat
from data import consts

def rad2deg(rad) : return 180./math.pi * rad
def deg2rad(deg) : return math.pi/180. * deg

def initCNCA(pts) :
    pts[0] = [0.,0.,0.]
    cn = consts.get("C_N")
    nca = consts.get("N_CA")
    cnca = consts.get("C_N_CA")
    pts[1] = [cn, 0., 0.,]
    pts[2] = [cn + nca * math.cos(deg2rad(180-cnca)), nca * math.sin(deg2rad(cnca)), 0.]

def initNCCA(pts) :
    pts[0] = [0.,0.,0.]
    cn = consts.get("C_N")
    nca = consts.get("CA_C")
    cnca = consts.get("CA_C_N")
    pts[1] = [cn, 0., 0.,]
    pts[2] = [cn + nca * math.cos(deg2rad(180-cnca)), nca * math.sin(deg2rad(cnca)), 0.]

def createRAT(resPhipsi, fwdfile, bwdfile) :
    fpts, bpts, IPs, OPs = [], [], [0,1,2], [3,4,5,6] # C,N,CA -> C,O,N,CA or N,C,CA -> N,C,O,CA
    for i in range(7) :
        fpts.append([0.,0.,0.])
        bpts.append([0.,0.,0.])
    initCNCA(fpts) ; initNCCA(bpts)
    fpts, bpts = VecVecFloat(fpts), VecVecFloat(bpts)
    IPs,OPs = VecInt(IPs), VecInt(OPs)
    fpb = PeptideBuilder( IPs, OPs, consts, "ALA", "ALAbuilder", None, None, 1 )
    bpb = PeptideBuilder( IPs, OPs, consts, "ALA", "ALAbuilder", None, None, 0 )
    ff = open(fwdfile, 'w')
    fb = open(bwdfile, 'w')
    for resn in resPhipsi.keys() :
        for phi,psi in resPhipsi[resn] :
            for omega in (-180, 0) :
                fpb.build1(fpts, phi, psi, omega)
                r = calcDist( VecFloat(fpts[2]), VecFloat(fpts[6]) )
                a = calcAngle( VecFloat(fpts[1]), VecFloat(fpts[2]), VecFloat(fpts[6]) )
                t = calcDihed( VecFloat(fpts[0]), VecFloat(fpts[1]), VecFloat(fpts[2]), VecFloat(fpts[6]) )
                print >> ff, resn, phi, psi, omega, int(round(r*10)), int(round(a/5.)*5), int(round(t/5.)*5)
                bpb.build1(bpts, phi, psi, omega)
                r = calcDist( VecFloat(bpts[2]), VecFloat(bpts[6]) )
                a = calcAngle( VecFloat(bpts[1]), VecFloat(bpts[2]), VecFloat(bpts[6]) )
                t = calcDihed( VecFloat(bpts[0]), VecFloat(bpts[1]), VecFloat(bpts[2]), VecFloat(bpts[6]) )
                print >> fb, resn, psi, phi, omega, int(round(r*10)), int(round(a/5.)*5), int(round(t/5.)*5)
    ff.close()
    fb.close()


import os
cbDatapath = os.environ["RTKROOT"] + "/data/"

if __name__ == "__main__" :
    lines = open(cbDatapath + "/PhipsiWeightedProp/top500.coil.count").readlines()
    for level in [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001] :
        resPhipsi = {}
        resns = lines[0].split()[2:]
        for r in resns : resPhipsi[r] = []
        for l in lines[1:] :
            flds = l.split()
            ps = (string.atoi(flds[0]), string.atoi(flds[1]))
            for i in range(2,len(flds)) :
                if string.atof(flds[i]) > level : resPhipsi[ resns[i-2] ].append(ps)
        ff = cbDatapath + "/rat/fwd.%f" % level
        fb = cbDatapath + "/rat/bwd.%f" % level
        createRAT(resPhipsi, ff, fb)
    import sys ; sys.exit(0)

    ratmap, revmap = {}, {}
    for l in open("ratmap", 'r').readlines() :
        flds = l.split()
        phi, psi, omega, r, a, t = \
            string.atof(flds[1]), string.atof(flds[2]), string.atof(flds[3]), \
            string.atof(flds[4]), string.atof(flds[5]), string.atof(flds[6]),
        if omega > -1 and omega < 1 : continue
        rn = round(r * 10)
        an = round( a/5. ) * 5
        tn = round( t/5. ) * 5

        phi = int(phi)
        psi = int(psi)
        omega = int(omega)
        rn = int(rn)
        an = int(an)
        tn = int(tn)
        #print "%d %d %d %3.2f %d %d" % (phi, psi, omega, rn, an, tn)
        #print phi, psi, omega, rn, an, tn
        if not (phi,psi,omega) in ratmap.keys() : ratmap[(phi,psi,omega)] = []
        ratmap[(phi,psi,omega)].append( (rn,an,tn) )
        if not (rn,an,tn) in revmap.keys() : revmap[(rn,an,tn)] = []
        revmap[(rn,an,tn)].append( (phi,psi,omega) )

    maxlen = -100
    for k,v in ratmap.items() :
        if maxlen < len(v) : maxlen = len(v)
    print "maxrat", maxlen
    maxlen = -100
    vlens = []
    for k,v in revmap.items() :
        if maxlen < len(v) : maxlen = len(v)
        vlens.append( len(v) )
    print "maxrev", maxlen
    vlens.sort()
    print vlens
        
