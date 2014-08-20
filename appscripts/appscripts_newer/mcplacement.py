import os, sys
cbDatapath = os.environ["RTKROOT"] + "/data/"
from data import resAtoms, consts
from samplers import PhipsiSampler, OmegaSampler
from builders import PeptideBuilder, VecVecFloat, VecInt, VecFloat
from geometry import calcDist, calcAngle, calcDihed

if __name__ == "__main__" :
    IPs, OPs = [], []

    res, resns, resnums, chids, inscodes = {}, {}, {}, {}, {}
    pts = []

    resns[0] = "THR" ; resnums[0] = "  92" ; chids[0] = "A" ; inscodes[0] = ' ' ; res[0] = {}
    pts.append( [1.715, -7.982, 5.876] ) ; res[0][' C  '] = len(pts)-1

    resns[1] = sys.argv[1] ; resnums[1] = "  93" ; chids[1] = "A" ; inscodes[1] = ' ' ; res[1] = {}
    pts.append([1.237, -6.819, 6.310]) ;  res[1][' N  '] = len(pts)-1
    pts.append([-0.163, -6.671, 6.702]) ; res[1][' CA '] = len(pts)-1
    pts.append([-0.803, -5.457, 6.037]) ; res[1][' C  '] = len(pts)-1
    pts.append([-0.846, -4.378, 6.621]) ; res[1][' O  '] = len(pts)-1
    print 'here'
    ps = PhipsiSampler( "%s/PhipsiWeightedProp/ps%s" % (cbDatapath,resns[1]) )
    print 'here'
    omegaSampler = OmegaSampler()
    print 'here'

    resns[2] = "ALA" ; resnums[2] = "  94" ; chids[2] = "A" ; inscodes[2] = ' ' ; res[2] = {}
    pts.append( [-1.312, -5.638, 4.822] ) ; res[2][' N  '] = len(pts)-1
    pts.append( [-1.312, -5.638, 4.822] ) ; res[2][' CA '] = len(pts)-1

    IPs = [ res[0][' C  '], res[1][' N  '], res[1][' CA '] ]
    OPs = [ res[1][' C  '], res[1][' O  '], res[2][' N  '], res[2][' CA '] ]

    pts = VecVecFloat(pts)
    IPs, OPs = VecInt(IPs), VecInt(OPs)
    print 'here'
    pb = PeptideBuilder(pts, IPs, OPs, consts, "PEPBUILDER", ps, omegaSampler, 1, 1)
    print 'here'

    from peptidebuild import ModelRenderer
    mr = ModelRenderer(res, resns, chids, resnums, inscodes, pts)

    for phi in range(-180, 180, 5) :
        for psi in range(-180, 180, 5) :
            if ps.findProb(phi,psi) < 1e-10 : continue
            print phi, psi
            for omega in (0,-180) :
                pb.build(phi, psi, omega)
                #mr.render()
                r = calcDist( VecFloat(pts[res[1][' CA ']]), VecFloat(pts[res[2][' CA ']]) )
                a = calcAngle( VecFloat(pts[res[1][' N  ']]), VecFloat(pts[res[1][' CA ']]), VecFloat(pts[res[2][' CA ']]) )
                t = calcDihed( VecFloat(pts[res[0][' C  ']]), VecFloat(pts[res[1][' N  ']]), VecFloat(pts[res[1][' CA ']]), VecFloat(pts[res[2][' CA ']]) )
                print "MAP", resns[1], int(phi), int(psi), int(omega), int(round(r*10)), int(round(a/5.)*5), int(round(t/5.)*5)
