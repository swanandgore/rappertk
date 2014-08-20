from data import resAtoms, nuConn
from pdbr import protein, line2resid, line2resn, line2atomname, line2crd, line2chid, line2resnum, line2resic
from samplers import NAsuiteSampler
from builders import NAsuiteBuilder, VecInt, VecVecFloat, VecFloat
from data import consts, vdwr
from strategy import makeClashExcl_1st2ndCovNbr, strategy
from geometry import MapIntIntFloat
from restraints import SphPosRestr
import prot2res

import os
cbDatapath = os.environ["RTKROOT"] + "/data/"

def nuMain() :
    blist, numtrials = [], []
    radii, knownPositions = [], []
    clashExclInds = []

    prot = protein('/tmp/nu.pdb', read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    radii = [0] * len(pts)
    # fill radii & ensure that all atomnames in res are indeed correct
    for index,n2i in res.items() :
        for aname,ai in n2i.items() :
            print resids[index], aname, ai
            assert aname in resAtoms[ resns[index] ]
            radii[ai] = vdwr['XXX'][aname[1]]

    covconn = {}
    for i in range(len(pts)) : covconn[i] = []
    keys = res.keys() ; keys.sort()
    for k in keys :
        for c in nuConn[ resns[k] ] :
            covconn[ res[k][c[0]] ].append( res[k][c[1]] ) ; covconn[ res[k][c[1]] ].append( res[k][c[0]] )
    for i in range( len(keys)-1 ) :
        ki, kj = keys[i], keys[i+1]
        n1, n2 = res[ki][' O3*'], res[kj][' P  ']
        covconn[n1].append(n2) ; covconn[n2].append(n1)

    # put 0th residue's P, O1P, O2P, O5*, C5*, C4*, C3* as known, and also some of last residue atoms
    names = [' P  ',' O1P',' O2P',' O5*',' C5*',' C4*',' C3*']
    keys = res.keys() ; keys.sort() ;
    k0, kl = keys[0], keys[ len(keys)-1 ]
    for aname in names : knownPositions.append( res[k0][aname] )
    for aname,ai in res[kl].items() :
        if not aname in names : knownPositions.append( res[kl][aname] )

    # make suite-builders
    nas = NAsuiteSampler("%s/rnaSuite" % cbDatapath, 1)
    pts = VecVecFloat(pts)
    irlist = []
    for ri in range(0, len(res)-1) :
        IPs, OPs = [], []
        IPs.append( res[ri][' C5*'] ) ; IPs.append( res[ri][' C4*'] ) ; IPs.append( res[ri][' C3*'] )

        OPs.append( res[ri][' O3*'] )
        OPs.append( res[ri+1][' P  '] )
        OPs.append( res[ri+1][' O5*'] )
        OPs.append( res[ri+1][' C5*'] )
        OPs.append( res[ri+1][' C4*'] )
        OPs.append( res[ri+1][' C3*'] )
        OPs.append( res[ri+1][' O1P'] )
        OPs.append( res[ri+1][' O2P'] )

        OPs.append( res[ri][' O4*'] ) ; OPs.append( res[ri][' C1*'] ) ; OPs.append( res[ri][' C2*'] ) ## sugar
        dna = 1
        if ' O2*' in res[ri].keys() : dna = 0
        if dna==0 : OPs.append( res[ri][' O2*'] ) ## extra oxygen for RNA

        if resns[ri] in [ '  T', '  C', '  U', ' DT', ' DC', ' DU' ] : ## TCU bases pyrimidines
            OPs.append(res[ri][' N1 ']) ; OPs.append(res[ri][' C2 ']) ; OPs.append(res[ri][' O2 ']); OPs.append(res[ri][' N3 '])
            OPs.append(res[ri][' C4 ']) ; OPs.append(res[ri][' C5 ']) ; OPs.append(res[ri][' C6 '])
            if resns[ri][2] == 'T' :
                OPs.append(res[ri][' O4 ']) ; OPs.append(res[ri][' C5M'])
            elif resns[ri][2] == 'C' :
                OPs.append(res[ri][' N4 '])
            elif resns[ri][2] == 'U' :
                OPs.append(res[ri][' O4 '])
            else : assert 0
        elif resns[ri] in [ '  A', '  G', ' DA', ' DG' ] : ## AG bases purines
            OPs.append( res[ri][' N1 ']) ; OPs.append( res[ri][' C2 ']) ; OPs.append( res[ri][' N3 '])
            OPs.append( res[ri][' C4 ']) ; OPs.append( res[ri][' C5 ']) ; OPs.append( res[ri][' C6 '])
            OPs.append( res[ri][' N7 ']) ; OPs.append( res[ri][' C8 ']) ; OPs.append( res[ri][' N9 '])
            if resns[ri][2] == 'A' :
                OPs.append( res[ri][' N6 '] )
            elif resns[ri][2] == 'G' :
                OPs.append( res[ri][' O6 '] ) ; OPs.append( res[ri][' N2 '] )
            else : assert 0
        else : assert 0
        blist.append( NAsuiteBuilder(pts, VecInt(IPs), VecInt(OPs), consts, "NAsuiteBuilder %s" % resids[ri], nas, resns[ri][2], dna, "C3'-endo", "C3'-endo") )
        #irlist.append( SphPosRestr(pts, VecInt([res[ri][' C1*']]), VecFloat( pts[res[ri][' C1*']] ), 5) )
        irlist.append( SphPosRestr(pts, VecInt([res[ri+1][' P  ']]), VecFloat( pts[res[ri+1][' P  ']] ), 0.5) )
        irlist.append( SphPosRestr(pts, VecInt([res[ri+1][' C4*']]), VecFloat( pts[res[ri+1][' C4*']] ), 0.5) )
        irlist.append( SphPosRestr(pts, VecInt([res[ri+1][' C3*']]), VecFloat( pts[res[ri+1][' C3*']] ), 0.5) )
        numtrials.append(500)

    #irlist = []
    #irlist.append( SphPosRestr(pts, VecInt([res[keys[len(keys)-1]][' P  ']]), VecFloat( pts[ res[keys[len(keys)-1]][' C4*'] ] ), 12) )

    overlapReductions = MapIntIntFloat()

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, pts)

    for b in blist : b.describe()

    clashExclInds = makeClashExcl_1st2ndCovNbr(covconn, blist)
    radii = VecFloat(radii)

    strategy(blist, pts, radii, overlapReductions, knownPositions, clashExclInds, numtrials, irlist, modelRenderer)

if __name__ == '__main__' :
    nuMain()
