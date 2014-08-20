from data import resAtoms, nuConn
from pdbr import protein, line2resid, line2resn, line2atomname, line2crd, line2chid, line2resnum, line2resic
from samplers import NAsampler
from builders import NAbuilder, VecInt, VecVecFloat, VecFloat
from data import consts, vdwr
from strategy import makeClashExcl_1st2ndCovNbr, strategy
from geometry import MapIntIntFloat
from restraints import SphPosRestr
import prot2res

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

    # put 1st residue and 2nd residue's P, O5* as known points. last res is also known except its P and O5*
    for k,v in res[0].items() : knownPositions.append(v)
    knownPositions.append(res[1][' P  '])
    knownPositions.append(res[1][' O5*'])
    keys = res.keys() ; keys.sort()
    for k,v in res[ keys[len(keys)-1] ].items() :
        if not k in [' P  ',' O5*'] : knownPositions.append(v)

    # make builders for 1st residue till secondlast.
    nas = NAsampler()
    pts = VecVecFloat(pts)
    for ri in range(1, len(res)-1) :
        IPs, OPs = [], []
        IPs.append( res[ri-1][' O3*'] ) ; IPs.append( res[ri][' P  '] ) ; IPs.append( res[ri][' O5*'] )

        OPs.append( res[ri][' O1P'] ) ; OPs.append( res[ri][' O2P'] ) ## backbone
        OPs.append( res[ri][' C5*'] ) ; OPs.append( res[ri][' C4*'] )
        OPs.append( res[ri][' C3*'] ) ; OPs.append( res[ri][' O3*'] )
        OPs.append( res[ri+1][' P  '] ) ; OPs.append( res[ri+1][' O5*'] )

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
        blist.append( NAbuilder(pts, VecInt(IPs), VecInt(OPs), consts, "NAbuilder %s" % resids[ri], nas, resns[ri][2], dna, "C3'-endo") )
        numtrials.append(10000)

    irlist = []
    irlist.append( SphPosRestr(pts, VecInt([res[keys[1]][' C4*']]), VecFloat( pts[ res[keys[1]][' C4*'] ] ), 1.0) )
    #irlist.append( SphPosRestr(pts, VecInt([res[keys[1]][' C1*']]), VecFloat( pts[ res[keys[1]][' C1*'] ] ), 1.0) )
    irlist.append( SphPosRestr(pts, VecInt([res[keys[len(keys)-1]][' P  ']]), VecFloat( pts[ res[keys[len(keys)-1]][' P  '] ] ), 1.5) )

    overlapReductions = MapIntIntFloat()

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, pts)

    for b in blist : b.describe()

    clashExclInds = makeClashExcl_1st2ndCovNbr(covconn, blist)
    radii = VecFloat(radii)

    strategy(blist, pts, radii, overlapReductions, knownPositions, clashExclInds, numtrials, irlist, modelRenderer)

if __name__ == '__main__' :
    nuMain()
