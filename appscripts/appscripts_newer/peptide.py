import math, os, sys

from geometry import MapIntIntFloat
from builders import Builder, VecInt, VecFloat, VecVecFloat, PeptideBuilder, CBbuilder, ChiBuilder
from samplers import OmegaSampler, PhipsiSampler, BBdepChiSampler
from restraints import DistanceRestraint

cbDatapath = os.environ["RTKROOT"] + "/data/"

from data import vdwr, resAtoms, consts, mcConn, scConn

##{{{ ModelRenderer class

def makePDBatomline(atomname, atomnum, resname, resnum, chid, inscode, x, y, z) :
    if atomnum > 9999 : atomnum = 9999
    if atomnum < -999 : atomnum = -999
    assert len(atomname) == 4
    assert len(resname) == 3
    assert len(chid) == 1
    assert len(inscode) == 1
    return "ATOM   %4d %4s %3s %1s%4s%1s   %8.3f%8.3f%8.3f" % (atomnum, atomname, resname, chid, resnum, inscode, x, y, z)          

class ModelRenderer :
    def __init__(self, res, resnames, chids, resnums, inscodes, dummies, fn) :
        '''dummies is a set of dummy indices which need not be rendered.'''
        self.modnum = 1
        self.res, self.resnames, self.chids, self.resnums, self.inscodes = res, resnames, chids, resnums, inscodes
        self.dummies = dummies # residue inds which need not be rendered
        try : self.f = open(fn, 'w')
        except :
            print "File [%s] cannot be opened for writing, fatal error" % fn
            sys.exit(1)

    def render(self, pts, ptinds=None) :
        print >> self.f, "MODEL     %4d" % self.modnum
        self.modnum = self.modnum + 1
        keys = self.res.keys() ; keys.sort()
        for i in keys :
            if i in self.dummies : continue
            for name,pti in self.res[i].items() :
                if ptinds and not pti in ptinds : continue
                #print self.resnames[i], self.resnums[i], self.inscodes[i], pts[pti][0], pts[pti][1], pts[pti][2], self.chids[i]
                print >> self.f, makePDBatomline(name, pti, self.resnames[i], self.resnums[i], self.chids[i], self.inscodes[i], pts[pti][0], pts[pti][1], pts[pti][2])
        print >> self.f, "ENDMDL"

    def __del__(self) :
        self.f.close()

##}}}

def init_C_N_CA(pts, a, b, c) :
    C_N = consts.get("C_N")
    N_CA = consts.get("N_CA")
    C_N_CA = consts.get("C_N_CA")
    pts[a] = VecFloat([0,0,0]) # C
    pts[b] = VecFloat([C_N,0,0]) # N
    pts[c] = VecFloat([C_N - N_CA*math.cos(C_N_CA), N_CA*math.sin(C_N_CA), 0]) # CA

if __name__ == "__main__" :
    print 'here 1'
    ## make PhipsiSamplers
    glyPP = PhipsiSampler(cbDatapath + "/PhipsiWeightedProp/psGLY")
    print 'here 2'
    omegaSampler = OmegaSampler()
    numtrials = []
    radii = VecFloat()
    ## residues ##################################################
    res = [ { ' C  ':0, } ] # 0th residue
    pts, knownPositions = VecVecFloat(), []
    pts.push_back( VecFloat([-999,-999,-999]) )
    radii.push_back( vdwr['GLY'][' C  '] )
    resseq, ac = ['VAL',]*20, 1 #seq to be built
    for i in range( len(resseq) ) :
        ares, resname = {}, resseq[i]
        for name in resAtoms[resname] : #[ ' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ] :
            ares[name] = ac
            ac = ac + 1
            pts.push_back( VecFloat([-999,-999,-999]) )
            radii.push_back( vdwr[resname][name] )
        res.append(ares)
    ares = { ' N  ':ac, ' CA ':ac+1} # last residue
    pts.push_back( VecFloat([-999,-999,-999]) )
    pts.push_back( VecFloat([-999,-999,-999]) )
    radii.push_back( vdwr['GLY'][' N  '] )
    radii.push_back( vdwr['GLY'][' CA '] )
    res.append(ares)

    print 'RESIDUES------------------------'
    for ares in res : print ares
    print '================================\n\n\n'

    ## bootstrapping
    init_C_N_CA( pts, res[0][' C  '], res[1][' N  '], res[1][' CA '] )
    knownPositions.append(res[0][' C  '])
    knownPositions.append(res[1][' N  '])
    knownPositions.append(res[1][' CA '])

# {{{ prepare builders, covalent-connectivity for mainchain
    covconn, blist = [], []
    for i in range(pts.size()) : covconn.append([])
    for i in range( 1, len(res)-1 ) : # PeptideBuilders
        IPs, OPs = [], []
        IPs.append( res[i-1][' C  '] )
        IPs.append( res[i][' N  '] )
        IPs.append( res[i][' CA '] )
        OPs.append( res[i][' C  '] )
        OPs.append( res[i][' O  '] )
        OPs.append( res[i+1][' N  '] )
        OPs.append( res[i+1][' CA '] )
        blist.append( PeptideBuilder(pts, VecInt(IPs), VecInt(OPs), consts, "PeptideBuilder for CA %d" % i, glyPP, omegaSampler) )
        numtrials.append(50)
    for i in range( 1, len(res)-1 ) : # covconn for mainchains of residues with full mainchain
        covconn[ res[i][' N  '] ].append( res[i-1][' C  '] ) ## prev c-n
        covconn[ res[i-1][' C  '] ].append( res[i][' N  '] )
        for conn in mcConn :
            covconn[ res[i][conn[0]] ].append( res[i][conn[1]] )
            covconn[ res[i][conn[1]] ].append( res[i][conn[0]] )
    n = len(res)-1
    covconn[ res[n][' N  '] ].append( res[n][' CA '] ) # covconn for res-last N,CA
    covconn[ res[n][' CA '] ].append( res[n][' N  '] )
    covconn[ res[n-1][' C  '] ].append( res[n][' N  '] )
    covconn[ res[n][' N  '] ].append( res[n-1][' C  '] )
    print "Mainchain builders and connectivity done\n\n\n"
# }}}

# {{{ prepare builders, covconn for CB and sidechain
    chisamplers = {}
    chisamplers["VAL"] = BBdepChiSampler( "%s/BBdepSCrot/sc%s" % (cbDatapath,"VAL") )
    chisamplers["ARG"] = BBdepChiSampler( "%s/BBdepSCrot/sc%s" % (cbDatapath,"ARG") )
    for i in range( 1, len(res)-1 ) :
        resname = resseq[i-1]
        if resname == 'GLY' : continue
        IPs, OPs = [], [] # CBbuilders and covconn
        IPs.append(res[i][' N  '])
        IPs.append(res[i][' C  '])
        IPs.append(res[i][' CA '])
        OPs.append(res[i][' CB '])
        blist.append( CBbuilder(pts, VecInt(IPs), VecInt(OPs), consts, "CBbuilder for %d" % i) )
        numtrials.append(1)
        covconn[ res[i][' CA '] ].append( res[i][' CB '] ) ## ca-cb
        covconn[ res[i][' CB '] ].append( res[i][' CA '] )
        if resname == 'ALA' : continue
        IPs, OPs = [], [] # sidechain-builders and covconn
        IPs.append(res[i-1][' C  '])
        IPs.append(res[i][' N  '])
        IPs.append(res[i][' CA '])
        IPs.append(res[i][' C  '])
        IPs.append(res[i+1][' N  '])
        IPs.append(res[i][' CB '])
        for name in resAtoms[resname]:
            if name in [ ' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ] : continue
            OPs.append(res[i][name])
        blist.append( ChiBuilder(pts, VecInt(IPs), VecInt(OPs), consts, "ChiBuilder for %d %s" % (i,resname), resname, chisamplers[resname]) )
        numtrials.append(10)
        for conn in scConn[resname] :
            covconn[ res[i][conn[0]] ].append( res[i][conn[1]] )
            covconn[ res[i][conn[1]] ].append( res[i][conn[0]] )
    print "Sidechain builders and connectivity done\n\n\n"
# }}}
    for b in blist : b.describe()
    print '\n\n\n'

    print "COVCONN----------------------"
    for ci in range(len(covconn)) : print "%d----" % ci, covconn[ci]
    print '-----------------------------\n\n\n'

    from strategy import makeClashExcl_1st2ndCovNbr
    clashExclInds = makeClashExcl_1st2ndCovNbr(covconn, blist)

# {{{ prepare restraints
    irlist = []
    #irlist.append( DistanceRestraint(pts, VecInt([res[1][' CA '],res[3][' CA ']]), 4, 6) )
# }}}

## app-specific code ends here and strategy begins
    resseq = ['GLY'] + resseq + ['GLY']
    modelRenderer = ModelRenderer(res, resseq, pts)
    from strategy import strategy
    overlapReductions = MapIntIntFloat()
    strategy(blist, pts, radii, overlapReductions, knownPositions, clashExclInds, numtrials, irlist, modelRenderer)
