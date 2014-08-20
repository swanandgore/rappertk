import math, os, sys, string

from builders import Builder, VecInt, VecFloat, VecVecFloat, PeptideBuilder, CBbuilder, ChiBuilder, CNCaBuilder, InformedPeptideBuilder
from geometry import MapIntIntFloat
from samplers import OmegaSampler, PhipsiSampler, BBdepChiSampler, InformedPhipsiSampler
from restraints import DistanceRestraint, SphPosRestr, CentroidPosRestraint
from pdbr import protein, line2resid, line2resn, line2crd, line2atomname, line2resnum, line2resic, line2chid, makeResid
import prot2res

cbDatapath = os.environ["RTKROOT"] + "/data/"

from data import vdwr, resAtoms, consts, mcConn, scConn
from cbutils import findBuilderOrder, findBuilderRestraintOrder, Build


class ChiSamplerProvider :
    def __init__(self) : self.samplers = {}
    def get(self, resn) :
        if not resn in self.samplers.keys() :
            self.samplers[resn] = BBdepChiSampler( "%s/BBdepSCrot/sc%s" % (cbDatapath,resn) )
        return self.samplers[resn]

class PhipsiSamplerProvider :
    def __init__(self) : self.samplers = {}
    def get(self, resn) :
        if not resn in self.samplers.keys() :
            self.samplers[resn] = PhipsiSampler( cbDatapath + "/PhipsiWeightedProp/ps%s" % resn )
        return self.samplers[resn]

## some sidechain rotamer states are so few that they can be sampled serially without weighted sampling, eg till TRP
## from ASN onwards, sessioinng is not reqd, there are enough states to yield a different samples in 20 trials
orderedSCsampling = [ "PRO", "CYS", "SER", "VAL", "THR", "PHE", "TYR", "ASP", "ILE", "HIS", "LEU", "TRP",]
numRot = {} # in Dunbrack's backbone dependent rotamer library
numRot["GLY"] = 0 ; numRot["ALA"] = 0 ; numRot["PRO"] = 2 ;
numRot["CYS"] = 3 ; numRot["SER"] = 3 ; numRot["VAL"] = 3 ; numRot["THR"] = 3 ;
numRot["PHE"] = 6 ; numRot["TYR"] = 6 ; numRot["ASP"] = 9 ;
numRot["ILE"] = 9 ; numRot["HIS"] = 9 ; numRot["LEU"] = 9 ; numRot["TRP"] = 9 ;
numRot["ASN"] = 18 ;
numRot["GLU"] = 27 ; numRot["MET"] = 27 ;
numRot["GLN"] = 36 ;
numRot["ARG"] = 81 ; numRot["LYS"] = 81 ;
useChiSessions = [ "PRO", "CYS", "SER", "VAL", "THR", "PHE", "TYR", "ASP", "ILE",
                "HIS", "LEU", "TRP", "ASN", "GLU", "MET", "GLN", "ARG", "LYS", "ALA", "GLY",]

pepfwd = 1
pepSession = 0

def covConnect(covconn, res, indi, namei, indj, namej) :
    i = res[indi][namei]
    j = res[indj][namej]
    if not i in covconn.keys() : covconn[i] = set()
    if not j in covconn.keys() : covconn[j] = set()
    covconn[i].add(j)
    covconn[j].add(i)

def removeSC(res, pts, forIndices) :
    newpts = []
    for index in forIndices :
        for an in res[index].keys() :
            if not an in [ ' N  ', ' CA ', ' C  ', ' O  ' ] : del res[index][an]
    for index in forIndices :
        for an,pi in res[index].items() :
            newpts.append( pts[pi] )
            res[index][an] = len(newpts)-1
    return res, newpts

if __name__ == "__main__" :
    ## make PhipsiSamplers
    glyPP = PhipsiSampler(cbDatapath + "/PhipsiWeightedProp/psGLY")
    omegaSampler = OmegaSampler()
    numtrials = []

    ## for a chain in pdbfile, rebuild using some CA-sphere restraint radius
    ## assume all else in pdbfile as just steric obstructions
    pdbfile, chid, caRad, scRad = sys.argv[1:5]
    assert len(chid) == 1
    caRad = string.atof(caRad)
    scRad = string.atof(scRad)
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    
    loop_started, index, startindex, endindex = None, 0, None, None
    keys = res.keys() ; keys.sort()
    for index in keys :
        if chids[index] == chid and not loop_started : loop_started, startindex = 1, index
        if chids[index] != chid and loop_started : loop_started, endindex = None, index-1
    if startindex != None and endindex == None : endindex = keys[len(keys)-1]
    assert startindex != None and endindex != None and startindex < endindex
    for index in range(startindex,endindex+1) : assert chid == chids[index]

    ## insert a dummy GLY residue with just atom C at the beginning for bootstrapping
    keys = list(res.keys())
    keys.sort()
    newkey = keys[0] - 1
    pts.append([0,0,0])
    res[newkey] = {' C  ': len(pts)-1}
    resnums[newkey] = '998' ; inscodes[newkey] = ' ' ; resns[newkey] = 'GLY' ; chids[newkey] = chid
    resids[newkey] = makeResid(resns[newkey], chids[newkey], resnums[newkey], inscodes[newkey])
    firstindex = newkey
    ## insert a dummy GLY residue with atoms N,CA at the end of traced chain for last CO and sidechain
    newkey = keys[len(keys)-1] + 1
    pts.append([0,0,0]) ; pts.append([0,0,0])
    res[newkey] = { ' N  ': len(pts)-2, ' CA ': len(pts)-1 }
    resnums[newkey] = '999' ; inscodes[newkey] = ' ' ; resns[newkey] = 'GLY' ; chids[newkey] = chid
    resids[newkey] = makeResid(resns[newkey], chids[newkey], resnums[newkey], inscodes[newkey])
    lastindex = newkey

    mconly = 1

    knownPositions = [] ## make known positions
    unknownIndices = [firstindex, lastindex] + range(startindex, endindex+1)
    for index, val in res.items() :
        for aname, pi in val.items() :
            if index not in unknownIndices : knownPositions.append(pi)

    if mconly == 1 : res, pts = removeSC(res, pts, unknownIndices)

    radii = [0] * len(pts) ## make radii
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else : radii[pi] = vdwr['XXX'][aname[1]]
    radii = VecFloat(radii)
    pts = VecVecFloat(pts)

    ## make CA-tracing restraints.
    ## this can be read in but currently i am just copying from given structure
    irlist = []
    for index in range(startindex, endindex+1) :
        print resids[index]
        ## CA restraint
        irlist.append( SphPosRestr(pts, VecInt([res[index][' CA ']]), VecFloat(pts[res[index][' CA ']]), caRad) )
        ## sidechain centroid restraints
        IPs, mean, count = [], [0.,0.,0.], 0
        for k,v in res[index].items() :
            if k in (' N  ',' CA ',' C  ',' O  ',' CB ') : continue
            IPs.append(v) ; count = count + 1
            mean[0] += pts[v][0] ; mean[1] += pts[v][1] ; mean[2] += pts[v][2]
        if count == 0 : continue
        mean[0] /= count ; mean[1] /= count ; mean[2] /= count
        irlist.append( CentroidPosRestraint(pts, VecInt(IPs), VecFloat(mean), scRad) )

    print "Residues and atom-numbers-----------------------"
    keys = res.keys() ; keys.sort()
    for k in keys : print resids[k], res[k]
    print "------------------------------------------------\n\n\n"

    blist = []
    ## make bootstrap-builder at beginning and CO builder in the end
    print "Making Backbone Builders"
    IPs, OPs = [], [ res[firstindex][' C  '], res[startindex][' N  '], res[startindex][' CA '] ]
    blist.append( CNCaBuilder( pts, VecInt(IPs), VecInt(OPs), consts, "CNCaBuilder", VecFloat(pts[res[startindex][' CA ']]), caRad ) )
    numtrials.append(5)
    psp = PhipsiSamplerProvider() ## make peptide builders
    isamps = []
    for index in range(startindex+1, endindex+1) : # peptide builders. i'th builder builds i'th CA
        IPs, OPs = [], []
        IPs.append(res[index-2][' C  ']) ; IPs.append(res[index-1][' N  ']) ; IPs.append(res[index-1][' CA '])
        OPs.append(res[index-1][' C  ']) ; OPs.append(res[index-1][' O  ']) ; OPs.append(res[index][' N  ']) ; OPs.append(res[index][' CA '])
        isamps.append( InformedPhipsiSampler(VecFloat(pts[res[index][' CA ']]), caRad, resns[index-1]) )
        blist.append( InformedPeptideBuilder(pts, VecInt(IPs), VecInt(OPs), consts, "InformedPeptideBuilder %s" % (resids[index]), isamps[len(isamps)-1], pepfwd, pepSession) )
        numtrials.append( int(caRad*5) )
    IPs = [ res[endindex-1][' C  '], res[endindex][' N  '], res[endindex][' CA '] ]
    OPs = [ res[endindex][' C  '], res[endindex][' O  '], res[lastindex][' N  '], res[lastindex][' CA '] ]
    blist.append( PeptideBuilder(pts, VecInt(IPs), VecInt(OPs), consts, "Final Peptide Builder", psp.get(resns[index-1]), omegaSampler, pepfwd, pepSession) )
    numtrials.append(20)
    print "Making Sidechain Builders"
    csp = ChiSamplerProvider() ## make CB / sidechain builders
    for index in range(startindex, endindex+1) : # CB and sidechain builders
        if mconly == 1 : continue
        resn = resns[index] ; resid = resids[index]
        if resn == "GLY" : continue
        IPs, OPs = [], [] # CB
        IPs.append(res[index][' N  ']) ; IPs.append(res[index][' C  ']) ; IPs.append(res[index][' CA '])
        OPs.append(res[index][' CB '])
        blist.append( CBbuilder(pts, VecInt(IPs), VecInt(OPs), consts, "CBbuilder %s" % (resid)) )
        numtrials.append(1)
        if resn == "ALA" : continue
        IPs, OPs = [], [] # sidechains
        IPs.append(res[index-1][' C  '])
        IPs.append(res[index][' N  ']) ; IPs.append(res[index][' CA ']) ; IPs.append(res[index][' C  '])
        if index == endindex : IPs.append(res[lastindex][' N  '])
        else : IPs.append(res[index+1][' N  '])
        IPs.append(res[index][' CB '])
        for name in resAtoms[resn]:
            if name in [ ' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ] : continue
            OPs.append(res[index][name])
        orsamp, useSession = 0, 0
        if resn in orderedSCsampling : orsamp = 1 
        if resn in useChiSessions : useSession = 1
        blist.append( ChiBuilder(pts, VecInt(IPs), VecInt(OPs), consts, "ChiBuilder %s" % (resid), resn, csp.get(resn), useSession, orsamp) )
        if resn in orderedSCsampling : numtrials.append( numRot[resn] + 1 )
        else : numtrials.append(20)

    print "Setting up covalent connections"
    covconn = {} ## generate covalent connectivity for residues start-1 till end+1, both included
    for index in range(startindex, endindex+1) :
        resn = resns[index]
        conlist = mcConn
        if mconly != 1 and (len(scConn[resn]) > 0 or resn == "ALA") : conlist = mcConn + scConn[resn] + [(' CA ',' CB ')]
        for c in conlist : covConnect(covconn, res, index, c[0], index, c[1])

    for index in range(startindex+1, endindex+1) : ## setup C-N conn between residues
        covConnect(covconn, res, index, ' N  ', index-1, ' C  ')
    covConnect(covconn, res, startindex, ' N  ', firstindex, ' C  ')
    covConnect(covconn, res, lastindex, ' N  ', endindex, ' C  ')
    covConnect(covconn, res, lastindex, ' N  ', lastindex, ' CA ')


    from strategy import makeClashExcl_1st2ndCovNbr, strategy
    clashExclInds = makeClashExcl_1st2ndCovNbr(covconn, blist)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, pts)

    overlapReductions = MapIntIntFloat()
    assert len(radii) == len(pts)
    strategy(blist, pts, radii, overlapReductions, knownPositions, clashExclInds, numtrials, irlist, modelRenderer)
