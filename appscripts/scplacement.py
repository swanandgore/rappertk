import math, os, sys, string

from builders import Builder, VecInt, VecFloat, VecVecFloat, PeptideBuilder, CBbuilder, ChiBuilder, CNCaBuilder, NanchorBuilder
from builders import VecBuilder, BuilderGroup
from geometry import MapIntIntFloat, CAtraceGH
from samplers import OmegaSampler, PhipsiSampler, BBdepChiSampler, PRLsampler
from restraints import DistanceRestraint, SphPosRestr, CentroidPosRestraint
from pdbr import protein, line2resid, line2resn, line2crd, line2atomname, line2resnum, line2resic, line2chid, makeResid
import prot2res

cbDatapath = os.environ["RTKROOT"] + "/data/"

from data import vdwr, resAtoms, consts, mcConn, scConn
from cbutils import findBuilderOrder, findBuilderRestraintOrder, Build


class PRL_BB_SC_SamplerProvider :
    def __init__(self) :
        self.samplers = {}
        self.aaMap = {
                'PHE':'F', 'TYR':'Y', 'TRP':'W',
                'GLY':'G', 'ALA':'A', 'VAL':'V', 'LEU':'L', 'ILE':'I', 'PRO':'P',
                'ASP':'D', 'GLU':'E', 'HIS':'H', 'LYS':'K', 'ARG':'R',
                'CYS':'C', 'MET':'M', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q'
            }
    def get(self, resn) :
        if not resn in self.samplers.keys() :
            self.samplers[resn] = PRLsampler( "%s/richardson.lib" % cbDatapath, self.aaMap[resn] )
        return self.samplers[resn]

class Dunbrack_BB_SC_SamplerProvider :
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

class OmegaSamplerProvider :
    def __init__(self) :
        self.omegaSampler = OmegaSampler(0)
        self.preproOmegaSampler = OmegaSampler(1)
    def get(self,nextresn) :
        if nextresn == 'PRO' : return self.preproOmegaSampler
        return self.omegaSampler

def makeCBbuilder(index, res, resns, resids) :
    IPs, OPs = [ res[index][' N  '],res[index][' C  '],res[index][' CA '] ], [ res[index][' CB '] ]
    return CBbuilder(VecInt(IPs), VecInt(OPs), consts, "CBbuilder %s" % (resids[index]), resns[index])

def makeChiBuilder(i, j, k, res, resns, resids) :
    IPs, OPs = [ res[i][' C  '],res[j][' N  '],res[j][' CA '],res[j][' C  '],res[k][' N  '],res[j][' CB '] ], []
    for name in resAtoms[resns[j]]:
        if name in [ ' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ] : continue
        OPs.append(res[j][name])
    return ChiBuilder(VecInt(IPs), VecInt(OPs), consts, "ChiBuilder %s" % (resids[j]), resns[j], csp.get(resns[j]), 0, 0)

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
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='template PDB file')
    parser.add_option("--chid", action='store', type='string', dest='chid', help='chain identifier', default=' ')
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid', default=2)
    parser.add_option("--num-models-wanted", action='store', type='int', dest='nmodels', help='number of models desired. number of attempts is generally 10 times this', default=100)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--outpdb", action='store', type='str', dest='outpdb', help='file to write the models to')
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=1)
    (options, args) = parser.parse_args()

    import misc
    misc.setVerbosity(options.verbose)
    #print options ; print args ; sys.exit(0)

    numtrials = []

    ## for a chain in pdbfile, rebuild using some CA-sphere restraint radius
    ## assume all else in pdbfile as just steric obstructions
    assert len(options.chid) == 1
    prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    
    loop_started, index, startindex, endindex = None, 0, None, None
    keys = res.keys() ; keys.sort()
    for index in keys :
        if chids[index] == options.chid and not loop_started : loop_started, startindex = 1, index
        if chids[index] != options.chid and loop_started : loop_started, endindex = None, index-1
    if startindex != None and endindex == None : endindex = keys[len(keys)-1]
    assert startindex != None and endindex != None and startindex < endindex
    for index in range(startindex,endindex+1) : assert options.chid == chids[index]

    keys = list(res.keys())
    keys.sort()
    ## insert a dummy GLY residue with atoms C,O,CA at the beginning of traced chain for last CO and sidechain
    newkey = keys[0]-1
    pts.append([0,0,0]) ; pts.append([0,0,0]) ; pts.append([0,0,0])
    res[newkey] = { ' C  ': len(pts)-1, ' O  ': len(pts)-2, ' CA ': len(pts)-3 }
    resnums[newkey] = '998' ; inscodes[newkey] = ' ' ; resns[newkey] = 'GLY' ; chids[newkey] = options.chid
    resids[newkey] = makeResid(resns[newkey], chids[newkey], resnums[newkey], inscodes[newkey])
    firstindex = newkey
    ## insert a dummy GLY residue with atoms N,CA at the end of traced chain for last CO and sidechain
    newkey = keys[len(keys)-1] + 1
    pts.append([0,0,0]) ; pts.append([0,0,0])
    res[newkey] = { ' N  ': len(pts)-2, ' CA ': len(pts)-1 }
    resnums[newkey] = '999' ; inscodes[newkey] = ' ' ; resns[newkey] = 'GLY' ; chids[newkey] = options.chid
    resids[newkey] = makeResid(resns[newkey], chids[newkey], resnums[newkey], inscodes[newkey])
    lastindex = newkey

    knownPositions = [] ## make known positions. all mainchain (except dummies) are known
    for index in range(startindex, endindex+1) :
        for aname, pi in res[index].items() :
            if aname in [' N  ',' CA ',' C  ',' O  '] :
                if index == startindex and aname == ' N  ' : pass
                elif index == endindex and aname == ' C  ' : pass
                elif index == endindex and aname == ' O  ' : pass
                else : knownPositions.append(pi)

    radii = [0] * len(pts) ## make radii
    for index, val in res.items() :
        for aname, pi in val.items() :
            if index in [firstindex,lastindex] : radii[pi] = 0
            elif resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else : radii[pi] = vdwr['XXX'][aname[1]]
    pts = VecVecFloat(pts)

    irlist = []
    for index in range(startindex, endindex+1) :
        if resns[index] == 'GLY' : continue
        IPs, mean, count = [], [0.,0.,0.], 0
        for k,v in res[index].items() :
            if k in (' N  ',' CA ',' C  ',' O  ',) : continue
            IPs.append(v) ; count = count + 1
            mean[0] += pts[v][0] ; mean[1] += pts[v][1] ; mean[2] += pts[v][2]
        assert len(IPs) > 0
        mean[0] /= count ; mean[1] /= count ; mean[2] /= count
        irlist.append( CentroidPosRestraint(VecInt(IPs), VecFloat(mean), options.scRad) )


    print "Residues and atom-numbers-----------------------"
    keys = res.keys() ; keys.sort()
    for k in keys : print resids[k], res[k]
    print "------------------------------------------------\n\n\n"

    blist = []
    pepTrials = 1000; scTrials = 50; cbTrials = 1 #XXX
    extraRestraintRad = .5
    psp = PhipsiSamplerProvider() ## make peptide builders
    omsp = OmegaSamplerProvider() ## make peptide builders
    #csp = Dunbrack_BB_SC_SamplerProvider() ## make CB / sidechain builders
    csp = PRL_BB_SC_SamplerProvider() ## make CB / sidechain builders

    IPs, OPs = [ res[startindex+1][' N  '], res[startindex][' C  '], res[startindex][' CA '] ], [ res[startindex][' N  '], res[firstindex][' C  '], res[firstindex][' O  '], res[firstindex][' CA '] ]
    blist.append( PeptideBuilder( VecInt(IPs), VecInt(OPs), consts, "RevPepBuilder", psp.get(resns[startindex]), omsp.get("GLY"), 0, 0 ) )
    numtrials.append(pepTrials)
    irlist.append( SphPosRestr(VecInt([res[startindex][' N  ']]), VecFloat(pts[res[startindex][' N  ']]), extraRestraintRad) )

    # for index, make CBbuilder and ChiBuilder that builds index'th CB and sidechain
    indices = range(startindex, endindex+1)
    for index in indices :
        i, j, k = index-1, index, index+1
        if index == startindex : i = firstindex
        if index == endindex : k = lastindex
        if resns[index] != 'GLY' :
            blist.append( makeCBbuilder(j, res, resns, resids) )
            numtrials.append(cbTrials)
            if resns[index] != 'ALA' :
                blist.append( makeChiBuilder(i,j,k, res, resns, resids) )
                numtrials.append(scTrials)

    IPs, OPs = [ res[endindex-1][' C  '], res[endindex][' N  '], res[endindex][' CA '] ], [ res[endindex][' C  '], res[endindex][' O  '], res[lastindex][' N  '], res[lastindex][' CA '] ]
    blist.append( PeptideBuilder( VecInt(IPs), VecInt(OPs), consts, "LastPepBuilder", psp.get(resns[endindex]), omsp.get("GLY"), 1, 0 ) )
    numtrials.append(pepTrials)
    irlist.append( SphPosRestr(VecInt([res[endindex][' C  ']]), VecFloat(pts[res[endindex][' C  ']]), extraRestraintRad) )
    irlist.append( SphPosRestr(VecInt([res[endindex][' O  ']]), VecFloat(pts[res[endindex][' O  ']]), extraRestraintRad) )


    for b in blist : b.describe() ; print

    from strategy import makeClashExcl_1st2ndCovNbr, strategy
    covconn = {} ;
    clashExclInds = makeClashExcl_1st2ndCovNbr(covconn, blist)

    ai2pos, aiSC = [-999] * len(pts), []
    for index,val in res.items() :
        for name, ai in val.items() :
            ai2pos[ai] = index
            if not name in [' N  ',' CA ',' C  ',' O  '] : aiSC.append(ai)
    scReduction = options.scReduction
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), options.scReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, options.outpdb)

    overlapReductions = MapIntIntFloat()
    assert len(radii) == len(pts)

    from BacktrackPopstrategy import BacktrackPopstrategy
    strategy = BacktrackPopstrategy(options.popsize, pts, radii, knownPositions, clashExclInds, overlapReductions, gridHelper,
        blist, numtrials, irlist, modelRenderer, options.nmodels*10, options.nmodels)
    #from PopulationStrategy import PopulationStrategy
    #PopulationStrategy(options.popsize, pts, radii, knownPositions, clashExclInds, overlapReductions, blist, numtrials, irlist, modelRenderer, options.nmodels*10, nmodels).execute()
    #from BasicStrategy import BasicStrategy
    #BasicStrategy(pts, radii, knownPositions, clashExclInds, overlapReductions, blist, numtrials, irlist, modelRenderer, 1000, 10).execute()
    #strategy(blist, pts, radii, overlapReductions, knownPositions, clashExclInds, numtrials, irlist, modelRenderer)
    strategy.execute()
