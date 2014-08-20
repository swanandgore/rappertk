from pdbr import protein
import prot2res
import re, string, sys
from geometry import calcDist, CAtraceGH, MapIntIntFloat
from builders import VecInt, VecFloat, VecBuilder, VecVecFloat, Rotator, TransRotator, BuilderGroup
from restraints import DistanceRestraint
from data import vdwr, consts, PROBE_DISULFIDE_OVERLAP_MARGIN
from prepareChain import  incompleteSCcorrection
import data
from data import vdwr, consts

def findLigDists(ligname, res, pts, resns) :
    ligdists = {}
    ligind = None
    for ind in res.keys() :
        if resns[ind] == ligname : ligind = ind ; break
    assert ligind
    for ind in res.keys() :
        if resns[ind] == ligname or 'HOH' in resns[ind]:
            continue
        for name,ai in res[ind].items() :
            mindist = 1e10
            for lan,lai in res[ligind].items() :
                ldist = calcDist( VecFloat(pts[lai]), VecFloat(pts[ai]) )
                #print "ldist", lan, ldist, pts[lai], pts[ai]
                if ldist < mindist : mindist = ldist
            if mindist < 1e9 : ligdists[ind] = mindist
    return ligdists

def parseAtomnames(l) :
    '''from line containing tokens of type [aaaa], extract 4-character atomnames and return the list'''
    i = 0 ; anames = []
    while i <= len(l)-4 :
        if l[i] == '[' :
            assert l[i+5] == ']'
            anames.append( l[i+1:i+5] )
            i = i + 6
        else : i = i + 1
    return anames

def makeLigBuilder(res, resns, pts, ligfile) :
    '''make a series of builders and group into a GroupBuilder.
        Currently recongnize only the lines starting with mindist, rotbond and init'''
    ligind, restraints, builders = None, [], []
    tol1, tol2, init1, init2 = None, None, None, None
    allAnames, ligname = set(), None
    for l in open(ligfile, 'r').readlines() :
        l = re.sub("\n", "", l)
        l = re.sub(" *#.*", "", l)
        if len(l) == 0 : continue
        flds = l.split()
        anames = parseAtomnames(l)
        for an in anames : allAnames.add(an)
        if l[0:7] == 'mindist' :
            for fi in range(1,len(flds)) :
                if flds[fi][0] == '[' : break
                tol = string.atof(flds[fi])
                restraints.append( DistanceRestraint(VecInt([res[ligind][anames[0]], res[ligind][anames[1]]]), "DistanceRestraint within ligand", tol, 1e10) )
        elif l[0:7] == 'ligname' :
            for k,v in resns.items() :
                if v == flds[1] : ligname = flds[1] ; ligind = k ; break
            assert ligind != None
        elif l[0:7] == 'rotbond' :
            min, max, step = string.atof(flds[1]), string.atof(flds[2]), string.atof(flds[3]),
            IPs, OPs = [ res[ligind][anames[0]],res[ligind][anames[1]] ], []
            for an in anames[2:] : OPs.append(res[ligind][an])
            builders.append( Rotator(VecInt(IPs), VecInt(OPs), None, "Rotator", min, max, step) )
        elif l[0:4] == 'init' :
            tol1, tol2 = string.atof(flds[1]), string.atof(flds[2])
            init1, init2 = anames[0], anames[1]
        else : print "Unrecognized file format for %s, exiting....." % ligfile ; assert None
    assert tol1
    IPs, OPs = [], [res[ligind][init1], res[ligind][init2]]
    allAnames.remove(init1) ; allAnames.remove(init2)
    for an in allAnames : OPs.append(res[ligind][an])
    builders = [ TransRotator(VecInt(IPs),VecInt(OPs),None,"TransRotator", VecFloat(pts[res[ligind][init1]]),tol1, VecFloat(pts[res[ligind][init2]]),tol2) ] + builders
    return builders, restraints, ligname

def main(pdbfile, ligfile, closeCutoff, caRad, scRad, scReduction, guidedSampling, outpdb, popsize, backtrack, nmodels, restrGen, mapfile,ranker) :

    prot = protein(pdbfile, read_hydrogens=0, read_waters=1, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    lignames = [] ; ligkey = []
    from data import three2one
    for k,v in resns.items():
        if v  in three2one.keys() or v == 'HOH':
            continue
        else :
            print v
            if v not in lignames:
                if len(res[k].values()) > 2 : 
                    if v in ["SO4","PO4","SF4"]:
                        continue
                    lignames.append(v)
                    ligkey.append(k)
    if len(lignames) == 0:
        print "# WARN: No ligand found"
        import sys; sys.exit()

    if len(lignames) > 1:
        print "#WARN: More than one ligand found"

    for x in ligkey :
        print "ligname ",resns[x]
        for a,b in res[x].items():
            print "[%s]"%a,
        print
    ## find ligand name
        











def callmain() :
    from commonOptions import makeParser, parseOptions, addXrayOptions
    parser = makeParser()
    parser = addXrayOptions(parser)
    parser.add_option("--ligfile", action='store', type='string', dest='ligfile', help='ligand description, see msq.ligdesc for description of MSQ in 1di9.pdb')
    parser.add_option("--around-ligand", action='store', type='float', dest='closeCutoff', help='min-dist between residue/ligand to be considered close to ligand', default=10)
    parser.add_option("--mapfile", action='store', type='string', dest='mapfile', help='ccp4 map restraining the shape of the molecule', default=None)
    parser.remove_option("--mconly")
    options = parseOptions(parser)

    import misc
    misc.setVerbosity(options.verbose)

    xrayRestGen = None
    from prepareChain import XrayRestraintsGenerator
    from xcheck import XrayScorer, XrayRanker
    ranker = None
    if options.mtzfn != None : 
        xrayRestGen = XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, options.sigXmin, options.sigXmax, options.sigXmean)



        esmin, esmax, esmean, rcmult = .000, 5., .0, 5 
        ranker = XrayRanker("phased.mtzFPFCPHIC_A.map", "map", "FC", "PHIC", "2F1-F2", esmin, esmax)
        ranker.rankChildren = rcmult ; ranker.rankRecurse = 1
        ranker.rankLeaderBuilderOnly = None ; ranker.rankGivenEnsemble = None
        
    main(options.pdbfile, options.ligfile, options.closeCutoff, options.caRad, options.scRad, options.scReduction, options.guidedSampling,options.outpdb, options.popsize, options.backtrack, options.nmodels,xrayRestGen, options.mapfile,ranker)
         

if __name__ == "__main__" : callmain()
