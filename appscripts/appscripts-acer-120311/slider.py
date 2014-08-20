import os, re
import misc
from pref import randomize
from data import sgtable
from pdbr import protein
from xray import sfall
from xcheck import XrayRanker
from peptidebuild import ModelRenderer
import prepareChain
from pref import cnsRefinement
from loopbuild import Multiloop
import prot2res
from multProtref import restoreChainids

def makeAllGLYmodel(pdbin, pdbout) :
    prot = protein(pdbin, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    newpts = []
    for ri in range(len(res)) :
        resns[ri] = "GLY"
        for an in [" N  "," CA "," C  "," O  "] :
            newpts.append( pts[ res[ri][an] ] )
            res[ri][an] = len(newpts)-1
        for an in res[ri].keys() :
            if an in [" N  "," CA "," C  "," O  "] : continue
            del res[ri][an]
    pts = newpts
    ModelRenderer(res, resns, chids, resnums, inscodes, [], pdbout).render(pts)

def main():
    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='to create all the files during refinement. it shdnt be already present.')
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='starting pdb containing a mainchain-only fragment')

    parser.add_option("--mtz", action='store', type='string', dest='mtz', help='str factors as mtz file')
    parser.add_option("--a", action='store', type='float', dest='a', help='cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='cell spacegroup, in CCP4 notation')
    parser.add_option("--resolution", action='store', type='float', dest='resolution', help='resolution of the data')

    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid', default=2)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=1)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='seed for randomizing', default=None)

    parser.add_option("--sequence", action='store', type='string', dest='sequence', help='seq to slide over the given fragment', default=None)
    parser.add_option("--slideres", action='store', type='string', dest='slideres', help='file containing resids over which sequence is to be slided', default=None)

    (options, args) = parser.parse_args()

    if not os.path.isdir(options.scratchdir) : os.mkdir(options.scratchdir)
    os.chdir(options.scratchdir)

    misc.setVerbosity(options.verbose)

    randomize(options.randomize)

    slideres = []
    for l in open(options.slideres, 'r').readlines() :
        l = re.sub("\n", "", re.sub("#.*", "", l))
        if len(l) == 0 : continue
        if not l[0] == "#" : slideres.append( l[1:len(l)-1] )
    print "SLIDERES", slideres
    options.slideres = slideres

    refprot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    refres, refresids, refresnums, refresns, refchids, refinscodes, refpts = prot2res.readProtRes(refprot)
    for ri in range(len(refres)) :
        assert refresns[ri] == "GLY"
        for an in refres[ri].keys() : assert an in [" N  "," CA "," C  "," O  "]

    from restraints import EDrestraint ; EDrestraint.setPenalty(0.) ; EDrestraint.setScatRad(1.) ;
    esmin, esmax, esmean, rcmult, xscoreCutoff = .000, 5., .0, 5, 0.9

    ## make a SA-weighted 2FoFc map based on given given mconly fragment
    cnsArgs = {}
    for cycle in range(20) : cnsArgs[cycle] = {} ; cnsArgs[cycle]["num_cycles"] = 1 ; cnsArgs[cycle]["temperature"] = 50 ; cnsArgs[cycle]["harmCA"] = 1
    cnsArgs[1]["num_cycles"] = 2 ; cnsArgs[1]["temperature"] = 3000
    sfall(options.pdbfile, options.mtz, "rfree.mtz")
    cnsRefinement("rfree.mtz", options.pdbfile, "phased0.mtz", "dontcare.pdb",
                options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg, options.resolution,
                cnsArgs, 0)
    restoreChainids("dontcare.pdb", options.pdbfile)

    #refprot = protein("dontcare.pdb", read_hydrogens=0, read_waters=0, read_hets=0)
    #refres, refresids, refresnums, refresns, refchids, refinscodes, refpts = prot2res.readProtRes(refprot)
        #res, resids, resnums, resns, chids, inscodes, pts = makeFrag( options.sequence[si:si+len(refres)] )

    for si in range( len(options.sequence)-len(options.slideres)+1 ) :
        prot = protein("dontcare.pdb", read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        mutmap = {}
        for resid in options.slideres :
            for k,v in resids.items() :
                if v == resid : mutmap[k] = options.sequence[ si + options.slideres.index(resid) ]
        from data import resAtoms
        from pdbr import makeResid
        from protinfo import AA13
        for ri,aa1 in mutmap.items() :
            print "mutate %5d [%s] -> %s" % (k, resids[ri], aa1)
            for an in resAtoms[ AA13[aa1] ] :
                if an in [" N  "," CA "," C  "," O  "] : continue
                assert not an in res[ri].keys()
                res[ri][an] = len(pts)
                pts.append( [-999., -999., -999.] )
            resns[ri] = AA13[aa1]
            resids[ri] = makeResid( resns[ri], chids[ri], resnums[ri], inscodes[ri] )
        badresids = []
        for ri,aa1 in mutmap.items() : badresids.append( resids[ri] )
        print "badresids", badresids

        ModelRenderer(res, resns, chids, resnums, inscodes, [], "start.pdb").render(pts)

        multiPrepC = prepareChain.PrepareChain("SCL1.0")
        xrayRestGen = [ prepareChain.XrayRestraintsGenerator("phased0.mtz2fofc.map", "map", "FC", "PHIC", "2F1-F2", esmin, esmax, esmean, ["SCLbuilder","ChiBuilder","CBbuilder"], ) ]
        ml = Multiloop("start.pdb", badresids, None, options.caRad, options.scRad, options.scReduction, None, options.popsize,
            options.backtrack, 1, "end.pdb", xrayRestGen, multiPrepC)
        ml.ranker = XrayRanker("phased0.mtz2fofc.map", "map", "FC", "PHIC", "2F1-F2", esmin, esmax)
        ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
        ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble = None
        ml.cellsym = [ options.a, options.b, options.c, options.alpha, options.beta, options.gamma, sgtable[options.sg][0] ]
        ml.run()
        sfall("end.pdb", options.mtz, "rfree.mtz")
        cnsRefinement("rfree.mtz", "end.pdb", "phased%d.mtz"%si, "end%d.pdb"%si,
                options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg, options.resolution,
                cnsArgs, 1)

## will make res, pts, resns etc for given sequence. chain id is ' ', residue numbers are 1, 2, .....
def makeFrag(seq) :
    from data import resAtoms
    from pdbr import makeResid
    from protinfo import AA13
    res, resids, resnums, resns, chids, inscodes, pts = {}, {}, {}, {}, {}, {}, []
    for aa1 in seq :
        curindex = len(res)
        res[curindex] = {}
        resnums[curindex] = curindex+1
        resns[curindex] = AA13[aa1]
        chids[curindex] = " "
        inscodes[curindex] = " "
        resids[curindex] = makeResid( resns[curindex], chids[curindex], resnums[curindex], inscodes[curindex] )
        for an in resAtoms[ AA13[aa1] ] :
            res[curindex][an] = len(pts)
            pts.append( [-999., -999., -999.] )
    return res, resids, resnums, resns, chids, inscodes, pts

if __name__=="__main__" :
    #makeAllGLYmodel("/home/sg363/RTK/xloop/1KX8.pdb", "1KX8gly.pdb")
    main()
