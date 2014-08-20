import sys, os, shutil, re
from xray import cif2mtz, uniqueify, sfall, cns_generate, cns_anneal, mtz2hkl, sgCCP4toCNS
from xcheck import XrayRanker, XrayScorer


def main() :
    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='to create all the files during refinement. it shdnt be already present.')
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='starting pdb')
    parser.add_option("--ensembleSize", action='store', type='int', dest='esize', help='number of conformers in a multiconformer model', default=1)

    parser.add_option("--sf", action='store', type='string', dest='sf', help='structure factors file')
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
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='randomize will produce a different refinement trajectory by seeding rtk randomly', default=None)
    
    (options, args) = parser.parse_args()

    from loopbuild import Multiloop
    from scplacement import SCplacement
    import prepareChain
    import misc
    misc.setVerbosity(options.verbose)

    pref.randomize(options.randomize)

    if not os.path.isdir(options.scratchdir) : os.mkdir(options.scratchdir)
    for mi in range(options.esize) :
        shutil.copyfile(options.pdbfile, "%s/c0m%d.pdb" % (options.scratchdir,mi))
    shutil.copyfile(options.sf, "%s/strfactors.mtz" % options.scratchdir)
    os.chdir(options.scratchdir)
    cif2mtz("strfactors.mtz", "base.mtz", options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg)
    uniqueify("base.mtz", "rfree.mtz")

    from restraints import EDrestraint ; EDrestraint.setPenalty(5.) ; EDrestraint.setScatRad(1.) ;
    esmin, esmax, esmean, rcmult = 0.0001, 5.0, 0.05, 10
    from builders import PeptideBridgeBuilder ; PeptideBridgeBuilder.setCTtrials(25); PeptideBridgeBuilder.setThetastep(5);
    from data import consts ; consts.set("TAU_QUALITY", 40.)

    mconly, guidedsampling = None, None

    xscorer = XrayScorer(None, 0.9)
    startcycle, endcycle = 0, 10
    for cycle in range(startcycle, endcycle) :
        badresids, xrayRestGen, xranker = None, None, None
        multiPrepC = prepareChain.PrepareChain("PRL")
        rtkmodel, joinpdbs = "r%d.pdb"%(cycle), []
        for mi in range(options.esize) :
            cnsin = "c%dm%d.pdb" % (cycle,mi) ; rtkout = "r%dm%d.pdb" % (cycle,mi)
            if cycle > 0 :
                badresids = xscorer.score(cnsin, phasedmtz, "FP", "FC", "PHIC", "2F1-F2") ## assess bad fit
                #xrayRestGen = [ prepareChain.XrayRestraintsGenerator(phasedmtz, "FP", "FC", "PHIC", "2F1-F2", esmin, esmax, esmean ) ]
                xranker = XrayRanker(phasedmtz, "FP", "FC", "PHIC", "2F1-F2", esmin, esmax)
                xranker.rankChildren = 10 ; xranker.rankRecurse = 1 ; xranker.rankGivenEnsemble = None
            ml = Multiloop(cnsin, badresids, mconly, options.caRad, options.scRad, options.scReduction, guidedsampling, options.popsize,
                options.backtrack, 1, "pre."+rtkout, xrayRestGen, multiPrepC)
            ml.ranker = xranker
            nb = ml.run()
            print mi, "multiloop built", nb, "model/s" ; assert nb > 0
            if cycle == 0 : os.rename("pre."+rtkout, rtkout)
            else :
                badresids = findChangedSC(cnsin, "pre."+rtkout)
                scPrepC = prepareChain.PrepareChain("PRL")
                SCplacement("pre."+rtkout, options.scReduction, rtkout, "dotfile", None, phasedmtz, "FP", "FC", "PHIC", "2F1-F2", esmin, esmax, None, 1, badresids, scPrepC).run()
                adjustBfac(rtkout, cnsin)
            removeRElines(rtkout, ["^MODEL", "^ENDMDL",])
            copyNonprotein(cnsin, rtkout)
            joinpdbs.append(rtkout)
        joinPDBs(rtkmodel, joinpdbs) ## all models pasted together, with occupancy fractionalized, extra ENDs removed, segids corrected
        sfall(rtkmodel, "rfree.mtz", "phased.mtz")

        phasedmtz, cnsout = "phased%d.mtz"%(cycle+1), "c%d.pdb"%(cycle+1),
        cnsRefinement("phased.mtz", rtkmodel, phasedmtz, cnsout,
            options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg, options.resolution, None,cycle)

        newpdbs = []
        for mi in range(options.esize) : newpdbs.append("c%dm%d.pdb" % (cycle+1,mi))
        splitPDBs(cnsout, newpdbs)
        for newpdb in newpdbs : restoreChainids(newpdb, "c0m0.pdb")

## the resids betn refpdb and pdbfn are identical except chid. change chid to that in refpdb
def restoreChainids(pdbfn, refpdb) :
    print "Restoring chain ids from", refpdb, "to", pdbfn
    from pdbr import protein, changeChid, line2resnum, line2resn, line2inscode, line2chid
    import prot2res
    prot = protein(pdbfn, read_hydrogens=0, read_waters=1, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    prot = protein(refpdb, read_hydrogens=0, read_waters=1, read_hets=1)
    refres, refresids, refresnums, refresns, refchids, refinscodes, refpts = prot2res.readProtRes(prot)
    #print len(chids) , len(refchids)
    #for v in resids.values() : print "TEST", v
    #for v in refresids.values() : print "REF ", v
    assert len(chids) == len(refchids)
    for ki in resids.keys() :
        print "RESTORE", resnums[ki] , refresnums[ki], resns[ki] , refresns[ki], inscodes[ki] , refinscodes[ki]
        assert resnums[ki] == refresnums[ki]
        assert resns[ki] == refresns[ki]
        assert inscodes[ki] == refinscodes[ki]
    curki = 0 ; newlines = []
    for l in open(pdbfn, 'r').readlines() :
        if not l[0:6] in ["ATOM  ","HETATM"] : newlines.append(l) ; continue
        if line2resnum(l) != refresnums[curki] : curki += 1

        #print refresnums
        #print line2resnum(l)
        #print refresnums[curki]
        #print line2resn(l)
        #print refresns[curki]
        #print line2inscode(l)
        #print refinscodes[curki] 

        assert line2resnum(l) == refresnums[curki]
        assert line2resn(l) == refresns[curki]
        assert line2inscode(l) == refinscodes[curki]
        newlines.append(changeChid(l, refchids[curki]))
    fp = open(pdbfn, 'w')
    curch = None
    for l in newlines :
        if curch == None : curch = line2chid(l)
        elif l[0:6] in ["ATOM  ","HETATM"] and curch != line2chid(l) : print >> fp, "TER" ; curch = line2chid(l)
        fp.write(l)
    print >> fp, "TER\nEND"
    fp.close()

## split based on segids
def splitPDBs(pdbfn, pdbs) :
    from pdbr import line2segid
    segilines, segord = {}, []
    for l in open(pdbfn,'r').readlines() :
        if not l[0:6] in ["ATOM  ","HETATM"] : continue
        segi = line2segid(l)
        if not segilines.has_key(segi) : segilines[segi] = []; segord.append(segi)
        segilines[segi].append(l)
    assert len(segord) == len(pdbs)
    for si in range(len(segord)) :
        segi, outpdb = segord[si], pdbs[si]
        ofp = open(outpdb, 'w')
        for l in segilines[segi] : ofp.write(l)
        ofp.close()

## remove MODEL,ENDMDL records, change segids
## assume there are only ATOM,HETATM,MODEL,ENDMDL,TER records in the pdbs and each pdb contains a model
def joinPDBs(finalpdb, pdbs) :
    from pdbr import changeOccupancy, changeSegid
    newlines = []
    for pi in range(len(pdbs)) :
        lines = open(pdbs[pi],'r').readlines()
        for l in lines :
            l = re.sub("\n", "", l)
            if l[0:5] == "MODEL" or l[0:3] == "END" : continue
            if l[0:6] in ["ATOM  ","HETATM"] : l = changeOccupancy(l, 1./len(pdbs))
            if l[0:6] in ["ATOM  ","HETATM"] : l = changeSegid(l, pi+1)
            newlines.append(l)
    fp = open(finalpdb, 'w')
    for l in newlines : print >> fp, l
    print >> fp, "END"
    fp.close()

if __name__ == "__main__" :
    main() ; sys.exit(0)
    adjustBfac("r0m0.pdb", "c0m0.pdb")
    copyNonprotein("c0m0.pdb", "r0m0.pdb")
