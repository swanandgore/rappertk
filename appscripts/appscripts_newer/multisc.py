import sys, os, shutil, re
from xray import cif2mtz, uniqueify, sfall, cns_generate, cns_anneal, mtz2hkl, sgCCP4toCNS
from xcheck import XrayRanker, XrayScorer
from multProtref import joinPDBs, splitPDBs, restoreChainids
from pref import adjustBfac, copyNonprotein, cnsRefinement, randomize, removeRElines
from screfine import keepMClinesOnly, buildCBs

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
    parser.add_option("--extraTOP", action='store', type='string', dest='extraTOP', help='ligand CNS topology file')
    parser.add_option("--extraPAR", action='store', type='string', dest='extraPAR', help='ligand CNS params file')

    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid', default=2)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=1)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='randomize will produce a different refinement trajectory by seeding rtk randomly', default=None)
    
    (options, args) = parser.parse_args()

    randomize(options.randomize)

    from loopbuild import Multiloop
    from scplacement import SCplacement
    import prepareChain
    import misc
    misc.setVerbosity(options.verbose)

    if not os.path.isdir(options.scratchdir) : os.mkdir(options.scratchdir)
    joinpdbs = []
    for mi in range(options.esize) :
        firstmodel = "%s/c0m%d.pdb" % (options.scratchdir,mi)
        if options.caRad < 0.5 :
            shutil.copyfile(options.pdbfile, firstmodel)
        else : # catrace
            mlPrepC = prepareChain.PrepareChain("PRL")
            Multiloop(options.pdbfile, None, None, options.caRad, options.scRad, options.scReduction, None, options.popsize, options.backtrack, 1, firstmodel, None, mlPrepC).run()
        keepMClinesOnly(firstmodel) ; buildCBs(firstmodel) ; joinpdbs.append(firstmodel)
    joinPDBs(options.scratchdir+"/c0.pdb", joinpdbs)
    #sys.exit(0)

    shutil.copyfile(options.sf, "%s/strfactors.mtz" % options.scratchdir)
    os.chdir(options.scratchdir)
    cif2mtz("strfactors.mtz", "base.mtz", options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg)
    uniqueify("base.mtz", "rfree.mtz")

    from restraints import EDrestraint ; EDrestraint.setPenalty(5.) ; EDrestraint.setScatRad(1.) ;
    esmin, esmax, esmean, rcmult = 0.01, 3.0, 0.05, 10

    mconly, guidedsampling = None, None

    xscorer = XrayScorer(None, 0.9)
    startcycle, endcycle = 0, 10
    for cycle in range(startcycle, endcycle) :
        phasedmtz = "phased%d.mtz" % cycle ; sfall("c%d.pdb"%cycle, "rfree.mtz", phasedmtz)
        newpdbs = []
        for mi in range(options.esize) : newpdbs.append("c%dm%d.pdb" % (cycle,mi))
        splitPDBs("c%d.pdb"%cycle, newpdbs)
        if cycle == 0 :
            for newpdb in newpdbs : restoreChainids(newpdb, "c0.pdb")
        else :
            for newpdb in newpdbs : restoreChainids(newpdb, options.pdbfile)
        rtkmodel, joinpdbs = "r%d.pdb"%(cycle), []
        for mi in range(options.esize) :
            cnsin = "c%dm%d.pdb" % (cycle,mi) ; rtkout = "r%dm%d.pdb" % (cycle,mi)
            if cycle > 0 : addSC, badresids = None, xscorer.score(cnsin, phasedmtz, "FP", "FC", "PHIC", "2F1-F2", "sc")
            else : addSC, badresids = 1, None
            if cycle > 5 : useGivenRot = 1
            else : useGivenRot = None
            scPrepC = prepareChain.PrepareChain("PRL")
            print "SC begin"
            dee = 1
            SCplacement(cnsin, options.scReduction, rtkout, "dotfile", dee, phasedmtz, "FP", "FC", "PHIC", "2F1-F2", esmin, esmax, addSC, useGivenRot, badresids, scPrepC).run()
            print "SC done"
            if cycle > 0 : adjustBfac(rtkout, cnsin)
            removeRElines(rtkout, ["^MODEL", "^ENDMDL",])
            if cycle == 0 : copyNonprotein(options.pdbfile, rtkout)
            else : copyNonprotein(cnsin, rtkout)
            joinpdbs.append(rtkout)
        joinPDBs(rtkmodel, joinpdbs) ## all models pasted together, with occupancy fractionalized, extra ENDs removed, segids corrected
        sfall(rtkmodel, "rfree.mtz", "phased.mtz")
        phasedmtz, cnsout = "phased%d.mtz"%(cycle+1), "c%d.pdb"%(cycle+1),
        cnsRefinement("phased.mtz", rtkmodel, phasedmtz, cnsout,
            options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg, options.resolution,
            None, cycle, options.extraTOP, options.extraPAR)


if __name__ == "__main__" :
    main() ; sys.exit(0)
    adjustBfac("r0m0.pdb", "c0m0.pdb")
    copyNonprotein("c0m0.pdb", "r0m0.pdb")
