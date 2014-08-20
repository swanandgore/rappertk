import os, shutil, re
from xray import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS

def removeZeroLines(filename) :
    lines = open(filename, 'r').readlines()
    fp = open(filename, 'w')
    for l in lines :
        if re.compile('ATOM.* 0.000.* 0.000.* 0.000').search(l) : continue
        if re.compile('HETATM.* 0.000.* 0.000.* 0.000').search(l) : continue
        if re.compile('HETATM.*9999.*9999.*9999').search(l) : continue
        if re.compile('ATOM.*9999.*9999.*9999').search(l) : continue
        fp.write(l)
    fp.close()

def main() :
    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='to create all the files during refinement. it shdnt be already present.')
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='starting pdb containing a model of pdb-ligand complex')
    parser.add_option("--sf", action='store', type='string', dest='sf', help='structure factors file')
    parser.add_option("--topCNS", action='store', type='string', dest='topCNS', help='ligand topology file for CNS', default=2)
    parser.add_option("--parCNS", action='store', type='string', dest='parCNS', help='ligand parameter file for CNS', default=2)
    parser.add_option("--a", action='store', type='float', dest='a', help='cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='cell spacegroup, in CCP4 notation')

    parser.add_option("--ligfile", action='store', type='string', dest='ligfile', help='ligand description, see msq.ligdesc for description of MSQ in 1di9.pdb')
    parser.add_option("--around-ligand", action='store', type='float', dest='closeCutoff', help='min-dist between residue/ligand to be considered close to ligand', default=10)
    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid', default=2)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=1)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)

    parser.add_option("--folabel", action='store', type='string', dest='folabel', help='Fo column label, defaults to FP', default='FP')
    parser.add_option("--fclabel", action='store', type='string', dest='fclabel', help='Fc column label, defaults to FC', default='FC')
    parser.add_option("--philabel", action='store', type='string', dest='philabel', help='Phase column label, defaults to PHIC', default='PHIC')
    parser.add_option("--maptype", action='store', type='string', dest='maptype', help='type of map. can be Fc or Fo or the default 2Fo-Fc', default='2Fo-Fc')

    (options, args) = parser.parse_args()

    import misc
    misc.setVerbosity(options.verbose)

    os.mkdir(options.scratchdir)
    shutil.copyfile(options.pdbfile, "%s/model0.pdb" % options.scratchdir)
    shutil.copyfile(options.sf, "%s/strfactors.mtz" % options.scratchdir)
    os.chdir(options.scratchdir)

    cif2mtz("strfactors.mtz", "base.mtz", options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg)
    uniqueify("base.mtz", "rfree.mtz")
    numRefCycles = 20
    for cycle in range(numRefCycles) :
        print "***********************************************"
        print "********************CYCLE %d*******************" % cycle
        print "***********************************************"
        from xray import sgCCP4toCNS
        modelIn = "model%d.pdb" % cycle
        phasedmtz = "phased%d.mtz" % cycle
        ## create a phased 2Fo-Fc map from current model and str factors
        sfall(modelIn, "rfree.mtz", phasedmtz)
        ## generate a rappertk model
        rtkmodel = "rtk%d.pdb" % cycle
        import prepareChain, protlig
        xrayRestGen = prepareChain.XrayRestraintsGenerator(phasedmtz, "FP", "FC", "PHIC", "2Fo-Fc", 0., 1.5, 1.)
        protlig.main(modelIn, options.ligfile, options.closeCutoff, options.caRad, options.scRad, options.scReduction, None,
            rtkmodel, options.popsize, "3X3", 1,
            xrayRestGen)
        #shutil.copyfile(modelIn, rtkmodel)
        ## run cns refinement
        mtz2hkl(phasedmtz, "cns.hkl")
        cns_generate(rtkmodel, "generate.mtf", "generate.pdb", options.topCNS, options.parCNS, "generate%d.log"%cycle)
        removeZeroLines("generate.pdb")
        cns_anneal(options.a, options.b, options.c, options.alpha, options.beta, options.gamma, sgCCP4toCNS[options.sg], "cns.hkl", "generate.mtf", "generate.pdb", options.parCNS, "anneal%d.log"%cycle)
        removeZeroLines("anneal_1.pdb")
        os.rename("anneal_1.pdb", "model%d.pdb"%(cycle+1))

if __name__ == "__main__" : main()
