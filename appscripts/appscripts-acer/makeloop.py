def main() :
    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='given pdb')
    parser.add_option("--outpdb", action='store', type='string', dest='outpdb', help='write model out here')

    parser.add_option("--a", action='store', type='float', dest='a', help='cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='cell spacegroup, in CCP4 notation')

    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid', default=2)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=1)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='seed for randomizing', default=None)
    parser.add_option("--loopres", action='store', type='string', dest='loopres', help='filename containing resids for starting perturbation', default=None)

    (options, args) = parser.parse_args()
    import misc
    misc.setVerbosity(options.verbose)

    loopresids = []
    for l in open(options.loopres, 'r').readlines() : loopresids.append(l[1:len(l)-2])
    options.loopres = loopresids

    from loopbuild import Multiloop
    from data import sgtable
    nmodels = 1
    ml = Multiloop(options.pdbfile, options.loopres, None, options.caRad, options.scRad, options.scReduction, None, options.popsize,
        options.backtrack, nmodels, options.outpdb, None, None)
    ml.cellsym = [ options.a, options.b, options.c, options.alpha, options.beta, options.gamma, sgtable[options.sg][0] ]
    nb = ml.run()

if __name__ == "__main__" : main()
