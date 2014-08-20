def makeParser() :
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='template PDB file')
    parser.add_option("--mconly", action='store', type='int', dest='mconly', help='build mainchain-only models. DEFAULT 0', default=0)
    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position. DEFAULT 1', default=1)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid. DEFAULT 2', default=2)
    parser.add_option("--num-models-wanted", action='store', type='int', dest='nmodels', help='number of models desired, 100 attempts per model. DEFAULT 100', default=100)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy. DEFAULT 100', default=100)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw radius in case of sidechain clashchecks. DEFAULT 1 ie no reduction', default=1)
    parser.add_option("--guided-sampling", action='store', type='int', dest='guidedSampling', help='whether to use guided phipsi-omega sampling or usual propensity-weighted. DEFAULT 0', default=0)
    parser.add_option("--buildN2C", action='store', type='int', dest='buildN2C', help='by default, build from N to C terminal. build C->N if 0.', default=1)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--outpdb", action='store', type='str', dest='outpdb', help='file to write the models to, must be set to valid filepath')
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    return parser


def parseOptions(parser) :
    (options, args) = parser.parse_args()

    ### DANGER commented
    #if options.mconly != 1 : options.mconly = None

    
    #if options.guidedSampling == 1 : options.guidedSampling = options.caRad
    #else : options.guidedSampling = None

    import misc
    misc.setVerbosity(options.verbose)

    return options


def addXrayOptions(parser) :
    parser.add_option("--mtz", action='store', type='string', dest='mtzfn', help='mtz file', default=None)
    parser.add_option("--f1label", action='store', type='string', dest='f1label', help='F1 column label DEFAULT FP', default='FP')
    parser.add_option("--f2label", action='store', type='string', dest='f2label', help='F2 column label DEFAULT FC', default='FC')
    parser.add_option("--philabel", action='store', type='string', dest='philabel', help='Phase column label DEFAULT PHIC', default='PHIC')
    parser.add_option("--maptype", action='store', type='string', dest='maptype', help='type of map. can be F1 or the default 2F1-F2', default='2F1-F2')
    parser.add_option("--sigXmin", action='store', type='float', dest='sigXmin', help='electron density cutoff below which score is negative. DEFAULT 0.25', default=0.25)
    parser.add_option("--sigXmax", action='store', type='float', dest='sigXmax', help='electron density cutoff above which score remains constant at max. DEFAULT 2.0', default=2.0)
    parser.add_option("--sigXmean", action='store', type='float', dest='sigXmean', help='mean electron density restraint cutoff. DEFAULT 1.0', default=1.0)
    return parser
