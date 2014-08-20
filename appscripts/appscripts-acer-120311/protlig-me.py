from pdbr import protein
import prot2res
import re, string, sys   , os , shutil , random
from geometry import calcDist, CAtraceGH, MapIntIntFloat
from builders import VecInt, VecFloat, VecBuilder, VecVecFloat, Rotator, TransRotator, BuilderGroup
from restraints import DistanceRestraint
from data import vdwr, consts, PROBE_DISULFIDE_OVERLAP_MARGIN
from prepareChain import  incompleteSCcorrection
import data
from data import vdwr, consts
from prefRapper import randomize
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
                if v == flds[1] :
                    ligname = flds[1] ; ligind = k ;
                    print "res-lig" , res[ligind]
                    break
            assert ligind != None
        elif l[0:7] == 'rotbond' :
            min, max, step = string.atof(flds[1]), string.atof(flds[2]), string.atof(flds[3]),
            print min,max , step
            IPs, OPs = [ res[ligind][anames[0]],res[ligind][anames[1]] ], []
            print "IPS", IPs

            for an in anames[2:] :
                OPs.append(res[ligind][an])
            builders.append( Rotator(VecInt(IPs), VecInt(OPs), None, "Rotator", min, max, step) )
            print "OPs" , OPs
        elif l[0:4] == 'init' :
            tol1, tol2 = string.atof(flds[1]), string.atof(flds[2])
            init1, init2 = anames[0], anames[1]
        else : print "Unrecognized file format for %s, exiting....." % ligfile ; assert None
    sys.exit()
    assert tol1
    IPs, OPs = [], [res[ligind][init1], res[ligind][init2]]
    allAnames.remove(init1) ;
    allAnames.remove(init2)
    for an in allAnames :
        OPs.append(res[ligind][an])
    builders = [ TransRotator(VecInt(IPs),VecInt(OPs),None,"TransRotator", VecFloat(pts[res[ligind][init1]]),tol1, VecFloat(pts[res[ligind][init2]]),tol2) ] + builders
    return builders, restraints, ligname

def main(pdbfile, ligfile, closeCutoff, caRad, scRad, scReduction, guidedSampling, outpdb, popsize, backtrack, nmodels, restrGen,ranker,scl="PRL") :

    prot = protein(pdbfile, read_hydrogens=0, read_waters=1, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    ## find ligand name
    ligname = None
    for l in open(ligfile,'r').readlines() :
        if l[0:7] == "ligname" : ligname = l[8:11]
    assert ligname


    ## apply missing sc correction
    scMissInds = incompleteSCcorrection(res, resns, pts)
    for i in scMissInds :
        print "SCMISSINDS", resids[i]

    ## find ligand's min-dist from other residues. find contiguous segments entirely within 10A distance of ligand
    leastDists = findLigDists(ligname, res, pts, resns)
    closeRes = [None] * len(leastDists)
    for ind,d in leastDists.items() :
        if d < closeCutoff :
            closeRes[ind] = 1

    segstart, segs = None, []
    for ind in range(len(closeRes)) :
        if closeRes[ind] and not segstart : segstart = ind
        elif closeRes[ind] == None and segstart :
            segs.append([segstart,ind-1])
            segstart = None
    if segstart : segs.append([segstart,len(closeRes)-1])

    ## correct segments that are separated by 1 residue
    #print segs
    segs2remove = []
    for si in range(1,len(segs)) :
        if segs[si][0] - segs[si-1][1] <= 2 :
            if segs[si][1]-segs[si][0] == 0 and segs[si-1][1]-segs[si-1][0] == 0 : # both lengths 1, merge segs
                segs[si][0] = segs[si-1][0]

                segs2remove.append(si-1)
                si += 1
            elif segs[si][1]-segs[si][0] == 0 : # shrink si-1
                segs[si-1][1] -= 1
            elif segs[si-1][1]-segs[si-1][0] == 0 : # shrink si
                segs[si][0] += 1
            else : segs[si][0] += 1
    newsegs = []
    for si in range(len(segs)) :
        if si in segs2remove : continue
        newsegs.append(segs[si])

    segs = newsegs
    badresids = []
    for seg in segs:
        for i in range(seg[0],seg[1]+1) :
            badresids.append(resids[i])


    if len(badresids) == 0:
        print "No residues in contact within specified cutoff %f Angstroms"% closeCut
        import sys; sys.exit()
            
    from loopbuild import locateRegionsRandomize
    loops, looptypes = locateRegionsRandomize(resids, chids, badresids) # locate loops and order them randomly


    blist, numtrials, bfoll = [], [], {}
    irlist,hydlist =  [], []
    dummies = {}
    

    bldrs, irlist, ligname = makeLigBuilder(res, resns, pts, ligfile)
    blist = [ BuilderGroup(VecBuilder(bldrs), "LigandBuilder") ]
    numtrials.append(1000)

    from prepareChain import  addNCdummyGly, removeSC, mergeBRlists, makeCApropRestraints , PrepareChain
    prepC = PrepareChain(scl)

        
    for li in range(len(loops)) :
        startindex, endindex = loops[li]
        chid = chids[startindex] ; assert chid == chids[endindex]

        bl,bf,nt,rl,dum, reorderBuilders = None,None,None,None,None, 1



        if looptypes[li] == "NtoC" :
            if resns[ loops[li][0] ] in [ "  A", "  T", "  C", "  G", "  U", ] :
                
                bl,nt,bf,rl,dum = prepareRNAchain(chid, res, resids, resnums, resns, chids, inscodes, pts, caRad, scRad)
            else : #### modified amk
                bl,nt,bf,rl,dum = prepC.preparePeptideChain(chid, res, resids, resnums, resns, chids, inscodes, pts, None , buildfwd, guidedSampling, caRad, scRad, scMissInds)
                    

        elif looptypes[li] == "Nter" : #### modified amk
            from prepareChain import addNdummyGly
            dum, firstindex = addNdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts)
            bl,bf,nt,rl = prepC.prepareChainTerminal("Nterm", startindex, endindex, firstindex, res, resids, resnums, resns, chids, inscodes, pts, None, caRad, scRad, guidedSampling, scMissInds)

        elif looptypes[li] == "Cter" : #### modified amk
            from prepareChain import addCdummyGly
            dum, lastindex = addCdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts)
            bl,bf,nt,rl = prepC.prepareChainTerminal("Cterm", startindex, endindex, lastindex, res, resids, resnums, resns, chids, inscodes, pts, None, caRad, scRad, guidedSampling, scMissInds)

        elif looptypes[li] == "loop" : #### modified amk
            if endindex-startindex < 7 and random.random() > 1.5 :
                bl,bf,nt,rl = prepC.preparePeptideLoop1(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, None, caRad, scRad, scMissInds, guidedSampling)


            else : #### modified amk
                bl,bf,nt,rl = prepC.preparePeptideMidloop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, None, caRad, scRad, scMissInds, guidedSampling)
                reorderBuilders = None
        else :
            print "Dont understand looptype", looptypes[li] ; assert None

        mergeBRlists(blist, bfoll, numtrials, irlist, bl, bf, nt, rl, dummies, dum)





    from prepareChain import printResAtoms
    printResAtoms(res,resids)

    ai2pos, aiSC = [-999] * len(pts), [-999] * len(pts) ## flag atoms
    for index,val in res.items() :
        for name, ai in val.items() :
            ai2pos[ai] = index
            if name == ' N  ' : aiSC[ai] = 0
            elif name == ' CA ' : aiSC[ai] = 1
            elif name == ' C  ' : aiSC[ai] = 2
            elif name == ' O  ' : aiSC[ai] = 3
            elif name == ' SG ' and resns[index] == 'CYS' : aiSC[ai] = 5 # for disulphide vdw reduction
            elif name == ' CD ' and resns[index] == 'PRO' : aiSC[ai] = 6 # PRO CD shdnt be clash-checked with prev res's C
            else : aiSC[ai] = 4
    for k in dummies.keys() :
        for an in dummies[k] : aiSC[res[k][an]] = -1

    print "---------------------------------------\n\n\n"        


    knownPositions = range(len(pts))
    for b in blist :
        bop = b.getOP()
        for i in range(bop.size()) :
            if bop[i] in knownPositions :
                knownPositions.remove(bop[i])




    radii = [0] * len(pts) ## make radii for all atoms including dummies
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else :
                try : radii[pi] = vdwr['XXX'][aname[1]]
                except KeyError : radii[pi] = vdwr['XXX']['C']
    assert len(radii) == len(pts)
    pts = VecVecFloat(pts)
    print "---------------------------------------\n\n\n"


    ssReduction = data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    if restrGen !=None :
        xrlist,xopt= restrGen.generate(blist, aiSC, res, resids, pts) 
        
    else :
        xrlist = [] ; xopt = [] 

        

    for r in irlist:
        xrlist.append(r)

    print "NOW PRINTING ALL RESTRAINTS:: "

    for r in xrlist:
        r.describe();
        print 
        


    from peptidebuild import ModelRenderer
    from PopulationStrategy import PopulationStrategy
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, dummies.keys(), outpdb)

    strategy = PopulationStrategy(backtrack,popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll,  None , xrlist , xopt , modelRenderer, nmodels*100, nmodels)

    strategy.ranker = ranker
    strategy.execute()
    import sys  ; sys.exit()


def callmain() :
    from commonOptions import makeParser, parseOptions, addXrayOptions
    parser = makeParser()
    parser = addXrayOptions(parser)
    parser.add_option("--ligfile", action='store', type='string', dest='ligfile', help='ligand description, see msq.ligdesc for description of MSQ in 1di9.pdb')
    parser.add_option("--around-ligand", action='store', type='float', dest='closeCutoff', help='min-dist between residue/ligand to be considered close to ligand', default=10)
    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='to create all the files during refinement. it shdnt be already present.',default = "rappertkout")
    parser.add_option("--mapfn1", action='store', type='string', dest='mapfn1', help='ccp4 map restraining the    shape of the molecule', default=None)

    parser.add_option("--a", action='store', type='float', dest='a', help='cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='cell spacegroup, in CCP4 notation')
    parser.add_option("--resolution", action='store', type='float', dest='resolution', help='resolution of the data')
    parser.add_option("--make-ed-optional", action='store', type='int', dest='edOpt', help='    If = 1, then the 0.0 (positive density) mainchain restraint will be made optional. If false, then the main chain will be unconditionally forced to lie in positive density. This is primarily useful when tracing through a structure with regions in very poor (non-existent) density', default=None)    
    
    parser.add_option("--sclib", action='store', type='string', dest='sclib', help='Sidechain library to use: PRL or SCL0.2 or SCL0.5 , SCL1.0 ', default="SCL1.0")
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='seed for randomizing', default=None)    

    parser.remove_option("--mconly")
    parser.remove_option("--buildN2C")
    parser.remove_option("--guided-sampling")
    options = parseOptions(parser)

    import misc
    misc.setVerbosity(options.verbose)
    from prefRapper import randomize
    randomize(options.randomize)
    



    if not os.path.isdir(options.scratchdir) :
        os.mkdir(options.scratchdir)

    if (os.path.isfile(options.pdbfile)==False) :
        print "Cannot find file %s "%options.pdbfile
        print "No file in directory ", os.getcwd()
        import sys ; 
        sys.exit()

    shutil.copyfile(options.pdbfile, "%s/%s" % (options.scratchdir,options.pdbfile))


    if (os.path.isfile(options.ligfile)==False) :
        print "Cannot find file %s "%options.ligfile
        print "No file in directory ", os.getcwd()
        import sys ; 
        sys.exit()

    shutil.copyfile(options.ligfile, "%s/%s" % (options.scratchdir,options.ligfile))



    ## PRINTING ALL PARAM values
    print "--ligfile = ",options.ligfile
    print "--around-ligand = ",options.closeCutoff
    print "--mapfn1 = ",options.mapfn1
    print "--pdb = ",options.pdbfile
    print "--ca-restraint-radius = ",options.caRad
    print "--sc-centroid-restraint-radius = ",options.scRad
    print "--num-models-wanted = ",options.nmodels
    print "--population-size = ",options.popsize
    print "--sidechain-vdw-reduction = ",options.scReduction
    #print "--guided-sampling = ",options.guidedSampling
    print "--backtrack = ",options.backtrack
    print "--outpdb = ",options.outpdb
    print "--verbose = ",options.verbose
    print "--randomize = ",options.randomize
    print "--mtz = ",options.mtzfn
    print "--f1label = ",options.f1label
    print "--f2label = ",options.f2label
    print "--philabel = ",options.philabel
    print "--maptype = ",options.maptype

    print "--a = ",options.a
    print "--b = ",options.b
    print "--c = ",options.c
    print "--alpha = ",options.alpha
    print "--beta = ",options.beta
    print "--gamma = ",options.gamma
    print "--sg = ",options.sg
    print "--make-ed-optional = ",options.edOpt
    print "--resolution = ",options.resolution
    print "--randomize = ",options.randomize
    
    print "--sigXmin = ",options.sigXmin
    print "--sigXmax = ",options.sigXmax
    print "--sigXmean = ",options.sigXmean    
    print "--sclib = ",options.sclib

    fp = open('parameters.txt', 'w')
    l =  "--ligfile = %s "%str(options.ligfile) ;      print >> fp, l 
    l =  "--around-ligand = %s "%str(options.closeCutoff) ;      print >> fp, l 
    l =  "--mapfn1 = %s "%str(options.mapfn1) ;      print >> fp, l 
    l =  "--pdb = %s "%str(options.pdbfile) ;      print >> fp, l 
    l =  "--mconly = %s "%str(options.mconly) ;      print >> fp, l 
    l =  "--ca-restraint-radius = %s "%str(options.caRad) ;      print >> fp, l 
    l =  "--sc-centroid-restraint-radius = %s "%str(options.scRad) ;      print >> fp, l 
    l =  "--num-models-wanted = %s "%str(options.nmodels) ;      print >> fp, l 
    l =  "--population-size = %s "%str(options.popsize) ;      print >> fp, l 
    l =  "--sidechain-vdw-reduction = %s "%str(options.scReduction) ;      print >> fp, l 
    #l =  "--guided-sampling = %s "%str(options.guidedSampling) ;      print >> fp, l 
    l =  "--backtrack = %s "%str(options.backtrack) ;      print >> fp, l 
    l =  "--outpdb = %s "%str(options.outpdb) ;      print >> fp, l 
    l =  "--verbose = %s "%str(options.verbose) ;      print >> fp, l 
    l =  "--mtz = %s "%str(options.mtzfn) ;      print >> fp, l 
    l =  "--f1label = %s "%str(options.f1label) ;      print >> fp, l 
    l =  "--f2label = %s "%str(options.f2label) ;      print >> fp, l 
    l =  "--philabel = %s "%str(options.philabel) ;      print >> fp, l 
    l =  "--maptype = %s "%str(options.maptype) ;      print >> fp, l 
    l =  "--sigXmin = %s "%str(options.sigXmin) ;      print >> fp, l 
    l =  "--sigXmax = %s "%str(options.sigXmax) ;      print >> fp, l 
    l =  "--sigXmean = %s "%str(options.sigXmean) ;      print >> fp, l 
    l =  "--sclib = %s "%str(options.sclib) ;      print >> fp, l 
    l =  "--randomize = %s "%str(options.randomize) ;      print >> fp, l 

    l =  "--a = %s "%str(options.a) ;      print >> fp, l 
    l =  "--b = %s "%str(options.b) ;      print >> fp, l 
    l =  "--c = %s "%str(options.c) ;      print >> fp, l 
    l =  "--alpha = %s "%str(options.alpha) ;      print >> fp, l 
    l =  "--beta = %s "%str(options.beta) ;      print >> fp, l 
    l =  "--gamma = %s "%str(options.gamma) ;      print >> fp, l 
    l =  "--sg = %s "%str(options.sg) ;      print >> fp, l 
    l =  "--make-ed-optional = %s "%str(options.edOpt) ;      print >> fp, l 
    l =  "--resolution = %s "%str(options.resolution) ;      print >> fp, l 


    fp.close()
    ###
    

    from prepareChain import XrayRestraintsGenerator
    from xcheck import XrayScorer, XrayRanker


    if options.mtzfn !=None or options.mapfn1 !=None : 
        if (options.a == None or options.b == None or options.c == None or options.alpha== None or options.beta==None or options.gamma == None or options.sg == None or options.resolution == None ):

            if options.resolution==None :
                print "Please check value of resolution = ", options.resolution
                import sys ; 
                sys.exit()
            
            print "Please check cell paramater a, b , c , alpha, beta , gamma = ",options.a , options.b , options.c , options.alpha , options.beta  , options.gamma 
            print "Please check space group " , options.sg
            import sys ; 
            sys.exit()




    if options.mapfn1 != None:
        if (os.path.isfile(options.mapfn1)==False) :
            print "Cannot find file %s "%options.mapfn1
            print "No file in directory ", os.getcwd()
            import sys ; 
            sys.exit()

        shutil.copyfile(options.mapfn1, "%s/%s" % (options.scratchdir,options.mapfn1))




    if options.mtzfn != None:
        if (os.path.isfile(options.mtzfn)==False) :
            print "Cannot find file %s "%options.mtzfn
            print "No file in directory ", os.getcwd()
            import sys ; 
            sys.exit()

        shutil.copyfile(options.mtzfn, "%s/%s" % (options.scratchdir,options.mtzfn))



    xrayRestGen = None ; ranker = None
    if options.mtzfn != None : 
        xrayRestGen = XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, options.sigXmin, options.sigXmax, options.sigXmean, [] , options.edOpt ) 
        shutil.copyfile(options.mtzfn, "%s/%s" % (options.scratchdir,options.mtzfn))

        #ranker = None
        #esmin, esmax, esmean, rcmult = , 5 
        #ranker = XrayRanker("phased.mtzFPFCPHIC_A.map", "map", "FC", "PHIC", "2F1-F2", esmin, esmax)
        #ranker.rankChildren = rcmult ; ranker.rankRecurse = 1
        #ranker.rankLeaderBuilderOnly = None ; ranker.rankGivenEnsemble = None


    elif options.mapfn1 != None  :
        xrayRestGen = XrayRestraintsGenerator(options.mapfn1, "map", options.f2label, options.philabel, options.maptype,  options.sigXmin,  options.sigXmax, options.sigXmean , [] , options.edOpt ) 
        shutil.copyfile(options.mapfn1, "%s/%s" % (options.scratchdir,options.mapfn1))


    else :
        print "WARNING !! No input mtz or map file , No electron density restraints will be used"



    os.chdir(options.scratchdir)

#    main(options.pdbfile, options.ligfile, options.closeCutoff, options.caRad, options.scRad, options.scReduction, options.guidedSampling,options.outpdb, options.popsize, options.backtrack, options.nmodels,xrayRestGen, options.mapfile,ranker,options.sclib)

    main(options.pdbfile, options.ligfile, options.closeCutoff, options.caRad, options.scRad, options.scReduction, None,options.outpdb, options.popsize, options.backtrack, options.nmodels,xrayRestGen,ranker,options.sclib)

if __name__ == "__main__" : callmain()
