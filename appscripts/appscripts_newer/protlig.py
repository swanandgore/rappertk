from pdbr import protein
import prot2res
import re, string, sys
from geometry import calcDist, CAtraceGH, MapIntIntFloat
from builders import VecInt, VecFloat, VecBuilder, VecVecFloat, Rotator, TransRotator, BuilderGroup
from restraints import DistanceRestraint
from data import vdwr, consts, PROBE_DISULFIDE_OVERLAP_MARGIN
from prepareChain import preparePeptideLoop, incompleteSCcorrection

def findLigDists(ligname, res, pts, resns) :
    ligdists = {}
    ligind = None
    for ind in res.keys() :
        if resns[ind] == ligname : ligind = ind ; break
    assert ligind
    for ind in res.keys() :
        if resns[ind] == ligname : continue
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

def main(pdbfile, ligfile, closeCutoff, caRad, scRad, scReduction, guidedSampling, outpdb, popsize, backtrack, nmodels, restrGen=None) :
    prot = protein(pdbfile, read_hydrogens=0, read_waters=1, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    ## find ligand name
    ligname = None
    for l in open(ligfile,'r').readlines() :
        if l[0:7] == "ligname" : ligname = l[8:11]
    assert ligname

    print "Residues and atom-numbers-----------------------"
    keys = res.keys() ; keys.sort()
    for k in keys : print "(%d)"%k, resids[k], res[k], chids[k]
    print "------------------------------------------------\n\n\n"

    ## apply missing sc correction
    scMissInds = incompleteSCcorrection(res, resns, pts)
    print "SCMISSINDS", scMissInds
    for i in scMissInds :
        print "SCMISSINDS", resids[i]

    ## find ligand's min-dist from other residues. find contiguous segments entirely within 10A distance of ligand
    leastDists = findLigDists(ligname, res, pts, resns)
    closeRes = [None] * len(leastDists)
    for ind,d in leastDists.items() :
        if d < closeCutoff : closeRes[ind] = 1
        print resids[ind], d
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
    print segs

    blist, numtrials, bfoll = [], [], {}
    knownPositions = range(len(pts))

    bldrs, irlist, ligname = makeLigBuilder(res, resns, pts, ligfile)
    blist = [ BuilderGroup(VecBuilder(bldrs), "LigandBuilder") ]
    for ind in res.keys() :
        if resns[ind] != ligname : continue
        for an,pi in res[ind].items() : knownPositions.remove(pi)
    numtrials.append(1000)

    ## for each close segment, setup loop sampling
    for seg in segs :
        startindex, endindex = seg[0], seg[1]
        if resns[startindex] == "HOH" : continue
        print "SEGMENT -------- [%s] [%s]" % (resids[startindex], resids[endindex]), seg
        if chids[startindex] != chids[endindex] :
            print "end-of-chain close-to-ligand regions are not yet supported" ; sys.exit(0)
        for ind in range(startindex-1,endindex+2) :
            for an,pi in res[ind].items() :
                if ind == startindex-1 and an in [' N  ',' CA '] : continue
                if ind == endindex+1 and an in [' CA ',' C  ',' O  '] : continue
                if pi in knownPositions : knownPositions.remove(pi)

        guidedRad = None
        if guidedSampling != None : guidedRad = caRad
        blist1, bfoll1, numtrials1, irlist1 = preparePeptideLoop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, None, caRad, scRad, scMissInds, guidedRad)
        blistsize = len(blist)
        for b in blist1 : blist.append(b)
        for bi in bfoll1.keys() :
            if not bfoll.has_key(bi+blistsize) : bfoll[bi+blistsize] = []
            for fi in bfoll1[bi] : bfoll[bi+blistsize].append(fi+blistsize)
        for n in numtrials1 : numtrials.append(n)
        for r in irlist1 : irlist.append(r)

    print "BUILDERS------------------------------"
    for b in blist :
        b.describe() ; print ''
    print "---------------------------------------\n\n\n"

    ## make radii
    radii = [0] * len(pts)
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else :
                try : radii[pi] = vdwr['XXX'][aname[1]]
                except KeyError : radii[pi] = vdwr['XXX']['C']
    pts = VecVecFloat(pts)

    ai2pos, aiSC = [-999] * len(pts), [-999] * len(pts)
    for index,val in res.items() :
        for name, ai in val.items() :
            ai2pos[ai] = index
            if name == ' N  ' : aiSC[ai] = 0
            elif name == ' CA ' : aiSC[ai] = 1
            elif name == ' C  ' : aiSC[ai] = 2
            elif name == ' O  ' : aiSC[ai] = 3
            elif name == ' S  ' and resns[index] == 'CYS' : aiSC[ai] = 5 # for disulphide vdw reduction
            else : aiSC[ai] = 4

    if restrGen :
        xrlist,optional = restrGen.generate(blist, aiSC)
        for oi in range(len(optional)) : optional[oi] += len(irlist)
        for xr in xrlist : irlist.append(xr)

    scReduction, ssReduction = scReduction, PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [], outpdb)

    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(backtrack, popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, 1, irlist, optional, modelRenderer, nmodels*100, nmodels)
    strategy.execute()

def callmain() :
    from commonOptions import makeParser, parseOptions, addXrayOptions
    parser = makeParser()
    parser = addXrayOptions(parser)
    parser.add_option("--ligfile", action='store', type='string', dest='ligfile', help='ligand description, see msq.ligdesc for description of MSQ in 1di9.pdb')
    parser.add_option("--around-ligand", action='store', type='float', dest='closeCutoff', help='min-dist between residue/ligand to be considered close to ligand', default=10)
    parser.remove_option("--mconly")
    options = parseOptions(parser)

    import misc
    misc.setVerbosity(options.verbose)

    xrayRestGen = None
    if options.mtzfn :
        from prepareChain import XrayRestraintsGenerator
        xrayRestGen = XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, options.sigXmin, options.sigXmax, options.sigXmean)

    main(options.pdbfile, options.ligfile, options.closeCutoff, options.caRad, options.scRad, options.scReduction, options.guidedSampling,
        options.outpdb, options.popsize, options.backtrack, options.nmodels,
        xrayRestGen)

if __name__ == "__main__" : callmain()
