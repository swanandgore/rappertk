from pdbr import protein
import prot2res
from prepareChain import preparePeptideLoop, prepareRNAchain, incompleteSCcorrection
from geometry import calcDist, CAtraceGH
from builders import VecFloat, VecVecFloat, VecInt
from data import vdwr, PROBE_DISULFIDE_OVERLAP_MARGIN


def findRNAdists(rnaChains, res, pts, chids) :
    '''returns min dist betn a CA and all C4*s'''
    rnaPis, leastDists = [], {}
    for ind in res.keys() :
        if chids[ind] in rnaChains : rnaPis.append(res[ind][' C4*'])
    for ind in res.keys() :
        if chids[ind] in rnaChains : continue
        if not res[ind].has_key(' CA ') : continue
        mindist = 1e10
        for rpi in rnaPis :
            dist = calcDist( VecFloat(pts[res[ind][' CA ']]), VecFloat(pts[rpi]) )
            if dist < mindist : mindist = dist
        leastDists[ind] = mindist
        assert mindist < 1e9
    return leastDists

def main(pdbfile, pchids, rchids, closeCutoff, rnaRad, caRad, scRad, scReduction, outpdb, nmodels, popsize, backtrack, restrGen) :
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    ## check that reqd rna chids are present
    rch, pch = set(list(rchids)), set(list(pchids))
    rchids, pchids = set(list(rchids)), set(list(pchids))
    uniqCh = set(chids.values())
    for ch in list(rch) :
        if ch in uniqCh : rch.remove(ch)
        else : print 'chain %s reqd but absent, fatal' % ch
    for ch in list(pch) :
        if ch in uniqCh : pch.remove(ch)
        else : print 'chain %s reqd but absent, fatal' % ch
    print "Reqd chains", rchids, pchids, "present"
    rch, pch = set(list(rchids)), set(list(pchids))

    ## apply missing sc correction
    scMissInds = incompleteSCcorrection(res, resns, pts)
    print "SCMISSINDS", scMissInds
    for i in scMissInds :
        print "SCMISSINDS", resids[i]

    ## find RNA's min-dist from other residues. find contiguous segments entirely within 10A distance of ligand
    leastDists = findRNAdists(rch, res, pts, chids)

    closeRes = [None] * len(leastDists)
    for ind,d in leastDists.items() :
        if chids[ind] in pch and d < closeCutoff : closeRes[ind] = 1
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

    knownPositions = range(len(pts))
    blist, numtrials, bfoll, irlist, optRests, dummies = [], [], {}, [], [], {}
    ## rna chains
    for chid in rch :
        for ind in res.keys() :
            if chids[ind] in rch :
                for pi in res[ind].values() : knownPositions.remove(pi)
        blist1, numtrials1, bfoll1, irlist1, dummies1 = prepareRNAchain(chid, res, resids, resnums, resns, chids, inscodes, pts, rnaRad)
        numbuilders = len(blist)
        for bi in range(len(blist1)) :
            blist.append( blist1[bi] )
            numtrials.append( numtrials1[bi] )
            if not bfoll1.has_key(bi) : continue #no followers for this builder
            assert not bfoll.has_key(numbuilders+bi)
            bfoll[ numbuilders+bi ] = []
            for fi in bfoll1[bi] :
                bfoll[ numbuilders+bi ].append( numbuilders+fi )
        for r in irlist1 : irlist.append(r)
        for k in dummies1.keys() :
            if not dummies.has_key(k) : dummies[k] = []
            for v in dummies1[k] : dummies[k].append(v)
    ## protein loops
    for seg in segs :
        startindex, endindex = seg[0], seg[1]
        print "SEGMENT -------- [%s] [%s]" % (resids[startindex], resids[endindex]), seg
        if chids[startindex] != chids[endindex] :
            print "end-of-chain close-to-ligand regions are not yet supported" ; sys.exit(0)
        for ind in range(startindex-1,endindex+2) :
            for an,pi in res[ind].items() :
                if ind == startindex-1 and an in [' N  ',' CA '] : continue
                if ind == endindex+1 and an in [' CA ',' C  ',' O  '] : continue
                if pi in knownPositions : knownPositions.remove(pi)
        blist1, bfoll1, numtrials1, irlist1 = preparePeptideLoop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, None, caRad, scRad, scMissInds)
        blistsize = len(blist)
        for b in blist1 : blist.append(b)
        for bi in bfoll1.keys() :
            if not bfoll.has_key(bi+blistsize) : bfoll[bi+blistsize] = []
            for fi in bfoll1[bi] : bfoll[bi+blistsize].append(fi+blistsize)
        for n in numtrials1 : numtrials.append(n)
        for r in irlist1 : irlist.append(r)

    print "Residues and atom-numbers-----------------------"
    keys = res.keys() ; keys.sort()
    for k in keys : print "(%d)"%k, resids[k], res[k], chids[k]
    print "------------------------------------------------\n\n\n"

    for b in blist :
        b.describe() ; print ''
    print "------------------------------------------------\n\n\n"

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
    for k in dummies.keys() :
        for an in dummies[k] : aiSC[ res[k][an] ] = -1 # dummy atoms

    if restrGen :
        xrlist, xop = restrGen.generate(blist, aiSC)
        for r in xrlist : irlist.append(r)
        for ri in xop : optRests.append(ri+len(irlist)-len(xrlist))

    ssReduction = PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [], outpdb)

    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(backtrack, popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, 1, irlist, optRests, modelRenderer, nmodels*100, nmodels)
    strategy.execute()

def callmain() :
    from commonOptions import makeParser, parseOptions, addXrayOptions
    parser = makeParser()
    parser = addXrayOptions(parser)
    parser.add_option("--protChains", action='store', type='string', dest='pchids', help='chain ids of protein chains', default='')
    parser.add_option("--rnaChains", action='store', type='string', dest='rchids', help='chain ids of rna chains', default='')
    parser.add_option("--rna-restraint-radius", action='store', type='float', dest='rnaRad', help='restraint sphere radius for rna backbone', default=2)
    parser.add_option("--around-rna", action='store', type='float', dest='closeCutoff', help='how close around rna is close', default=10)
    parser.remove_option("--guided-sampling")
    parser.remove_option("--buildN2C")
    parser.remove_option("--mconly")
    options = parseOptions(parser)

    xrayRestGen = None
    if options.mtzfn :
        from prepareChain import XrayRestraintsGenerator
        xrayRestGen = XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, options.sigXmin, options.sigXmax, options.sigXmean)

    main(options.pdbfile, options.pchids, options.rchids, options.closeCutoff,
        options.rnaRad, options.caRad, options.scRad, options.scReduction,
        options.outpdb, options.nmodels, options.popsize, options.backtrack,
        xrayRestGen)

if __name__ == "__main__" : callmain()
