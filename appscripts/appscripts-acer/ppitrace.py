import math, os, sys, string

from builders import VecInt, VecFloat, VecVecFloat
from geometry import MapIntIntFloat, CAtraceGH, calcDist
from pdbr import protein, line2resid, line2resn, line2crd, line2atomname, line2resnum, line2resic, line2chid, makeResid
import prot2res
import data
from restraints import EDrestraint

cbDatapath = os.environ["RTKROOT"] + "/data/"

from data import vdwr, resAtoms, consts, mcConn, scConn
from cbutils import findBuilderOrder, findBuilderRestraintOrder, Build

from prepareChain import removeSC, mergeBRlists, prepareSConly

def findStartEnd(chids, chid) :
    indices = chids.keys()
    indices.sort()
    start,end = None,None
    for k in indices :
        if chids[k] == chid and not start : start = k
        elif chids[k] == chid and start : end = k
    return start,end

def findCloseCA(start1,end1, start2,end2,cutoff, res,pts) :
    close = [None] * (end1+1)
    for ind1 in range(start1,end1+1) :
        for ind2 in range(start2,end2+1) :
            if calcDist( VecFloat(pts[res[ind1][' CA ']]), VecFloat(pts[res[ind2][' CA ']]) ) > 10 : continue
            for an1 in res[ind1].keys() :
                for an2 in res[ind2].keys() :
                    if calcDist( VecFloat(pts[res[ind1][an1]]), VecFloat(pts[res[ind2][an2]]) ) > cutoff : continue
                    close[ind1] = 1 ; break
                if close[ind1] : break
            if close[ind1] : break
    return close

def makeCloseSegs(closeRes) :
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
    return segs

def main(pdbfile, chid12, mconly, sconly, closeCut, caRad, scRad, scReduction, restrGen, outpdb, popsize, backtrack, nmodels) :
    assert len(chid12) == 2
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    if mconly == 1 : res, pts = removeSC(res, pts, res.keys())

    print "Residues and atom-numbers-----------------------"
    keys = res.keys() ; keys.sort()
    for k in keys : print resids[k], res[k]
    print "------------------------------------------------\n\n\n"
    
    start1, end1 = findStartEnd(chids, chid12[0])
    start2, end2 = findStartEnd(chids, chid12[1])
    assert start1 and end1 and start2 and end2

    # which interface loops need rebuilding?
    close1 = findCloseCA(start1,end1, start2,end2,closeCut, res,pts)
    close2 = findCloseCA(start2,end2, start1,end1,closeCut, res,pts)

    segs = makeCloseSegs(close1)
    segs += makeCloseSegs(close2)

    blist, numtrials, bfoll, rlist = [], [], {}, []
    for startindex,endindex in segs :
        print "SEGMENT -------- [%s] [%s]" % (resids[startindex], resids[endindex])
        if chids[startindex] != chids[endindex] :
            print "end-of-chain close-to-ligand regions are not yet supported" ; sys.exit(0)
        if sconly == 1 :
            blist1, bfoll1, numtrials1, rlist1 = prepareSConly(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, scRad, [])
        else :  
            blist1, bfoll1, numtrials1, rlist1 = preparePeptideLoop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, [])
        mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1)

    print "BUILDERS------------------------------"
    for b in blist :
        b.describe() ; print ''
    print "---------------------------------------\n\n\n"

    knownPositions = range(len(pts)) ## make known positions by removing built points
    for b in blist :
        bop = b.getOP()
        for i in range(bop.size()) : knownPositions.remove(bop[i])

    radii = [0] * len(pts) ## make radii
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
            elif name == ' CD ' and resns[index] == 'PRO' : aiSC[ai] = 6 # for disulphide vdw reduction
            else : aiSC[ai] = 4

    optRestraints = []
    if restrGen :
        xrlist, optional = restrGen.generate(blist, aiSC)
        irsize = len(rlist)
        for xr in xrlist : rlist.append(xr)
        for ori in optional : optRestraints.append(irsize + ori)

    scReduction, ssReduction = scReduction, data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [], outpdb)

    overlapReductions = MapIntIntFloat()
    assert len(radii) == len(pts)

    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(backtrack, popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, 1, rlist, optRestraints, modelRenderer, nmodels*100, nmodels)
    strategy.execute()

def callmain() :
    from commonOptions import makeParser, parseOptions, addXrayOptions
    parser = makeParser()
    parser = addXrayOptions(parser)
    parser.add_option("--chids", action='store', type='string', dest='chids', help='chains to be traced', default=' ')
    parser.add_option("--sconly", action='store', type='int', dest='sconly', help='modify sidechains only', default=None)
    parser.add_option("--close-cutoff", action='store', type='float', dest='closeCut', help='how close is close for CAs across interface', default=5)
    parser.remove_option("--buildN2C")
    options = parseOptions(parser)

    xrayRestGen = None
    if options.mtzfn :
        from prepareChain import XrayRestraintsGenerator
        xrayRestGen = XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, options.sigXmin, options.sigXmax, options.sigXmean)

    assert not (options.mconly == 1 and options.sconly == 1)
    main(options.pdbfile, options.chids, options.mconly, options.sconly, options.closeCut,
        options.caRad, options.scRad, options.scReduction, xrayRestGen,
        options.outpdb, options.popsize, options.backtrack, options.nmodels)

if __name__ == "__main__" : callmain()
