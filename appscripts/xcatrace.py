import math, os, sys, string
from restraints import EDrestraint
from builders import VecInt, VecVecFloat
from pdbr import protein
import prot2res
from prepareChain import preparePeptideChain, removeSC, mergeBRlists, printResAtoms
import data, misc
from data import vdwr
from geometry import CAtraceGH

def main(pdbfile, traceChids, caRad, scRad, scReduction, mconly, outpdb, popsize, nmodels,
            restrGen=None, backtrack=None, buildfwd=1, guidedSamplingRadius=None) :
    assert mconly == 1 or mconly == None

    ## for a chain in pdbfile, rebuild using some CA-sphere restraint radius
    ## assume all else in pdbfile as just steric obstructions
    assert len(traceChids) > 0
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    knownPositions, unknownIndices = [], [] ## make known positions
    for k in chids.keys() :
        if chids[k] in traceChids : unknownIndices.append(k)
    for index, val in res.items() :
        for aname, pi in val.items() :
            if index not in unknownIndices : knownPositions.append(pi)

    if mconly == 1 : res, pts = removeSC(res, pts, unknownIndices)

    import checkProtChains
    mcmiss, scmiss, chainBreaks = checkProtChains.check(res, resids, resnums, resns, chids, inscodes, pts, 5)
    newTraceChids, newpossiblechids = [], [ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ' ',]
    for ch in set(chids.values()) : newpossiblechids.remove(ch)
    reskeys = res.keys() ; reskeys.sort()
    start = None
    print chainBreaks
    for chid in traceChids : # change chid of each continuous fragment in chains to be traced
        start,end = None,None
        for k in reskeys :
            if chids[k] == chid and start == None : start = k
            if start != None :
                if k in chainBreaks : end = k+1
                elif reskeys[len(reskeys)-1] == k : end = k+1
                elif chids[k+1] != chid : end = k+1
            if start != None and end != None :
                newch = newpossiblechids[0] ; newTraceChids.append(newch) ; newpossiblechids.remove(newch)
                print "renaming", start, k+1, newch
                for ri in range(start,k+1) : chids[ri] = newch
                start, end = None, None
    for k in reskeys :
        print chids[k], "-------", resids[k]

    
    blist, numtrials, bfoll, rlist, dummies = [], [], {}, [], {}
    for chid in newTraceChids :
        blist1, numtrials1, bfoll1, rlist1, dummies1 = preparePeptideChain(chid, res, resids, resnums, resns, chids, inscodes, pts, mconly, buildfwd, guidedSamplingRadius, caRad, scRad)
        mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1)
        for k in dummies1.keys() : dummies[k] = dummies1[k]

    radii = [0] * len(pts) ## make radii for all atoms including dummies
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else :
                try : radii[pi] = vdwr['XXX'][aname[1]]
                except KeyError : radii[pi] = vdwr['XXX']['C']

    ##set all unknown positions to (0,0,0) just to make sure that we'r not copying anything from given structure
    for index in res.keys() :
        for aname,pi in res[index].items() :
            if pi in knownPositions : continue
            pts[pi] = (0.,0.,0.)
    pts = VecVecFloat(pts)


    if misc.verbose(6) : printResAtoms(res, resids)

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

    optRestraints = []
    if restrGen :
        xrlist, optional = restrGen.generate(blist, aiSC)
        irsize = len(rlist)
        for xr in xrlist : rlist.append(xr)
        for ori in optional : optRestraints.append(irsize + ori)

    scReduction, ssReduction = scReduction, data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, dummies.keys(), outpdb)

    assert len(radii) == len(pts)

    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(backtrack, popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, 1, rlist, optRestraints, modelRenderer, nmodels*100, nmodels)
    strategy.execute()

def callmain() :
    from commonOptions import makeParser, parseOptions, addXrayOptions
    parser = makeParser()
    parser = addXrayOptions(parser)
    parser.add_option("--chids", action='store', type='string', dest='chids', help='chains to be traced', default=' ')
    parser.remove_option("--buildN2C")
    options = parseOptions(parser)

    xrayRestGen = None
    if options.mtzfn :
        from prepareChain import XrayRestraintsGenerator
        xrayRestGen = XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, options.sigXmin, options.sigXmax, options.sigXmean)

    main(options.pdbfile, options.chids, options.caRad, options.scRad, options.scReduction,
        options.mconly, options.outpdb, options.popsize, options.nmodels,
        xrayRestGen, backtrack=options.backtrack, buildfwd = 1, guidedSamplingRadius = options.guidedSampling)

if __name__ == "__main__" : callmain()
