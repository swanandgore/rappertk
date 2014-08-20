import math, os, sys, optparse

from prepareChain import removeSC, preparePeptideLoop
from pdbr import protein
import prot2res
from geometry import CAtraceGH
import data
from data import vdwr
from builders import VecVecFloat, VecInt

def main(pdbfile, startResid, endResid, mconly, caRad, scRad, scReduction, guidedSampling, popsize, backtrack, nmodels, outpdb, restrGen=None) :
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)

    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    startindex, endindex = None, None
    for index in range(len(resids)) :
        print "----%s-----%s------" % (resids[index], startResid)
        if resids[index] == startResid : startindex = index
        if resids[index] == endResid : endindex = index

    assert startindex <= endindex and startindex != None and endindex != None
    assert chids[startindex] == chids[endindex] == chids[startindex-2] == chids[endindex+2]

    if mconly : res, pts = removeSC(res, pts, res.keys())

    knownPositions = range(len(pts))
    for ind in range(startindex-1,endindex+2) :
        for an,pi in res[ind].items() :
            knownPositions.remove(pi)
    knownPositions.append(res[startindex-1][' N  '])
    knownPositions.append(res[startindex-1][' CA '])
    knownPositions.append(res[endindex+1][' CA '])
    knownPositions.append(res[endindex+1][' C  '])
    knownPositions.append(res[endindex+1][' O  '])

    guidedRadius = None
    if guidedSampling : guidedRadius = caRad
    blist, bfoll, numtrials, irlist = preparePeptideLoop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, [], guidedRadius)

    print "Residues and atom-numbers-----------------------"
    keys = res.keys() ; keys.sort()
    for k in keys : print resids[k], res[k]
    print "------------------------------------------------\n\n\n"

    for r in irlist : print r.name()
    for b in blist : b.describe() ; print ''
    print "------------------------------------------------\n\n\n"

    radii = [0] * len(pts) ## make radii for all atoms including dummies
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else :
                try : radii[pi] = vdwr['XXX'][aname[1]]
                except KeyError : radii[pi] = vdwr['XXX']['C']

    for index in res.keys() :
        for aname,pi in res[index].items() :
            if pi in knownPositions : continue
            pts[pi] = (0.,0.,0.)
    pts = VecVecFloat(pts)

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

    optRestraints = []
    if restrGen :
        xrlist, optional = restrGen.generate(blist, aiSC)
        irsize = len(irlist)
        for xr in xrlist : irlist.append(xr)
        for ori in optional : optRestraints.append(irsize + ori)

    scReduction, ssReduction = scReduction, data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [], outpdb)

    assert len(radii) == len(pts)

    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(backtrack, popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, 1, irlist, optRestraints, modelRenderer, nmodels*100, nmodels)
    strategy.execute()

def callmain() :
    from commonOptions import makeParser, parseOptions, addXrayOptions
    parser = makeParser()
    parser = addXrayOptions(parser)
    parser.add_option("--startres", action='store', type='string', dest='startResid', help='resid of start loop residue')
    parser.add_option("--endres", action='store', type='string', dest='endResid', help='resid of start loop residue')
    parser.remove_option("--buildN2C")
    options = parseOptions(parser)

    xrayRestGen = None
    if options.mtzfn :
        from prepareChain import XrayRestraintsGenerator
        xrayRestGen = XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, options.sigXmin, options.sigXmax, options.sigXmean)

    main(options.pdbfile, options.startResid, options.endResid, options.mconly,
            options.caRad, options.scRad, options.scReduction, options.guidedSampling,
            options.popsize, options.backtrack, options.nmodels, options.outpdb,
            xrayRestGen)

if __name__ == "__main__" : callmain()
