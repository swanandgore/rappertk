from pdbr import protein, makeResid
import prot2res
from builders import VecVecFloat, VecInt
import data
from geometry import CAtraceGH
import prepareChain
from data import vdwr
from misc import verbose

def covConnect(covconn, res, indi, namei, indj, namej) :
    i = res[indi][namei]
    j = res[indj][namej]
    if not i in covconn.keys() : covconn[i] = set()
    if not j in covconn.keys() : covconn[j] = set()
    covconn[i].add(j)
    covconn[j].add(i)


def main(pdbfile, chid, caRad, scRad, scReduction, mconly, outpdb, popsize, nmodels,
            restraintAdditionFunction=None, backtrack=None, buildfwd=1, guidedSamplingRadius=None) :
    assert mconly == 1 or mconly == None

    ## for a chain in pdbfile, rebuild using some CA-sphere restraint radius
    ## assume all else in pdbfile as just steric obstructions
    assert len(chid) == 1
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    knownPositions, unknownIndices = [], [] ## make known positions
    for k in chids.keys() :
        if chids[k] == chid : unknownIndices.append(k)
    for index, val in res.items() :
        for aname, pi in val.items() :
            if index not in unknownIndices : knownPositions.append(pi)

    if mconly == 1 : res, pts = prepareChain.removeSC(res, pts, unknownIndices)

    blist, numtrials, bfoll, irlist, dummies = prepareChain.preparePeptideChain(chid, res, resids, resnums, resns, chids, inscodes, pts, mconly, buildfwd, guidedSamplingRadius, caRad, scRad)
    optRestraints = []

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

    if verbose(6) : prepareChain.printResAtoms(res, resids)
    keys = res.keys() ; keys.sort()

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

    if restraintAdditionFunction :
        xrlist, optional = restraintAdditionFunction(blist, aiSC)
        irsize = len(irlist)
        for xr in xrlist : irlist.append(xr)
        for ori in optional : optRestraints.append(irsize + ori)

    scReduction, ssReduction = scReduction, data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, dummies.keys(), outpdb)

    assert len(radii) == len(pts)

    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(backtrack, popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, 1, irlist, optRestraints, modelRenderer, nmodels*100, nmodels)
    strategy.execute()

def callmain() :
    from commonOptions import makeParser, parseOptions
    parser = makeParser()
    parser.add_option("--chid", action='store', type='string', dest='chid', help='chain identifier. DEFAULT blank', default=' ')
    options = parseOptions(parser)

    main(options.pdbfile, options.chid, options.caRad, options.scRad, options.scReduction, options.mconly, options.outpdb, options.popsize, options.nmodels,
        guidedSamplingRadius = options.guidedSampling, backtrack = options.backtrack, buildfwd=options.buildN2C)

if __name__ == "__main__" : callmain()
