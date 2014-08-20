from pdbr import protein
import prot2res
from prepareChain import prepareChainTerminal
from builders import VecFloat, VecInt, VecVecFloat, DihedBuilder, OptFragPlacer
from geometry import calcDist, calcAngle, calcDihed, CAtraceGH
from samplers import FragSampler
import data
from data import vdwr

def linkermain() :
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='template PDB file')
    parser.add_option("--startres", action='store', type='string', dest='startres', help='starting residue identifier', default=None)
    parser.add_option("--endres", action='store', type='string', dest='endres', help='ending residue identifier', default=None)
    parser.add_option("--num-models-wanted", action='store', type='int', dest='nmodels', help='number of models desired. number of attempts is generally 10 times this', default=100)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--outpdb", action='store', type='str', dest='outpdb', help='file to write the models to')
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    (options, args) = parser.parse_args()

    import misc
    misc.setVerbosity(options.verbose)

    prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    reskeys = list(res.keys()) ; reskeys.sort()
    startindex, endindex = None, None
    for ri in reskeys :
        if resids[ri] == options.startres : startindex = ri
        elif startindex != None and resids[ri] == options.endres : endindex = ri
    assert startindex and endindex
    print "Found startindex,endindex [%s][%s]----[%s][%s]" % (resids[startindex-1],resids[startindex], resids[endindex],resids[endindex+1])

    knownPositions = []
    for ri in range(0,startindex-1) :
        for name,ai in res[ri].items() : knownPositions.append(ai)
    knownPositions.append( res[startindex-1][' N  '] )
    knownPositions.append( res[startindex-1][' CA '] )

    blist, bfoll, numtrials, rlist = prepareChainTerminal("Cterm", startindex, endindex, endindex+1, res, resids, resnums, resns, chids, inscodes, pts, None, 1e6, 1e6, None, [])

    ips, ops = [], []
    pC = VecFloat(pts[res[endindex][' C  ']]) ; ips.append(res[endindex][' C  '])
    pN = VecFloat(pts[res[endindex+1][' N  ']]) ; ips.append(res[endindex+1][' N  '])
    pCA = VecFloat(pts[res[endindex+1][' CA ']]) ; ips.append(res[endindex+1][' CA '])
    pC1 = VecFloat(pts[res[endindex+1][' C  ']]) ; ops.append(res[endindex+1][' C  '])
    blist.append( DihedBuilder( VecInt(ips), VecInt(ops), "end+1 phi builder", calcDist(VecFloat(pCA),VecFloat(pC1)),
        calcAngle(VecFloat(pN),VecFloat(pCA),VecFloat(pC1)), calcDihed(VecFloat(pC),VecFloat(pN),VecFloat(pCA),VecFloat(pC1)) ) )
    numtrials.append(1)

    ## make a fragments sampler and builder for later domain. frag builder samples the fragment and places the sampled fragment to have optimal overlap between its inputs points and corr.points in fragment
    ## in this case those points are e+1's N,CA,C
    frag, fpis, fpos = [], [], []
    fpis.append(res[endindex+1][' N  ']) ; fpis.append(res[endindex+1][' CA ']) ; fpis.append(res[endindex+1][' C  '])
    for i in fpis : frag.append(list(pts[i]))
    for ri in range(endindex+1,len(reskeys)) :
        for an,ai in res[ri].items() :
            if ai in fpis : continue
            fpos.append(ai) ; frag.append(list(pts[ai]))
    frag = VecVecFloat(frag)

    fs = FragSampler(1, frag)
    blist.append( OptFragPlacer(VecInt(fpis), VecInt(fpos), fs, VecInt(range(3)), VecInt(range(3,len(frag))), "OptFragPlacer") )
    numtrials.append(1)

    radii = [0] * len(pts) ## make radii for all atoms including dummies
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else :
                try : radii[pi] = vdwr['XXX'][aname[1]]
                except KeyError : radii[pi] = vdwr['XXX']['C']

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
    scReduction, ssReduction = 1, data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [], options.outpdb)

    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(options.popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, 1, rlist, [], modelRenderer, options.nmodels*100, options.nmodels)
    strategy.execute()

if __name__ == "__main__" : linkermain()
