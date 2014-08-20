from pdbr import protein
import prot2res
from prepareChain import prepareBetaSheet, addNCdummyGly, prepareAlphaHelix, preparePeptideLoop, prepareChainTerminal, removeSC, mergeBRlists, makeCApropRestraints
from builders import NanchorBuilder, VecFloat, VecVecFloat
import data
from restraints import EnvelopeRestraint
from data import vdwr, consts
from geometry import CAtraceGH, Grid
from builders import VecInt

def bootstrapCAprop(bsA,bsB, caRad, res,resids,pts,ncDummies) :
    bs1, bs2, rlist = bsA, bsB, []
    print bs1, bs2
    assert bs1+1 == bs2
    if res.has_key(bs1-1) and not ncDummies.has_key(bs1-1) : bs1 -= 1
    if res.has_key(bs2+1) and not ncDummies.has_key(bs2+1) : bs2 += 1
    rlist += makeCApropRestraints(bsA, 1, bs1,bsB, caRad, res, resids, pts)
    rlist += makeCApropRestraints(bsB, 1, bsA,bs2, caRad, res, resids, pts)
    return rlist

def betamain() :
    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='pdb containing a model of pdb-ligand complex')
    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=0.75)
    parser.add_option("--outpdb", action='store', type='str', dest='outpdb', help='file to write the models to')
    parser.add_option("--num-models-wanted", action='store', type='int', dest='nmodels', help='number of models desired. number of attempts is generally 10 times this', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    (options, args) = parser.parse_args()

    import misc
    misc.setVerbosity(options.verbose)

    prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    starts, ends = [ "THR    2 ", "GLU   49 ", "GLU   37 " ], [ "VAL    8 ", "THR   58 ", "ILE   46 " ]
    corrRes = ["TYR    4 ", "LEU   53 ", "LYS   42 "]
    directions = ["bkwd", "bkwd", "fwd"]
    for index in res.keys() :
        #print "[%s]" % resids[index]
        if resids[index] in starts : starts[ starts.index(resids[index]) ] = index
        if resids[index] in ends : ends[ ends.index(resids[index]) ] = index
        if resids[index] in corrRes : corrRes[ corrRes.index(resids[index]) ] = index
    print starts, ends, corrRes
    cres = [(corrRes[0],corrRes[1]), (corrRes[1],corrRes[2])]

    blist, numtrials, bfoll = [], [], {}
    # make bootstraps before fwd and after bkwd strands
    for si in range(len(starts)) :
        if directions[si] == "fwd" : ind = starts[si]-2
        else : ind = ends[si]+1
        IPs, OPs = [], [ res[ind][' CA '], res[ind][' C  '], res[ind][' O  '], res[ind+1][' N  '], res[ind+1][' CA '], ]
        blist.append( NanchorBuilder( VecInt(IPs), VecInt(OPs), consts, "NanchorBuilder [%s][%s]"%(resids[ind],resids[ind+1]),
            VecFloat(pts[res[ind][' CA ']]), options.caRad, VecFloat(pts[res[ind+1][' CA ']]), options.caRad ) )
        numtrials.append(100)

    blist1, bfoll1, numtrials1, irlist = prepareBetaSheet(res, resids, resnums, resns, chids, inscodes, pts, options.caRad, starts, ends, directions, cres, mconly, options.caRad)
    blistsize = len(blist)
    for bi in range(len(blist1)) :
        blist.append( blist1[bi] )
        numtrials.append( numtrials1[bi] )
        if bfoll1.has_key(bi) :
            bfoll[blistsize+bi] = []
            for bi1 in bfoll1[bi] : bfoll[blistsize + bi].append( blistsize + bi1 )

    print "-----------------BUILDERS------------------"
    for b in blist : print b.name()
    print "-----------------RESTRAINTS----------------"
    for r in irlist : print r.name()
    print "-------------------------------------------"

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

    ## everything that is not built is known and dummy, just for the time being
    knownPositions = range(len(pts))
    for b in blist :
        bop = b.getOP()
        for i in range(bop.size()) : knownPositions.remove(bop[i])
    for ki in knownPositions : aiSC[ki] = -1 ##XXX

    radii = [0] * len(pts) ## make radii for all atoms including dummies
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else :
                try : radii[pi] = vdwr['XXX'][aname[1]]
                except KeyError : radii[pi] = vdwr['XXX']['C']

    ssReduction = data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), options.scReduction, ssReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, options.outpdb)

    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(None, options.popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, None, irlist, [], modelRenderer, options.nmodels*100, options.nmodels)
    strategy.execute()

def parseResids(l) :
    start, resids = None, []
    for i in range(len(l)) :
        if l[i] == '[' : start = i
        elif l[i] == ']' :
            resids.append( l[start+1:i] )
    return resids


def parseSSfile(ssfile, resids) :
    helices, sheets = [], []
    resids2index = {}
    for k,v in resids.items() : resids2index[v] = k
    lines = open(ssfile, 'r').readlines()
    import re
    for li in range(len(lines)) :
        l = re.sub("\n", "", lines[li])
        if l == '' or l[0] == '#' : continue
        if l[0:6] == "helix " :
            helixtype = l.split()[1]
            assert helixtype in ["4_13", "3_10"]
            rids = parseResids(l)
            for i in range(0,len(rids),2) :
                helix = (resids2index[rids[i]],resids2index[rids[i+1]])
                if helix[1] < helix[0] : helix = (helix[1],helix[0])
                helices.append([helixtype, helix])
        if l[0:6] == "sheet " :
            rids = parseResids(l)
            strands = []
            for i in range(0,len(rids),2) :
                strand = (resids2index[rids[i]],resids2index[rids[i+1]])
                if strand[1] < strand[0] : strand = (strand[1],strand[0])
                strands.append(strand)
            sheet = [strands]
            li += 1 ; l = lines[li]
            assert l[0:11] == "strand-pair"
            l = re.sub("  *", " ", re.sub("\n", "", l))
            sheet.append(l[12:].split())
            li += 1 ; l = lines[li]
            assert l[0:21] == "strand-correspondence"
            rids = parseResids(l)
            corrs = []
            for i in range(0,len(rids),2) :
                corr = (resids2index[rids[i]],resids2index[rids[i+1]])
                corrs.append(corr)
            sheet.append(corrs)
            sheets.append(sheet)
    return helices, sheets

def createEnvelopeRestraints(blist, envgrid, envrad) :
    erlist = []
    for b in blist :
        bop = b.getOP()
        rstr = "EnvelopeRestraint on op of %s" % b.name()
        erlist.append( EnvelopeRestraint(bop, rstr, envgrid, envrad) )
    return erlist

## main routine for doing building in EM scenario
## secondary structure CA restraint radii are 3A. no other positional restraints.
## 3A spheres on all atoms simulate the EM envelope
def main(pdbfile, mconly, ssfile, caRad, scReduction, guidedSampling, popsize, outpdb, nmodels, restrGen) :
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    assert mconly == 1 or mconly == None
    if mconly : res, pts = removeSC(res, pts, res.keys())

    numres = len(resids)
    numpts = len(pts)
    ncDummies, firstindex, lastindex = addNCdummyGly(chids[0], res, resids, resnums, resns, chids, inscodes, pts)
    assert numres + 2 == len(resids) ; assert numpts < len(pts)

    helices, sheets = parseSSfile(ssfile, resids)

    freeResids = {}
    for k in resids.keys() : freeResids[k] = "loop"
    bootstraps = set()

    newsheets = []
    for strands, interStrand, corrInds in sheets :
        dirs = ["fwd"] * len(strands)
        sheetBS = []
        for i in range(len(interStrand)) : ## arbitrarily assign building directions
            if interStrand[i] == "parallel" : dirs[i+1] = "fwd"
            else : dirs[i+1] = "bkwd"
        print "SHEET-----------------------------------------------------------------------------------------------"
        for stri in range(len(strands)) :
            strand = strands[stri]
            print "strand   [%s][%s]" % (resids[strand[0]], resids[strand[1]]), dirs[stri]
            for i in range(strand[0], strand[1]+1) : freeResids[i] = "strand"
            b1, b2 = strand[0]-2, strand[0]-1
            if dirs[stri] == "bkwd" : b1, b2 = strand[1]+1, strand[1]+2
            if ncDummies.has_key(b1) or ncDummies.has_key(b2) : ## cant bootstrap on dummy
                print "Cant bootstrap on dummies [%s] [%s]" % (resids[b1], resids[b2]) ; assert None
            if freeResids[b1] == freeResids[b2] == "loop" :
                freeResids[b1] = "bootstrap" ; freeResids[b2] = "bootstrap" ; sheetBS.append((b1,b2))
            elif freeResids[b1] == freeResids[b2] == "bootstrap" :  sheetBS.append(None)
            else : assert None
        newsheets.append( [strands, dirs, interStrand, corrInds, sheetBS] )
        for i in range(len(interStrand)) :
            print "%8s [%s][%s]" % (interStrand[i], resids[corrInds[i][0]], resids[corrInds[i][1]])
        print "----------------------------------------------------------------------------------------------------"
    sheets = newsheets
    newhelices = []
    for htype,helix in helices :
        print "HELIX-----------------------------------------------------------------------------------------------"
        print htype, "[%s][%s]" % (resids[helix[0]], resids[helix[1]])
        for i in range(helix[0], helix[1]+1) : freeResids[i] = "helix"
        assert not(freeResids[helix[0]-2] == freeResids[helix[0]-1] == freeResids[helix[1]+1] == freeResids[helix[1]+2] == "bootstrap")
        bsdir, bootstrapReqd, helixBS = None, None, None
        ## if already bootstrapped, use it; else try making Nter as bootstrap; else try Cterm bootstrap
        if freeResids[helix[0]-2] == freeResids[helix[0]-1] == "bootstrap" : bsdir = "fwd"
        elif freeResids[helix[1]+1] == freeResids[helix[1]+2] == "bootstrap" : bsdir = "bkwd"
        elif not ncDummies.has_key(helix[0]-2) and not ncDummies.has_key(helix[0]-1) and freeResids[helix[0]-2] == freeResids[helix[0]-1] == "loop" :
            freeResids[helix[0]-2], freeResids[helix[0]-1], bsdir, bootstrapReqd = "bootstrap", "bootstrap", "fwd", 1
        elif not ncDummies.has_key(helix[1]+1) and not ncDummies.has_key(helix[1]+2) and freeResids[helix[1]+1] == freeResids[helix[1]+2] == "loop" :
            freeResids[helix[1]+1], freeResids[helix[1]+2], bsdir, bootstrapReqd = "bootstrap", "bootstrap", "bkwd", 1
        else : assert None
        if bsdir == "fwd" :
            newhelix = (helix[0],helix[1])
            if bootstrapReqd : helixBS = (helix[0]-2,helix[0]-1)
        else :
            newhelix = (helix[1],helix[0])
            if bootstrapReqd : helixBS = (helix[1]+1,helix[1]+2)
        newhelices.append([htype, newhelix, helixBS])
        print "----------------------------------------------------------------------------------------------------"
    helices = newhelices

    freeResids[firstindex] = "Nterm"
    freeResids[lastindex] = "Cterm"
    for i in range(0, len(freeResids)-2) :
        if freeResids[i] == "loop" : freeResids[i] = "Nterm"
        else : nter = i-1 ; break
    for i in range(len(freeResids)-3, -1, -1) :
        if freeResids[i] == "loop" : freeResids[i] = "Cterm"
        else : cter = i+1 ; break
    loops, loopstart = [], None
    for k in range(len(freeResids)-2) :
        if freeResids[k] == "loop" and not loopstart : loopstart = k
        elif loopstart and freeResids[k] != "loop" : loops.append([loopstart,k-1]) ; loopstart = None

    for i in resids.keys() : print resids[i], freeResids[i]
    blist, numtrials, bfoll, rlist = [], [], {}, []

    for strands, dirs, interStrand, corrInds, sheetBS in sheets : # sheets
        print "BOOTSTRAP", sheetBS
        starts, ends = [], []
        for beg,end in strands : starts.append(beg) ; ends.append(end)
        blist1, bfoll1, numtrials1, rlist1 = prepareBetaSheet(res, resids, resnums, resns, chids, inscodes, pts, caRad, starts, ends, sheetBS, dirs, corrInds, mconly, guidedSampling)
        mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1)
        for strandBS in sheetBS :
            if strandBS : rlist += bootstrapCAprop(strandBS[0],strandBS[1], caRad,res,resids,pts,ncDummies)
    for htype,helix,helixBS in helices : # helices
        print "BOOTSTRAP", helixBS
        blist1, bfoll1, numtrials1, rlist1 = prepareAlphaHelix(htype, res, resids, resnums, resns, chids, inscodes, pts, caRad, helix[0], helix[1], helixBS, mconly, guidedSampling)
        mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1)
        if helixBS : rlist += bootstrapCAprop(helixBS[0],helixBS[1], caRad,res,resids,pts,ncDummies)
    for loop in loops : # loops
        blist1, bfoll1, numtrials1, rlist1 = preparePeptideLoop(loop[0], loop[1], res, resids, resnums, resns, chids, inscodes, pts, mconly, 1e10, 1e10, [])
        mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1)
    blist1, bfoll1, numtrials1, rlist1 = prepareChainTerminal("Nterm", nter, 0, firstindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, 1000, 1000, [])
    mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1) # N term, bkwd build
    blist1, bfoll1, numtrials1, rlist1 = prepareChainTerminal("Cterm", cter, len(resids)-3, lastindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, 1000, 1000, [])
    mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1) # C term, fwd build

    #envpts = VecVecFloat(pts)
    #envelopeGrid = Grid(3, envpts, VecFloat([1.]*len(pts)), None)
    #envelopeGrid.justAdd(VecInt(range(len(pts))))
    #rlist += createEnvelopeRestraints(blist, envelopeGrid, options.envelopeRadius)

    print "-----------------POINTSET------------------"
    for k,v in res.items() : print k,v
    print "-----------------BUILDERS------------------"
    for b in blist : print b.name()
    print "-----------------RESTRAINTS----------------"
    for r in rlist : print r.name()
    print "-------------------------------------------"

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

    ## everything that is not built is known and dummy, just for the time being
    knownPositions = range(len(pts))
    for b in blist :
        bop = b.getOP()
        for i in range(bop.size()) :
            if bop[i] in knownPositions : knownPositions.remove(bop[i])
    for ki in knownPositions : aiSC[ki] = -1 ##XXX

    if restrGen : rlist += restrGen.generate(blist, aiSC)

    radii = [0] * len(pts) ## make radii for all atoms including dummies
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else :
                try : radii[pi] = vdwr['XXX'][aname[1]]
                except KeyError : radii[pi] = vdwr['XXX']['C']

    ssReduction = data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, ncDummies.keys(), outpdb)

    for pi in range(len(pts)) :
        if pi in knownPositions : continue
        pts[pi] = (0.,0.,0.)

    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(None, popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, None, rlist, [], modelRenderer, nmodels*100, nmodels)
    strategy.execute()

def callmain() :
    from commonOptions import makeParser, parseOptions, addXrayOptions
    parser = makeParser()
    parser.add_option("--ssfile", action='store', type='string', dest='ssfile', help='list of secondary structures, see applications/ssfile for example')
    parser.add_option("--mapfile", action='store', type='string', dest='mapfile', help='ccp4 map restraining the shape of the molecule', default=None)
    parser.add_option("--envelope-radius", action='store', type='float', dest='envelopeRadius', help='radius of spheres centered on given coordinates, union of which defines the ebvelope', default=5)
    parser.remove_option("--backtrack")
    parser = addXrayOptions(parser)
    options = parseOptions(parser)

    xrayRestGen = None
    if options.mapfile :
        from prepareChain import XrayRestraintsGenerator
        xrayRestGen = XrayRestraintsGenerator(options.mapfile, "map", "dontcare", "dontcare", "dontcare", options.sigXmin, options.sigXmax, options.sigXmean)

    main(options.pdbfile, options.mconly, options.ssfile, options.caRad, options.scReduction, options.guidedSampling, options.popsize, options.outpdb, options.nmodels, xrayRestGen)

if __name__ == "__main__" : callmain()
