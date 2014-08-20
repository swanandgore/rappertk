from pdbr import protein
import prot2res
from prepareChain import  addNCdummyGly, removeSC, mergeBRlists, makeCApropRestraints
from builders import NanchorBuilder, VecFloat, VecVecFloat
import data
from restraints import EnvelopeRestraint
from data import vdwr, consts
from geometry import CAtraceGH, Grid
from builders import VecInt
import sys,os,shutil
from xcheck import XrayScorer, XrayRanker








def bootstrapCAprop(bsA,bsB, caRad, res,resids,pts,ncDummies) :
    bs1, bs2, rlist = bsA, bsB, []
    print bs1, bs2
    assert bs1+1 == bs2
    if res.has_key(bs1-1) and not ncDummies.has_key(bs1-1) : bs1 -= 1
    if res.has_key(bs2+1) and not ncDummies.has_key(bs2+1) : bs2 += 1
    rlist += makeCApropRestraints(bsA, 1, bs1,bsB, caRad, res, resids, pts)
    rlist += makeCApropRestraints(bsB, 1, bsA,bs2, caRad, res, resids, pts)
    return rlist




def parseResids(l) :
    start, resids = None, []
    for i in range(len(l)) :
        if l[i] == '[' : start = i
        elif l[i] == ']' :
            resids.append( l[start+1:i] )
    return resids


def parseResnums(l) :
    start, resnums = None, []
    for i in range(len(l)) :
        if l[i] == '[' : start = i
        elif l[i] == ']' :
            resnums.append( l[start+1:i] )
#print resnums
    #sys.exit()
    return resnums












def parseSSfile(ssfile, resids) :
    helices, sheets = [], []
    resids2index = {}
    for k,v in resids.items() :
        resids2index[v] = k



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
                if helix[1] < helix[0] :
                    helix = (helix[1],helix[0])
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
                print rids[i],
                corr = (resids2index[rids[i]],resids2index[rids[i+1]])
                corrs.append(corr)
            sheet.append(corrs)
            sheets.append(sheet)
    return helices, sheets




def createXEnvelopeRestraintsCA(blist, mapfile, aiSC) :
    from restraints import EDrestraint
    erlist = []
    for b in blist :
        bop = b.getOP()
        edrIP = []
        if "Pept" in b.name() or "anchor" in b.name(): 
            for i in range(bop.size()) :
                print bop[i]
                #if aiSC[bop[i]] != -1 :
                edrIP.append(bop[i])
                rstr = "XEnvelopeRestraint on op of %s" % b.name()
            esmin, esmax, esmean, rcmult = .000, 5., .0, 5
            erlist.append( EDrestraint.makeEDrestraintFromMap(VecInt(edrIP), rstr, mapfile,esmin, esmax, esmean) )
        else :
            print "DDD", b.name() 
            continue

    return erlist



def parseSSfile2(ssfile, resids , resnums) :

    helices, sheets = [], []
    resnums2index = {}
    for k,v in resnums.items() :
        resnums2index[int(v)] = int(k)
    lines = open(ssfile, 'r').readlines()
    import re
    for li in range(len(lines)) :
        l = re.sub("\n", "", lines[li])
        if l == '' or l[0] == '#' : continue
        if l[0:6] == "helix " :
            helixtype = l.split()[1]
            assert helixtype in ["4_13", "3_10"]
            rids = parseResnums(l)
            print rids

            for i in range(0,len(rids),2) :
                helix = (resnums2index[int(rids[i])],resnums2index[int(rids[i+1])])
                if helix[1] < helix[0] : helix = (helix[1],helix[0])
                helices.append([helixtype, helix])

        if l[0:6] == "sheet " :
            rids = parseResnums(l)
            print rids
            strands = []
            for i in range(0,len(rids),2) :
                strand = (resnums2index[int(rids[i])],resnums2index[int(rids[i+1])])
                if strand[1] < strand[0] : strand = (strand[1],strand[0])
                strands.append(strand)
            sheet = [strands]
            li += 1 ; l = lines[li]
            assert l[0:11] == "strand-pair"
            l = re.sub("  *", " ", re.sub("\n", "", l))
            sheet.append(l[12:].split())
            li += 1 ; l = lines[li]
            assert l[0:21] == "strand-correspondence"
            rids = parseResnums(l)
            print rids

            corrs = []
            for i in range(0,len(rids),2) :
                print rids[i],
                corr = (resnums2index[int(rids[i])],resnums2index[int(rids[i+1])])
                corrs.append(corr)
            sheet.append(corrs)
            sheets.append(sheet)
    return helices, sheets






def parseSSfile3(ssfile, resids , resnums,chids) :
    print "chids = " , chids

    helices, sheets = [], []
    resnums2index = {}
    for k,v in resnums.items() :
        nk = str(int(v)) + chids[k]
        resnums2index[nk] = k
        print "roar",k , v , chids[k] , nk

    

    lines = open(ssfile, 'r').readlines()
    import re
    for li in range(len(lines)) :
        l = re.sub("\n", "", lines[li])
        if l == '' or l[0] == '#' : continue
        if l[0:6] == "helix " :
            helixtype = l.split()[1]
            assert helixtype in ["4_13", "3_10"]
            rids = parseResnums(l)


            for i in range(0,len(rids),2) :
                helix = (resnums2index[rids[i]],resnums2index[rids[i+1]])
                if helix[1] < helix[0] :
                    helix = (helix[1],helix[0])
                helices.append([helixtype, helix])

        if l[0:6] == "sheet " :
            rids = parseResnums(l)

            strands = []
            for i in range(0,len(rids),2) :
                strand = (resnums2index[rids[i]],resnums2index[rids[i+1]])
                if strand[1] < strand[0] : strand = (strand[1],strand[0])
                strands.append(strand)
            sheet = [strands]
            li += 1 ; l = lines[li]
            assert l[0:11] == "strand-pair"
            l = re.sub("  *", " ", re.sub("\n", "", l))
            sheet.append(l[12:].split())
            li += 1 ; l = lines[li]
            assert l[0:21] == "strand-correspondence"
            rids = parseResnums(l)
            
            corrs = []
            for i in range(0,len(rids),2) :

                corr = (resnums2index[rids[i]],resnums2index[rids[i+1]])
                corrs.append(corr)
            sheet.append(corrs)
            sheets.append(sheet)
    return helices, sheets










def parseSSfileAK(ssfile,resids,resnums) :
    helices, sheets = [], []
    resids2index = {}
    resnums2index = {}
    for k,v in resids.items() :
        resids2index[v] = k


    for k,v in resnums.items() :
        resnums2index[int(v)] = int(k)


    lines = open(ssfile, 'r').readlines()
    import re
    for li in range(len(lines)) :
        l = re.sub("\n", "", lines[li])
        if l == '' or l[0] == '#' : continue
        if l[0:6] == "helix " :
            helixtype = l.split()[1]
            assert helixtype in ["4_13", "3_10"]
            rids = parseResnums(l)
            for i in range(0,len(rids),2) :
                helix = (resnums2index[int(rids[i])],resnums2index[int(rids[i+1])])
                #if helix[1] < helix[0] :
                #    helix = (helix[1],helix[0])

                helices.append([helixtype, helix])
        if l[0:6] == "sheet " :
            rids = parseResnums(l)
            strands = []
            for i in range(0,len(rids),2) :
                strand = (resnums2index[int(rids[i])],resnums2index[int(rids[i+1])])
                if strand[1] < strand[0] : strand = (strand[1],strand[0])
                strands.append(strand)
            sheet = [strands]
            li += 1 ; l = lines[li]
            assert l[0:11] == "strand-pair"
            l = re.sub("  *", " ", re.sub("\n", "", l))
            sheet.append(l[12:].split())
            li += 1 ; l = lines[li]
            assert l[0:21] == "strand-correspondence"
            rids = parseResnums(l)
            corrs = []
            for i in range(0,len(rids),2) :
                print rids[i],
                corr = (resnums2index[int(rids[i])],
                        resnums2index[int(rids[i+1])])
                corrs.append(corr)
            sheet.append(corrs)
            sheets.append(sheet)
    print sheets
    #import  sys; sys.exit()
    return helices, sheets


def createXEnvelopeRestraints(blist, mapfile, aiSC) :
    from restraints import EDrestraint
    erlist = []
    for b in blist :
        bop = b.getOP()
        edrIP = []
        for i in range(bop.size()) :
            if aiSC[bop[i]] != -1 : edrIP.append(bop[i])
        rstr = "XEnvelopeRestraint on op of %s" % b.name()
        esmin, esmax, esmean, rcmult = .000, 5., .0, 5
        erlist.append( EDrestraint.makeEDrestraintFromMap(VecInt(edrIP), rstr, mapfile,esmin, esmax, esmean) )
    return erlist


def createXEnvelopeRestraints2(blist, mapfile, aiSC,ncDummies) :
    from restraints import EDrestraint
    erlist = []
    
    print "aiSC",aiSC
    print "NCD",ncDummies
    for b in blist :
        bop = b.getOP()
        edrIP = []
        rstr = "XEnvelopeRestraint on op of %s" % b.name()
        print rstr
        for i in range(bop.size()) :
            print i , bop[i]

            if aiSC[bop[i]] != -1 :
                edrIP.append(bop[i])

        esmin, esmax, esmean, rcmult = .000, 5., .0, 5
        erlist.append( EDrestraint.makeEDrestraintFromMap(VecInt(edrIP), rstr, mapfile,esmin, esmax, esmean) )
    return erlist


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


def main() :
    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='pdb containing a model of pdb-ligand complex')
    parser.add_option("--outpdb", action='store', type='str', dest='outpdb', help='file to write the models to')
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--ssfile", action='store', type='string', dest='ssfile', help='list of secondary structures, see applications/ssfile for example')
    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--loopca-restraint-radius", action='store', type='float', dest='loopCaRad', help='Radius for loop residues', default=1)

    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='to create all the files during refinement. it shdnt be already present.')

    
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=500)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=0.75)
    parser.add_option("--num-models-wanted", action='store', type='int', dest='nmodels', help='number of models desired. number of attempts is generally 10 times this', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--mconly", action='store', type='int', dest='mconly', help='mainchain only', default=None)
    parser.add_option("--mapfile", action='store', type='string', dest='mapfile', help='ccp4 map restraining the shape of the molecule', default=None)
    

    parser.add_option("--guided-sampling", action='store', type='int', dest='informed', help='whether to use guided sampling. No by deafult', default=None)

    parser.add_option("--make-optional", action='store', type='int', dest='makeOpt', help='whether to make SS restraints optional . Yes by default :1 / 0: Enforce ', default=1)

    (options, args) = parser.parse_args()

    if options.informed : options.informed = options.caRad

    import misc
    misc.setVerbosity(options.verbose)

    if not os.path.isdir(options.scratchdir) : os.mkdir(options.scratchdir)

    if options.ssfile !=None:
        shutil.copyfile(options.ssfile,"%s/%s" % (options.scratchdir,options.ssfile))
    shutil.copyfile(options.pdbfile, "%s/%s" % (options.scratchdir,options.pdbfile))
    shutil.copyfile(options.mapfile, "%s/%s" % (options.scratchdir,options.mapfile))


    os.chdir(options.scratchdir)

    prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    
    ptID = {}
    for k,v in res.items():
        for a,b in v.items():
            ptID[b] = "%s-%s"%(resids[k],a) 
        
    #import sys; sys.exit()
    assert options.mconly == 1 or options.mconly == None
    if options.mconly : res, pts = removeSC(res, pts, res.keys())

    numres = len(resids)
    numpts = len(pts)
    ncDummies, firstindex, lastindex = addNCdummyGly(chids[0], res, resids, resnums, resns, chids, inscodes, pts)
    assert numres + 2 == len(resids) ; assert numpts < len(pts)

    helices, sheets = parseSSfileAK(options.ssfile, resids,resnums)

    freeResids = {}
    for k in resids.keys() : freeResids[k] = "loop"
    bootstraps = set()

    newsheets = []
    for strands, interStrand, corrInds in sheets :
        dirs = ["fwd"] * len(strands)
        sheetBS = []
        print len(strands) , len(interStrand), strands
        for i in range(len(interStrand)) : ## arbitrarily assign building directions
            if interStrand[i] == "parallel" :
                dirs[i+1] = "fwd"
            else :
                print i,dirs
                dirs[i+1] = "bkwd"
            
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
                freeResids[b1] = "bootstrap"
                freeResids[b2] = "bootstrap"
                sheetBS.append((b1,b2))
            elif freeResids[b1] == freeResids[b2] == "bootstrap" :  sheetBS.append(None)
            else :
                print "FR",freeResids[b1] ,freeResids[b2],resids[b1],resids[b2]
                for k,v in freeResids.items() :
                    print resnums[k],v
                assert None
        newsheets.append( [strands, dirs, interStrand, corrInds, sheetBS] )


        for k,v in freeResids.items() :
            print "ZZ",resnums[k],v

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
    nter = firstindex
    cter = lastindex
    
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

    for i in resids.keys() :
        print "Resids",resids[i], freeResids[i]

    
    blist, numtrials, bfoll, rlist,hydlist = [], [], {}, [], []
    
    from prepareChain import PrepareChain
    pc = PrepareChain("PRL")

    for strands, dirs, interStrand, corrInds, sheetBS in sheets : # sheets
        print "BOOTSTRAP", sheetBS
        starts, ends = [], []
        for beg,end in strands : starts.append(beg) ; ends.append(end)

        blist1, bfoll1, numtrials1, rlist1, hydlistB = pc.prepareBetaSheet(res, resids, resnums, resns, chids, inscodes, pts, options.caRad, starts, ends, sheetBS, dirs, corrInds, options.mconly, options.informed)
        
        
        mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1)
        for h in hydlistB : hydlist.append(h)

        for strandBS in sheetBS :
            if strandBS : rlist += bootstrapCAprop(strandBS[0],strandBS[1], options.caRad,res,resids,pts,ncDummies)


    for htype,helix,helixBS in helices : # helices
        print "BOOTSTRAP", helixBS
        blist1, bfoll1, numtrials1, rlist1,hydlistH = pc.prepareAlphaHelix(htype, res, resids, resnums, resns, chids, inscodes, pts, options.caRad, helix[0], helix[1], helixBS, options.mconly, options.informed)

        mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1)
        for h in hydlistH : hydlist.append(h)

        if helixBS :
            rlist += bootstrapCAprop(helixBS[0],helixBS[1], options.caRad,res,resids,pts,ncDummies)



    for loop in loops : # loops
        blist1, bfoll1, numtrials1, rlist1 = pc.preparePeptideLoop1(loop[0], loop[1], res, resids, resnums, resns, chids, inscodes, pts, options.mconly, options.loopCaRad, 1e10, [])

        
        mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1)


    
    nstart = 0
    nstop = nter 
    
    if nter < 0 :
        nstart = nter 
        nstop =  0
    blist1, bfoll1, numtrials1, rlist1 = pc.prepareChainTerminal("Nterm", nstart,nstop, firstindex, res, resids, resnums, resns, chids, inscodes, pts, options.mconly, 1000, 1000, [])
    mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1) # N term, bkwd build

    
    
    blist1, bfoll1, numtrials1, rlist1 = pc.prepareChainTerminal("Cterm", cter, len(resids)-3, lastindex, res, resids, resnums, resns, chids, inscodes, pts, options.mconly, 1000, 1000, [])
    mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1) # C term, fwd build


    print "-----------------POINTSET------------------"
    for k,v in res.items() :
        for a, b in v.items():
            print resids[k],pts[b]

    print "-----------------BUILDERS------------------"
    for b in blist :
        print "----------------------------------------"
        print b.name()
        bop = b.getOP()
        bip = b.getIP()
        #for x in range(bop.size()):
        #    if bop[x] in ptID.keys():
        #        print "OP",bop[x],ptID[bop[x]]

       # for x in range(bip.size()):
       #     if bip[x] in ptID.keys():
       #         print "IP",bip[x],ptID[bip[x]]
    ##sys.exit()
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
            if bop[i] in knownPositions :
                knownPositions.remove(bop[i])
    for ki in knownPositions : aiSC[ki] = -1 ##XXX

    if options.mapfile : rlist += createXEnvelopeRestraints(blist, options.mapfile, aiSC)

    ## Add a ranker
    
    esmin, esmax, esmean, rcmult = .000, 5., .0, 5
    ranker = XrayRanker(options.mapfile, "map", "FC", "PHIC", "F1", esmin, esmax)
    ranker.rankChildren = rcmult ; ranker.rankRecurse = 1
    ranker.rankLeaderBuilderOnly = None ; ranker.rankGivenEnsemble = None
    


    
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
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, ncDummies.keys(), options.outpdb)

    #for pi in range(len(pts)) :
    #    if pi in knownPositions : continue
    #    pts[pi] = (0.,0.,0.)


    print "known ..",
    for k in knownPositions:
        print k



    
    extraRestraints, extraRestrOpt = [], []
    optional = []
    

    for r in hydlist :
        extraRestraints.append(r)
        optional.append(len(extraRestraints)-1)
        

    rsize = len(rlist)

    for xr in extraRestraints :
        rlist.append(xr)

    for ori in optional :
        extraRestrOpt.append(rsize + ori)

    if options.makeOpt == 0:
        extraRestrOpt = []


    natt = 25
    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(options.backtrack,options.popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll,None,rlist ,extraRestrOpt, modelRenderer, natt, options.nmodels,None)
    strategy.snapRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [], "snap.pdb")
    strategy.ranker = ranker




    strategy.execute()


if __name__ == "__main__" : main()
