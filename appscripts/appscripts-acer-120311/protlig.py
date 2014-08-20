from pdbr import protein
import prot2res
import re, string, sys, os, shutil
from geometry import calcDist, CAtraceGH, MapIntIntFloat
from builders import VecInt, VecFloat, VecBuilder, VecVecFloat, Rotator, TransRotator, BuilderGroup
from restraints import DistanceRestraint
from data import vdwr, consts, PROBE_DISULFIDE_OVERLAP_MARGIN
from prepareChain import  incompleteSCcorrection
import data
from data import vdwr, consts
from data import sgtable , long2shortHM
import random

def findLigDists(ligname, res, pts, resns,hetind) :
    ligdists = {}
    ligind = None
    print "ligname",ligname,len(ligname)
#    import sys ; sys.exit()
    print res.keys()
    for ind in res.keys() :
        if resns[ind] == ligname : ligind = ind ; break
        print resns[ind]
    if ligind == None:
        print "%s not found in ip PDB file"%ligname
        sys.exit()
    for ind in res.keys() :
        if resns[ind] == ligname or 'HOH' in resns[ind] or hetind[ind] == 1:
            print "Continuing for" , ligname
            continue
        for name,ai in res[ind].items() :
            mindist = 1e10
            for lan,lai in res[ligind].items() :
                ldist = calcDist( VecFloat(pts[lai]), VecFloat(pts[ai]) )
                #print "ldist", lan, ldist, pts[lai], pts[ai]
                if ldist < mindist : mindist = ldist
            if mindist < 1e9 :
                print "indie",ind,resns[ind]
                ligdists[ind] = mindist
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
        if 'mindist' in l:
            for fi in range(1,len(flds)) :
                if flds[fi][0] == '[' : break
                tol = string.atof(flds[fi])
                restraints.append( DistanceRestraint(VecInt([res[ligind][anames[0]], res[ligind][anames[1]]]), "DistanceRestraint within ligand", tol, 1e10) )
        elif  'ligname' in l :
            for k,v in resns.items() :
                if v == flds[1] : ligname = flds[1] ; ligind = k ; break
            if ligind == None:
                print "ligname entry not found in ",ligfile
                import sys; sys.exit()
            assert ligind != None
        elif  'rotbond'  in l :
            min, max, step = string.atof(flds[1]), string.atof(flds[2]), string.atof(flds[3]),
            IPs, OPs = [ res[ligind][anames[0]],res[ligind][anames[1]] ], []
            for an in anames[2:] : OPs.append(res[ligind][an])
            builders.append( Rotator(VecInt(IPs), VecInt(OPs), None, "Rotator", min, max, step) )
        elif 'init' in l :
            tol1, tol2 = string.atof(flds[1]), string.atof(flds[2])
            init1, init2 = anames[0], anames[1]
        else :
            print "Unrecognized file format for %s, exiting....." % ligfile ;
            import sys; sys.exit()
            assert None
    assert tol1
    IPs, OPs = [], [res[ligind][init1], res[ligind][init2]]
    allAnames.remove(init1) ; allAnames.remove(init2)
    for an in allAnames : OPs.append(res[ligind][an])
    builders = [ TransRotator(VecInt(IPs),VecInt(OPs),None,"TransRotator", VecFloat(pts[res[ligind][init1]]),tol1, VecFloat(pts[res[ligind][init2]]),tol2) ] + builders
    return builders, restraints, ligname

def main(pdbfile, ligfile, closeCutoff, caRad, scRad, scReduction, guidedSampling, outpdb, popsize, backtrack, nmodels, restrGen, mapfile,ranker,sclib,edopt,allopt) :

    prot = protein(pdbfile, read_hydrogens=0, read_waters=1, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts, hetind = prot2res.readProtRes3(prot)

    ## find ligand name
    ligname = ""
    for l in open(ligfile,'r').readlines() :
        if l[0:7] == "ligname" :
            for xl in l[8:]:
                if xl in [" ","\t","\n"] :
                    continue
                else:

                    ligname = ligname + xl
            
    assert ligname


    ## apply missing sc correction
    scMissInds = incompleteSCcorrection(res, resns, pts)
    for i in scMissInds :
        print "SCMISSINDS", resids[i]

    ## find ligand's min-dist from other residues. find contiguous segments entirely within 10A distance of ligand
    leastDists = findLigDists(ligname, res, pts, resns,hetind)
    print leastDists , len(leastDists.keys())
    closeRes = [None] * (len(leastDists))
    print leastDists , len(leastDists.keys()), len(closeRes)
    for ind,d in leastDists.items() :
        print "indi",ind,d
        if d < closeCutoff :
            print ind
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

    #segs = newsegs
    badresids = []
    for seg in segs:
        for i in range(seg[0],seg[1]+1) :
            badresids.append(resids[i])
            
    from loopbuildV5 import locateRegionsRandomize2
    loops, looptypes = locateRegionsRandomize2(resids, chids, badresids,hetind) # locate loops and order them randomly
    for x in range (len(loops)):
        print "nl", resids[loops[x][0]],resids[loops[x][1]],looptypes[x]

    blist, numtrials, bfoll = [], [], {}
    blist, numtrials, bfoll, rlist,hydlist = [], [], {}, [], []
    dummies = {}
        
    bldrs, irlist, ligname = makeLigBuilder(res, resns, pts, ligfile)
    blist = [ BuilderGroup(VecBuilder(bldrs), "LigandBuilder") ]
    numtrials.append(1000)
    

    ## DANGER make ligand known
    #for ind in res.keys() :
    #    if resns[ind] != ligname : continue
    #    for an,pi in res[ind].items() : knownPositions.remove(pi)
    

    ## for each close segment, setup loop sampling
    #for seg in segs :
    #    startindex, endindex = seg[0], seg[1]
     #   if resns[startindex] == "HOH" : continue
    #    print "SEGMENT -------- [%s] [%s]" % (resids[startindex], resids[endindex]), seg
    #    if chids[startindex] != chids[endindex] :
    #        print "end-of-chain close-to-ligand regions are not yet supported" ; sys.exit(0)
    #    for ind in range(startindex-1,endindex+2) :
    #        for an,pi in res[ind].items() :
    #            if ind == startindex-1 and an in [' N  ',' CA '] : continue
    #            if ind == endindex+1 and an in [' CA ',' C  ',' O  '] : continue
    #            if pi in knownPositions : knownPositions.remove(pi)#

#        guidedRad = None
#        if guidedSampling != None : guidedRad = caRad
#        from prepareChain import PrepareChain
        
 #       pc  = PrepareChain("PRL")
#        
  #      blist1, bfoll1, numtrials1, irlist1 = pc.preparePeptideLoop1(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, None, caRad, scRad, scMissInds, guidedRad)
   #     blistsize = len(blist)
   #     for b in blist1 : blist.append(b)
   #     for bi in bfoll1.keys() :
   #         if not bfoll.has_key(bi+blistsize) : bfoll[bi+blistsize] = []
   #         for fi in bfoll1[bi] : bfoll[bi+blistsize].append(fi+blistsize)
    #    for n in numtrials1 : numtrials.append(n)
        ## DANGER
        #for r in irlist1 : irlist.append(r)

    ## WORKS !
    #XUcaRad = [] ; loop = [] ; looptypes = [] ; midloop = [] ;  badresids = []
    #import prepareChain
    #from loopbuild import Multiloop

    #for s in segs:
    #    XUcaRad.append(caRad)
    #    loop.append(s)
    #    looptypes.append("loop")
    #    midloop.append(None)
    #badresids.append(loop)
    #badresids.append(looptypes)
    #badresids.append(midloop)
    
    #multiPrepC = prepareChain.PrepareChain("PRL")
    #ml = Multiloop(pdbfile, badresids , None, caRad ,scRad, scReduction, None, popsize, backtrack, 1, outpdb , None , multiPrepC, 0 , XUcaRad)
    #nb,notbuilt = ml.run('CA',[],0)
    ## WORKS

    from prepareChainV3 import PrepareChain
    from prepareChain import  addNCdummyGly, removeSC, mergeBRlists, makeCApropRestraints
    prepC = PrepareChain(sclib)

        
    for li in range(len(loops)) :
        startindex, endindex = loops[li]
        print "Start , stop " , startindex , endindex
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
                looptypes[li] = "loop"
                reorderBuilders = None
        else :
            print "Dont understand looptype", looptypes[li] ; assert None

        mergeBRlists(blist, bfoll, numtrials, rlist, bl, bf, nt, rl, dummies, dum)








       # blist1, bfoll1, numtrials1, rlist1 = pc.preparePeptideLoop1(loop[0], loop[1], res, resids, resnums, resns, chids, inscodes, pts,  None ,caRad, scRad, [])
        #mergeBRlists(blist, bfoll, numtrials, rlist, blist1, bfoll1, numtrials1, rlist1)


    print "BUILDERS------------------------------"
    for b in blist :
        b.describe() ; print ''
    print "---------------------------------------\n\n\n"


    print "Residues ------------------------------"
    for  k , v in res.items():
        print k , resids[k] , v
    print "---------------------------------------\n\n\n"

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
        
    knownPositions = range(len(pts))
    for b in blist :
        bop = b.getOP()
        for i in range(bop.size()) :
            if bop[i] in knownPositions :
                knownPositions.remove(bop[i])
    for ki in knownPositions : aiSC[ki] = -1 ##XXX

    radii = [0] * len(pts) ## make radii for all atoms including dummies
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else :
                try : radii[pi] = vdwr['XXX'][aname[1]]
                except KeyError : radii[pi] = vdwr['XXX']['C']
        
    ssReduction = data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    if restrGen !=None :
        xrlist,xopt= restrGen.generate(blist, aiSC,res,resids,pts)
    else :
        xrlist = []


    extraRestraints = [];    optional = []

    print "NOW PRINTING ALL RESTRAINTS:: "
    for r in rlist:
        extraRestraints.append(r)
        r.describe() ; print 
        if allopt == 1:
            optional.append(len(extraRestraints)-1)

        
    for r in xrlist:
        extraRestraints.append(r)
        r.describe(); print 
        if allopt == 1 or edopt == 1 :
            optional.append(len(extraRestraints)-1)


    for r in irlist:
        extraRestraints.append(r)
        r.describe(); print
        if allopt == 1:
            optional.append(len(extraRestraints)-1)

        

    nbuilt = 0 ; natt = 25
    from peptidebuild import ModelRenderer
    from PopulationStrategy import PopulationStrategy

    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, dummies.keys(), outpdb)
    strategy = PopulationStrategy(backtrack,popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll,None,extraRestraints , optional , modelRenderer, natt, nmodels,None)
    strategy.ranker = ranker
    strategy.execute()
    sys.exit()


def callmain() :
    from commonOptions import makeParser , parseOptions
    import optparse ; parser = optparse.OptionParser()
    import sys


    parser.add_option("--dir-xyzout", action='store', type='string', dest='dir_xyzout', help='Name of directory to create all the files during refinement. complete PATH needed')
    
    parser.add_option("--xyzin", action='store', type='string', dest='pdbfile', help='Name of PDB input structure')
    parser.add_option("--xyzout", action='store', type='string', dest='pdbout', help='Name of PDB output file',default="modelout.pdb")

    parser.add_option("--mapin", action='store', type='string', dest='mapfn', help='Input map file',default=None)
    parser.add_option("--hklin", action='store', type='string', dest='mtzfn', help='Input phased MTZ file',default=None)
    parser.add_option("--use-ca-restraints", action='store', type='int', dest='caRes', help='[0/1], 1:Yes, 0:No. Apply positional restraints on the C-alpha atom',default=1)
    parser.add_option("--use-sc-restraints", action='store', type='int', dest='scRes', help='[0/1], 1:Yes, 0:No. Apply positional restraints on the centroid of the sidechain atoms',default=1)

    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='Radius of spherical restraint on the CAlpha atom position', default=1.0)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='Radius of spherical restraint on the centroid of the sidechain atoms', default=2.0)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='Factor by which to reduce effective Van der Waals distance for sidechain atoms', default=1.0)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='Population size', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='Level of printed output [0-10], 0:Concise output log, 10: Detailed output log', default=0)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='[numsteps]X[stepsize] e.g. 4X5 will set backtrack to numsteps and stepsize to 4,5 respectively. For detailed help see Rappertk wiki/webpage', default="2X3")
    parser.add_option("--mconly", action='store', type='int', dest='mconly', help='Build mainchain only', default=None)
    parser.add_option("--sconly", action='store', type='int', dest='sconly', help='Build side chains only', default=0)
    parser.add_option("--opsax", action='store', type='int', dest='opsax', help='Reassign side chains with OPSAX [0/1] 1:Yes, 0;No', default=0)

    parser.add_option("--use-freer", action='store', type='int', dest='usefreer', help='Use FreeR set ? ', default=0)

    parser.add_option("--rotamerlib", action='store', type='string', dest='rotLib', help='Rotamer library to use when building side chains', default='PRL')        
    parser.add_option("--num-models", action='store', type='int', dest='nmodels', help='Number of models wanted ', default=1)

    ### Added params to be RAPPER-like
    parser.add_option("--a", action='store', type='float', dest='a', help='Cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='Cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='Cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='Cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='Cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='Cell spacegroup, in CCP4 notation')
    parser.add_option("--resolution", action='store', type='float', dest='resolution', help='Resolution of the data')

    parser.add_option("--FP", action='store', type='string', dest='f1label', help='Column label for FP in MTZ file', default=None)
    parser.add_option("--SIGFP", action='store', type='string', dest='sigf1label', help='Column label for sigFP in MTZ file', default=None)
    parser.add_option("--FC", action='store', type='string', dest='f2label', help='Column label for FC in MTZ file', default=None)
    parser.add_option("--PHIC", action='store', type='string', dest='phiclabel', help='Column label for PHIC in MTZ file', default=None)
    parser.add_option("--FREER", action='store', type='string', dest='freerlabel', help='Column label for FreeR in MTZ file', default=None)
    parser.add_option("--m", action='store', type='int', dest='m', help='scale FP to generate maps when MTZ file is given', default=2)
    parser.add_option("--n", action='store', type='int', dest='n', help='scale FC to generate maps when MTZ file is given', default=1)


    parser.add_option("--default-mainchain-b-factor", action='store', type='float', dest='mcBfac', help='The B-factor assigned to the newly built main chain atoms', default=20.)
    parser.add_option("--default-sidechain-b-factor", action='store', type='float', dest='scBfac', help='The B-factor assigned to the newly built side chain atoms', default=30.)
    parser.add_option("--models-get-native-bfactors", action='store', type='int', dest='nativeBfac', help='[0/1] 1:Yes 0:No. Assign B-factors of remodelled atoms to original values', default=1)

    parser.add_option("--minimum-sig", action='store', type='float', dest='minXSig', help='Minimum sigma ', default=0.25)
    parser.add_option("--detect-breaks", action='store', type='int', dest='detectBreaks', help='Detect chain breaks: [1,0] Yes:1, No:0', default=1)
    parser.add_option("--maximum-sig", action='store', type='float', dest='maxXSig', help='Maximum sigma ', default=2.0)
    parser.add_option("--cacaCutoff", action='store', type='float', dest='cacaCutoff', help='Minimum distance between adjacent Calpha atoms in order to detect a chain-break', default=5.)

    parser.add_option("--make-ed-optional", action='store', type='int', dest='edOpt', help='[0/1] 1:Yes 0:No. If 1 then positive density mainchain restraint will be made optional. If 0, then the main chain will be unconditionally forced to lie in positive density. This is primarily useful when tracing through a structure with regions in very poor (non-existent) density', default=0)
    parser.add_option("--make-all-restraints-optional", action='store', type='int', dest='allOpt', help='[0/1] 1:Yes 0:No. If 1, then all  restraints will be made optional', default=0)    


    
    parser.add_option("--ligfile", action='store', type='string', dest='ligfile', help='ligand description, see msq.ligdesc for description of MSQ in 1di9.pdb')
    parser.add_option("--around-ligand", action='store', type='float', dest='closeCutoff', help='min-dist between residue/ligand to be considered close to ligand', default=10)

    (options, args) = parser.parse_args()
    import misc;    misc.setVerbosity(options.verbose)

    if options.dir_xyzout == None :
        options.dir_xyzout = os.getcwd()+"/rappertkout"

    if not os.path.isdir(options.dir_xyzout) : os.mkdir(options.dir_xyzout)
    if (os.path.isfile(options.pdbfile)==False) :
        print "Cannot find file %s "%options.pdbfile
        import sys ; 
        sys.exit()
    shutil.copyfile(options.pdbfile, "%s/init.pdb" % (options.dir_xyzout))

    
    if (os.path.isfile(options.ligfile)==False) :
        print "Cannot find file %s "%options.ligfile
        import sys ; sys.exit()
    shutil.copyfile(options.ligfile, "%s/lig.desc" % (options.dir_xyzout))
    ligfilepath = "%s/lig.desc" % (options.dir_xyzout)

    if options.mapfn != None  :
        if (os.path.isfile(options.mapfn)==False) :
            print "Cannot find file %s "%options.mapfn
            import sys ; sys.exit()
        shutil.copyfile(options.mapfn, "%s/init.map" % (options.dir_xyzout))

    if options.mtzfn != None : 
        if (os.path.isfile(options.mtzfn)==False) :
            print "Cannot find file %s "%options.mtzfn
            print "No file in directory ", os.getcwd()
            import sys ; sys.exit()
        shutil.copyfile(options.mtzfn, "%s/init.mtz" % (options.dir_xyzout))


    os.chdir(options.dir_xyzout)
    fp = open('parameters.txt', 'w')
    fp.close()

    if (options.mconly == 1 and  options.sconly == 1) :
        print "Mainchain only and sidechain only modes cannot be used together"
        import sys ; sys.exit()

    xrayRestGen = None ; mapfilepath = None ; ranker = None
    
    if options.mtzfn != None or options.mapfn !=None : 
        if (options.a == None):
            from stump import getCRYST , getRESO
            options.a,options.b,options.c,options.alpha , options.beta , options.gamma, sg  = getCRYST(options.pdbfile)
        if (options.a == None or options.b == None or options.c == None or options.alpha== None or options.beta==None or options.gamma == None ):
            print "Please input cell paramater a, b , c , alpha, beta , gamma = ",options.a , options.b , options.c , options.alpha , options.beta  , options.gamma 

            sys.exit()
        if (options.sg == None):
            options.sg = sg
            if (options.sg == None):
                print "Please input space group " , options.sg

        from data import sgtable , long2shortHM
        if options.sg  in long2shortHM.keys():
            shortsg = long2shortHM[options.sg]
            options.sg = shortsg
        if options.sg not in sgtable.keys():
            print "Check --sg , Not recognised [%s][%d]"%( options.sg, len(options.sg))
            sys.exit()

        if (options.resolution == None):
            options.resolution = getRESO(options.pdbfile)
            if (options.resolution == None):
                print "Please input resolution " , options.resolution
                sys.exit()
            else :
                print "Resolution", options.resolution
            

        if (options.f1label == None):
            print "Please specify FP label  " , options.f1label
            sys.exit()
        if (options.sigf1label == None):
            print "Please specify SIGFP label  " , options.sigf1label
            sys.exit()
            
        xrayRestGen = None ;        ranker = None ;     cellsym = []
        from prepareChain import XrayRestraintsGenerator
        from xcheck import XrayScorer, XrayRanker
        from data import sgtable , long2shortHM
        esmin, esmax, esmean = .000, 5., .0, 
        cellsym = [ options.a, options.b, options.c, options.alpha, options.beta, options.gamma, sgtable[options.sg][0] ]

    if options.mapfn != None  :
        mapfilepath  = "%s/init.map" % options.dir_xyzout
        xrayRestGen =  XrayRestraintsGenerator(mapfilepath, "map", options.f2label, options.phiclabel, "", esmin, esmax, esmean,[], options.edOpt ) 
            
    
    elif options.mtzfn != None : 
        mtzfilepath  = "%s/init.mtz" % options.dir_xyzout
        from xray import mtzdump
        check = mtzdump(mtzfilepath,options.f1label, options.sigf1label, options.freerlabel,options.f2label, options.phiclabel)
        if check == None :
            print "One or more specified column labels not present in mtz file"
            import sys ; sys.exit()
        ccp4map = "%s%dFo-%dFC.map"%(mtzfilepath,options.m,options.n)
        mapcoeff = "%dF1-%dF2"%(options.m, options.n)


        if options.usefreer == 1 : 
            if options.freerlabel == None :
                print "If you would like to use the freer set, please specify column name for freeR"
                sys.exit()
                
        if (options.f2label == None or options.phiclabel == None ):
            print "Column labels for  FC and PHIC  not specified, will use input PDB structure to obtain FC and PHIC"
            from xray import sfall, fft
            sfall(options.pdbfile, mtzfilepath, "phased.mtz",options.resolution, options.usefreer, options.f1label,options.sigf1label,options.freerlabel)

            mtzfilepath  = "phased.mtz"
            options.f2label  = "FC"
            options.phiclabel = "PHIC"

        if options.usefreer == 1 : 
            fft(mtzfilepath, ccp4map, options.f1label,options.f2label,options.phiclabel,options.sigf1label,options.m,options.n,options.freerlabel)
        else :
            fft(mtzfilepath, ccp4map, options.f1label,options.f2label,options.phiclabel,options.sigf1label,options.m,options.n,None)
            
        xrayRestGen =  XrayRestraintsGenerator(ccp4map, "map", options.f2label, options.phiclabel, mapcoeff, esmin, esmax, esmean,[], options.edOpt ) 
        mapfilepath  = "%s/%s"%(options.dir_xyzout,ccp4map)

    else :
        print "WARNING !! No input mtz or map file , No electron density restraints will be used"

    
    if options.detectBreaks == 1:
        caca = options.cacaCutoff
    else :
        caca = None

#    if options.mtzfn != None : 
#
#        xrayRestGen = XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, options.sigXmin, options.sigXmax, options.sigXmean)
##        esmin, esmax, esmean, rcmult = .000, 5., .0, 5 

#        ranker = XrayRanker("phased.mtzFPFCPHIC_A.map", "map", "FC", "PHIC", "2F1-F2", esmin, esmax)
#        ranker.rankChildren = rcmult ; ranker.rankRecurse = 1
#        ranker.rankLeaderBuilderOnly = None ; ranker.rankGivenEnsemble = None
        
    main(options.pdbfile, ligfilepath, options.closeCutoff, options.caRad, options.scRad, options.scReduction, None,options.pdbout, options.popsize, options.backtrack, options.nmodels,xrayRestGen, mapfilepath ,ranker,options.rotLib,options.edOpt,options.allOpt)

if __name__ == "__main__" : callmain()
