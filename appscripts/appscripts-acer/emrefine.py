
import os, shutil, re,sys,random
import geometry
from xray import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS, fft, omitmap, mapman
from procrun import proc_run_exitOnError as execCmd
from xcheck import XrayScorer, XrayRanker
from builders import VecInt,VecVecFloat
from data import sgtable
from evalCAtrace import comparePhiPsiOmegaChi
from restraints import DistanceRestraint #,DihedralRestraint
##from scplacement import SCplacement
#from loopbuildFixEM import Multiloop
import prepareChain
#from prepareChain import removeSC,makeDihedralRestraints
from pdbr import protein, isAAres
import prot2res

from restraints import SphPosRestr, CentroidPosRestraint, DistanceRestraint, AngleRestraint, EDrestraint #, #DihedRestraint#,DihedralRestraintD
from prepareChainV5 import removeSC, mergeBRlists, addNdummyGly, addCdummyGly, printResAtoms, PrepareChain, prepareRNAchain, incompleteSCcorrection , incompleteMCSCcorrection  , incompleteMCcorrection , incompleteSCcorrection2



def makeHbondRestraints(ni, oi, res, resids) :
    rlist = []
    ONhbd = (1.5, 3.5) # distance
    CONhbd = (100, 200) # angle limits
    print res[ni]
    print resids[oi],res[oi]
    rlist.append( DistanceRestraint(VecInt([res[ni][' N  '],res[oi][' O  ']]), "betaHbond-NO [%s] [%s]" % (resids[ni],resids[oi]), ONhbd[0], ONhbd[1]) )
    rlist.append( AngleRestraint(VecInt([res[ni][' N  '],res[oi][' O  '],res[oi][' C  ']]), "betaHbond-NOC [%s] [%s]" % (resids[ni],resids[oi]), CONhbd[0], CONhbd[1]) )
    return rlist







def findWorseFits(XcorrZ, cutoff=0.9) :
    badkeys = []
    #vals = list(XcorrZ.values())
    #vals.sort()

    #cutoff = vals[ (len(vals)-1)/4 ]
    for k in XcorrZ.keys() :
        if XcorrZ[k] < cutoff and k[0:3] != "HOH" : badkeys.append(k)
    #print "BADKEYS", cutoff, len(badkeys), badkeys
    return badkeys

def main() :

    import optparse ; parser = optparse.OptionParser()

    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='to create all the files during refinement. it shdnt be already present.')
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='starting pdb containing a model of pdb-ligand complex')

    parser.add_option("--debug", action='store', type='int', dest='debug', help='starting pdb containing a model of pdb-ligand complex')
    parser.add_option("--hbdfile", action='store', type='string', dest='hbdfile', help='Hbond restraints to enfirce SS structure')
    parser.add_option("--buildN2C", action='store', type='int', dest='buildN2C', help='by default, build from N to C terminal. build C->N if 0.', default=1)
    parser.add_option("--map", action='store', type='string', dest='map', help='structure factors file')
    parser.add_option("--init-map", action='store', type='string', dest='inimap', help='structure factors file',default=None)
    parser.add_option("--resolution", action='store', type='float', dest='resolution', help='resolution of the data')
    parser.add_option("--test", action='store', type='float', dest='test', help='resolution of the data',default =1.)
    parser.add_option("--startcycle", action='store', type='int', dest='startcycle', help='Total number of RTK cycles',default=1)
    parser.add_option("--stopcycle", action='store', type='int', dest='stopcycle', help='Total number of RTK cycles',default=50)
    parser.add_option("--num-models-wanted", action='store', type='int', dest='nmodels', help='number of models desired, 100 attempts per model. DEFAULT 100', default=100)

    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=1)
    parser.add_option("--optional", action='store', type='int', dest='opt', help='Optional:0 , Not optional:1',default=0)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid', default=2)
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='seed for randomizing', default=None)

    parser.add_option("--start", action='store', type='int', dest='start', help='start', default=None)
    parser.add_option("--stop", action='store', type='int', dest='stop', help='stop', default=None)
    parser.add_option("--mconly", action='store', type='int', dest='mconly', help='mconly:1', default=None)
    parser.add_option("--ssfile", action='store', type='string', dest='ssfile', help='list of secondary structures, see applications/ssfile for example')
    parser.add_option("--modelAll", action='store', type='int', dest='modelAll', help='Rebuiild only loop,set tp 1',default=0)
    parser.add_option("--cutoff", action='store', type='float', dest='cutoff', help='CC cut off for flagginf mis-fit regions',default=None)
    parser.add_option("--loopOnly", action='store', type='int', dest='loopOnly', help='Rebuiild only loop,set tp 1',default=0)
    parser.add_option("--helicesOnly", action='store', type='int', dest='helicesOnly', help='Rebuiild only helices',default=0)
    parser.add_option("--pre-process", action='store', type='int', dest='preProcess', help='Preprocessing bad regions to overcome possible local minima True:1 , False:0',default=0)
    parser.add_option("--restoreSS", action='store', type='int', dest='restoreSS', help='Restore SS ? ',default=0)
    parser.add_option("--use-midloop", action='store', type='int', dest='midloop', help='Use bridge building, True:1 , False:0',default=1)


    parser.add_option("--rank", action='store', type='int', dest='ranker', help='Rank, True:1 , False:0',default=1)

### Additional homo infor

    parser.add_option("--pir", action='store', type='string', dest='pirfile', help='pir file')
    parser.add_option("--ensemble", action='store', type='string', dest='ensemble', help='make Ensemble',default=1)
    parser.add_option("--clr", action='store', type='string', dest='clrfile', help='clr file')
    parser.add_option("--use-choral", action='store', type='int', dest='choral', help='use choral',default=0)
    parser.add_option("--pdbext", action='store', type='str', dest='pdbext', help='extension of template files [atm/ brk / pdb ] ?? ',default="pdb")
    parser.add_option("--pdb-path", action='store', type='str', dest='pdbpath', help='directory containning the pdb files ',default='/home/anjum/')
    parser.add_option("--band", action='store', type='int', dest='band', help='band ss elements ? ? ', default=None)
    parser.add_option("--short_sidechain_restraint_threshold", action='store', type='float', dest='ssc_threshold', help='Restraint radii on short side chain atoms', default=1.5)
    parser.add_option("--mainchain_restraint_threshold", action='store', type='float', dest='mc_threshold', help='Restraint radii on MC atoms', default=1.5)

    parser.add_option("--params", action='store', type='string', dest='paramfile', help='parameter file')
    parser.add_option("--modelfrag", action='store', type='string', dest='modelfrag', help='modelfrag?')

    (options, args) = parser.parse_args()
    import misc
    misc.setVerbosity(options.verbose)
    randomize(options.randomize)
    if not os.path.isdir(options.scratchdir) : os.mkdir(options.scratchdir)
    if options.ssfile !=None:
        shutil.copyfile(options.ssfile,"%s/ssfile" % (options.scratchdir))
    if options.hbdfile !=None:
        shutil.copyfile(options.hbdfile,"%s/%s" % (options.scratchdir,options.hbdfile))



    shutil.copyfile(options.pdbfile, "%s/model%d.pdb" % (options.scratchdir,options.startcycle-1))

    ssList = []  ; resList  = [] 
    if options.clrfile !=None or options.pirfile !=None :
        from homologyEM import main as homologymain
        print options.clrfile
        res, resids, resnums, resns, chids, inscodes, pts , resList ,start ,stop, params , msa  = homologymain(options.scratchdir, options.pirfile , options.clrfile, options.paramfile, None, None , options.choral, None, None , options.verbose , options.mconly , options.pdbext , options.pdbpath , None , None , None , None , None , options.mc_threshold,options.ssc_threshold,None, 1)

    else:
        os.chdir(options.scratchdir)
    ### Values of params and i/o file names 
    guidedsampling = None
    multiPrepC = prepareChain.PrepareChain("PRL")
    popsize = options.popsize
    #totcycles = options.totcycles
    if options.cutoff == None:
        if options.resolution == 10.0:
            options.cutoff = 0.98
        if options.resolution == 6.0:
            options.cutoff = 0.97
            

    from restraints import EDrestraint ; EDrestraint.setPenalty(5.) ; EDrestraint.setScatRad(1.) ;
    esmin, esmax, esmean, rcmult, xscoreCutoff = .000, 5., .0, 1, options.cutoff
    
    ranker = XrayRanker(options.map, "map", "FC", "PHIC", "F1", esmin, esmax)
    xrayRestGen = []
    #xrayRestGen.append( prepareChain.XrayRestraintsGenerator(options.map, "map", "FC", "PHIC", "F1", esmin, esmax, esmean, ["SCLbuilder","ChiBuilder","CBbuilder"], ) )
    xrayRestGen.append( prepareChain.XrayRestraintsGenerator(options.map, "map", "FC", "PHIC", "F1", esmin, esmax, esmean, [], ) )
    
    ## make working directory and copy files over

        
    

    ### Read in SS file and make helical restraints
    helices , sheets = [], []
    helixind = []
    if options.ssfile !=None:

        modelIn = "model%d.pdb" % (options.startcycle-1) 
        prot = protein(modelIn, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        rev_resnum  = {}
        for k,v in resnums.items() : rev_resnum[int(v)] = int(k)
        from betaBuilder import parseSSfileAK
        helices , sheets = [], [] 
        helices , sheets = parseSSfileAK(options.ssfile,resids,resnums)

        temp = [] ;     helixind = []
        mcMissInds = []; scMissInds = [] ; caMissInds = [];  badids = []
        import checkProtChainsV4
        mcmiss, scmiss, chainBreaks1 = checkProtChainsV4.check(res, resids, resnums, resns, chids, inscodes, pts ,  4.0 ,[],0)
        for helix in helices:
            hstartkey = int(helix[1][0]-2)  
            hstopkey = int(helix[1][1]+2)  
            for hs in range(hstartkey,hstopkey+1):
                helixind.append(hs)
        if options.mconly != 1:
            for sc in scmiss :
                if sc in helixind: 
                    scMissInds.append(sc)
                
        for mc in mcmiss :
            if mc in helixind : 
                if ' CA ' not in res[mc].keys() :
                    caMissInds.append(mc)

                mcMissInds.append(mc)
                if options.helicesOnly == 1 :
                    misspts = incompleteMCcorrection(res, resns, pts , mc,resids)

        from peptidebuild import ModelRenderer
        modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [], "model0.pdb")
        modelRenderer.render(pts)

        for k,v in res.items():
            print resids[k], res[k]
#        import sys ; sys.exit()

        temp = [] ;     helixind = []
        for helix in helices:
            hstartkey = int(helix[1][0])
            hstopkey = int(helix[1][1])
            temp.append(makeHelix(hstartkey,hstopkey,res,resids))
            for hs in range(hstartkey,hstopkey+1):
                helixind.append(hs)
        for k in temp:
            for j in k:
                ssList.append(j)


                            


    ### Read in SS file and hnd file and make sheet restraints
        if options.hbdfile != None:
            NOpartners = getNOpartners(options.hbdfile)
            for strands, dirs, corrInds in sheets : 
                starts, ends = [], []
                for beg,end in strands :
                    starts.append(beg) ;
                    ends.append(end)
                            
                #for si in range(len(starts)) :
                    #temp = makeDihedralRestraints(starts[si],ends[si],res,resids,"beta")
                    #for restr in temp:
                    #    ssList.append(restr)
                
                print "NO",NOpartners
                print "res" , rev_resnum.keys()
                for k,v in NOpartners.items():
                    if  k in rev_resnum.keys() : #and v in rev_resnum.keys():
                        npartner = rev_resnum[k]
                        print "N" , npartner
                        for vv in v:
                            if vv in rev_resnum.keys():
                                opartner = rev_resnum[vv]
                                for si in range(len(starts)) :
                                    if npartner in range(starts[si],ends[si]+1):
                                        for pi in range(len(starts)) :
                                            if opartner in range(starts[pi],ends[pi]+1) and si != pi:
                                                hh = makeHbondRestraints(npartner, opartner, res, resids) 
                                                for h in hh :
                                                    ssList.append(h)
                                                    print "APP" , npartner , opartner



    print "options.start",options.startcycle
    for cycle in range(options.startcycle-1 , options.stopcycle):
        print "cyc",cycle
        rtkmodel = "model%d.pdb" % (cycle+1) 
        modelIn = "model%d.pdb" % (cycle) 
        FCEMmap =  modelIn +".FC.map"
        if options.helicesOnly != 1 :
            continue
        elif not  os.path.isfile(modelIn) :
            print "No file"
            print "No file", modelIn,cycle
            print os.getcwd()
            sys.exit()
        else :
            FCEMmap =  modelIn +".FC.map"
            FCEMmap = makeFCEMmap(modelIn,options.resolution,options.map)
            
    
            prot = protein(modelIn, read_hydrogens=0, read_waters=0, read_hets=0)
            res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
            if options.mconly : res, pts = removeSC(res, pts, res.keys())
            rev_resnum  = {}
            for k,v in resnums.items() : rev_resnum[int(v)] = int(k)
            
            badresids = None

            if options.helicesOnly == 1:
                loops , looptypes, midloops, badloops = [],[],[],[]
                for hh in  helices:
                    hstart = hh[1][0]
                    hstop = hh[1][1]
                    loops.append([hstart,hstop])
                    if options.modelfrag == "1" :
                        looptypes.append("NtoC")
                    else :
                        looptypes.append("loop")

                if options.midloop == 1 :
                    midloops = []
                    for ll in range(0,len(loops)):
                        midloop = int(loops[ll][0]) + ( (int(loops[ll][1])-int(loops[ll][0]) ) /2)
                        midloops.append(midloop)
                        
                else :
                    midloops = []
                    for ll in range(0,len(loops)):
                        midloops.append(None)



                    
                badloops.append(loops)
                badloops.append(looptypes)
                badloops.append(midloops)
                XUcaRad = [] ; XUscRad = []
                for l in looptypes:
                    XUcaRad.append(options.caRad)
                    XUscRad.append(options.scRad)
    
    
                ml = Multiloop(modelIn, badloops , options.mconly, options.caRad ,options.scRad, options.scReduction, guidedsampling, popsize, options.backtrack,options.nmodels,rtkmodel, xrayRestGen, multiPrepC,options.buildN2C,XUcaRad , XUscRad)
                
                ml.ranker = XrayRanker(options.map, "map", "FC", "PHIC", "F1", esmin, esmax)
                ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = None

                if options.ranker == 1 :
                    ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble =  1
                
                ml.cellsym = None
                restraintList = ssList+resList
                nb,notbuilt = ml.run('CA',restraintList,0)
    if options.helicesOnly  == 1:
        sys.exit()
    for cycle in range(options.startcycle -1 , options.stopcycle):
        if cycle == 0 and options.test == 2. :
            oldcarad = options.caRad
            oldscrad = options.scRad
            options.caRad = 10000.
            options.scRad = 10000.
        elif cycle > 0 and options.test == 2. :
            options.caRad = oldcarad
            options.scRad = oldscrad
        rtkmodel = "model%d.pdb" % (cycle+1) 
        modelIn = "model%d.pdb" % (cycle) 
        FCEMmap =  modelIn +".FC.map"
        if not  os.path.isfile(modelIn) :
            print "No file", modelIn
            print os.getcwd()
            sys.exit()
        else :
            FCEMmap =  modelIn +".FC.map"
            FCEMmap = makeFCEMmap(modelIn,options.resolution,options.map)
            
    
            prot = protein(modelIn, read_hydrogens=0, read_waters=0, read_hets=0)
            res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
            if options.mconly : res, pts = removeSC(res, pts, res.keys())
            rev_resnum  = {}
            for k,v in resnums.items() : rev_resnum[int(v)] = int(k)
            
            badresids = None
            
            if options.helicesOnly == 1:
                loops , looptypes, midloops, badloops = [],[],[],[]
                for hh in  helices:
                    hstart = hh[1][0]
                    hstop = hh[1][1]
                    loops.append([hstart,hstop])
                    looptypes.append("loop")
                    midloops.append(None)
                    print "HELIX ",hstart , hstop
                    
                badloops.append(loops)
                badloops.append(looptypes)
                badloops.append(midloops)
                XUcaRad = [] ; XUscRad = []
                for l in looptypes:
                    XUcaRad.append(options.caRad)
                    XUscRad.append(options.scRad)
    
    
                ml = Multiloop(modelIn, badloops , options.mconly, options.caRad ,options.scRad, options.scReduction, guidedsampling, popsize, options.backtrack,options.nmodels,rtkmodel, xrayRestGen, multiPrepC,options.buildN2C,XUcaRad , XUscRad)
                
                ml.ranker = XrayRanker(options.map, "map", "FC", "PHIC", "F1", esmin, esmax)
                ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1

                if options.ranker == 1 :
                    ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble =  1
                
                ml.cellsym = None
                restraintList = ssList+resList
                nb,notbuilt = ml.run('CA',restraintList,0)                
            elif options.modelAll == 1:
                badresids = []
                for residue in resids.values():
                    badresids.append(residue)
                    
    
            elif options.start != None and options.stop!=None :
                loop = []
                for k,v in resids.items():
                    if int(resnums[k]) > options.start and int(resnums[k]) < options.stop:
                        loop.append(v)
                badresids = loop
    
#            elif options.loopOnly == 0  and  (cycle <= 19 or cycle >= 31):
            elif options.loopOnly ==  0 :
                xscorer = XrayScorer(None, xscoreCutoff)
                badmc = xscorer.score(modelIn, options.map, FCEMmap, "Anything", "Anything", "Anything", "mc")
                badpept = xscorer.score(modelIn, options.map, FCEMmap, "Anything", "Anything", "Anything", "pept")
                badmcsc = xscorer.score(modelIn, options.map,FCEMmap, "Anything", "Anything", "Anything", "mcsc")
                badscs = xscorer.score(modelIn, options.map, FCEMmap, "Anything", "Anything", "Anything", "sc")
        
                #badresids = list ( set(badmc+badpept+badscs) )
                #badresids = list ( set(badmc) )
                #badresids = list ( set(badmcsc) )
                badresids = list ( set(badmcsc+badmc+badpept) )
                if badresids == []:
                    print "Nothing flagged to be fixxed ..exiting"
                #sys.exit(1)
                    
            elif options.loopOnly == 1 or  (cycle > 19 and cycle < 31): ## ie. 20 - 30 
                ## only residues in loop regions are rebuilt
                badresids = []
                starts, ends = [], []
                for strands, dirs, corrInds in sheets : 
                    for beg,end in strands :
                        starts.append(beg) ;
                        ends.append(end)
                for residue in resids.values():
                    badresids.append(residue)
        
                for hh in  helices:
                    hstart = hh[1][0]
                    hstop = hh[1][1]
                    for hresidue in range(hstart,hstop+1):
                        badresids.remove(resids[hresidue])
                        print "Removing helical", resids[hresidue]
                for si in range(len(starts)):
                    for sresidue in range(starts[si],ends[si]):
                        badresids.remove(resids[sresidue])
                        print "Removing strand", resids[sresidue]
                        
            else :
                badresids = None

            print "------------------------------------------------------------------------------"

            
            print "residue to be rebuit",badresids
            
        
        
            ## check band end points  are not in the middle of a helix
            from loopbuild import locateRegionsRandomize
            loops, looptypes = locateRegionsRandomize(resids, chids, badresids) # locate loops and order them randomly
    
        
        
            if options.preProcess  == 1:
                loops.sort()
                for ll in range(0,len(loops)-1):
                    end1 = int(loops[ll][1])
                    start2 = int(loops[ll+1][0])
                    if int(start2) - int(end1) < 10:
                        for xx in range((end1)+1,start2):
                            badresids.append(resids[xx])
                loops, looptypes = locateRegionsRandomize(resids, chids, badresids) # locate loops and order them randomly
            if options.midloop == 1 :
                midloops = []
                for ll in range(0,len(loops)):
                    midloop = int(loops[ll][0]) + ( (int(loops[ll][1])-int(loops[ll][0]) ) /2)
                    if midloop in helixind : 
                        midloops.append(None)
                        #for fib in range((int(loops[ll][1])-int(loops[ll][1]))/2):
                        #    if (midloop+fib) not in helixind:
                        #        midloop =  midloop+fib
                        #        break
                        #    if (midloop-fib) not in helixind:
                        #        midloop =  midloop-fib
                        #        break
                    else:
                        midloops.append(midloop)
                        
            else :
                midloops = []
                for ll in range(0,len(loops)):



                    midloops.append(None)
        
        
            ### score each band to choose Ca restraint radii for that band
            XUcaRad = [] ; XUscRad = []
            badloops = []
            loopind  = []

            if cycle == 0 :
                origCarad = options.caRad
                origScrad = options.scRad
            elif cycle > 15 and cycle < 30:
                options.caRad =  10 *  origCarad
                options.scRad =  10 *  origScrad 
            else :
                options.caRad =   origCarad
                options.scRad =  origScrad 
                
            for ll in loops:
                loopall = []
                crds = []
                alist = []
                lstart= ll[0]
                lstop= ll[1]
                loopind.append(lstart)
                loopind.append(lstop)
                count = 0
                for lx in range (lstart,lstop+1): 
                    loopall.append(lx)

                XUcaRad.append(options.caRad)
                XUscRad.append(options.scRad)



            badloops.append(loops)
            badloops.append(looptypes)
            badloops.append(midloops)

            print "badloops",badloops[0]
            print "badloops",badloops[1]
            print "badloops",badloops[2]
            
            ml = Multiloop(modelIn, badloops, options.mconly, options.caRad ,options.scRad, options.scReduction, guidedsampling, popsize, options.backtrack, 1,rtkmodel, xrayRestGen, multiPrepC,options.buildN2C,XUcaRad,XUscRad,options.modelfrag)
            
            ml.ranker = XrayRanker(options.map, "map", "FC", "PHIC", "F1", esmin, esmax)
            ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = None
            
            if options.ranker == 1 and options.test != 2.0 :
                ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble =  1

            elif options.ranker == 1 and options.test == 2.0 and cycle > 0:
                ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble =  1

            ml.cellsym = None

            #ml.ranker = None
            
            if cycle >= 14 and cycle <= 19 : 
                restraintList = ssList

            if cycle > 19 and cycle <= 30 : 
                restraintList = []
                
            else : ## for cycle > 30 , loose cm will come into play
                restraintList = ssList + resList
                



            #if cycle >= 15 and cycle <=25: 
            #    ml.ranker.rankGivenEnsemble =  None
            #else :
            #    ml.ranker.rankGivenEnsemble =  1
            if options.debug == 1:
                erlist = []
                for k, v in res.items():
                    resmcpts = []  ;                    resscpts = []
                    for a , b in v.items():
                        if a in [' CA ', ' N  ',' C  ',' O  '] : 
                            resmcpts.append(b)
                            descr1 = "ED on %s ca"%resids[k] 
                        
                        else:
                            resscpts.append(b)
                            descr2 = "ED on %s sc"%resids[k]
                    if len(resmcpts) !=0 :
                        erlist.append( EDrestraint.makeEDrestraintFromMap(VecInt(resmcpts), descr1, options.map, esmin, esmax, esmean) )
                    if len(resscpts) !=0 :                    
                        erlist.append( EDrestraint.makeEDrestraintFromMap(VecInt(resscpts), descr2, options.map, esmin, esmax, esmean) )
                        
                for er in erlist :
                    print er.describe();
                    print er.satisfied(VecVecFloat(pts))
                    print "========================="
                for sr in resList :
                    print sr.describe();sr.satisfied(VecVecFloat(pts))
                    
                for sr in ssList :
                    print sr.describe();sr.satisfied(VecVecFloat(pts))

                sys.exit()

            nb,notbuilt = ml.run('CA',restraintList,0)




            
    ### preFinal round to rebuild helices that have become loopxs

                
    loops , looptypes, midloops, badloops = [],[],[],[]
    for hh in  helices:
        hstart = hh[1][0]
        hstop = hh[1][1]
        loops.append([hstart,hstop])
        looptypes.append("loop")
        midloops.append(None)
                
    badloops.append(loops)
    badloops.append(looptypes)
    badloops.append(midloops)
    XUcaRad = [];     XUscRad = []; 
    for l in looptypes:
        XUcaRad.append(5.0)
        XUscRad.append(5.0)

    if os.path.isfile("pre.ensemble.pdb") or options.restoreSS == 0:
        if options.restoreSS == 0:
            shutil.copyfile(rtkmodel,"pre.ensemble.pdb")
        print "Done..passing"
        
    else :
        modelIn = rtkmodel
        rtkmodel = "pre.ensemble.pdb"

        ml = Multiloop(modelIn, badloops, options.mconly, options.caRad ,options.scRad, options.scReduction, guidedsampling, popsize, options.backtrack,1 ,rtkmodel, xrayRestGen, multiPrepC,options.buildN2C,XUcaRad,XUscRad)
        if options.loopOnly != 1 :
            ml.ranker = XrayRanker(options.map, "map", "FC", "PHIC", "F1", esmin, esmax)
            ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
        if options.ranker == 1 :
            ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble =  1

        ml.cellsym = None
        nb,notbuilt = ml.run('CA',ssList,1)
    
    
        #tryRebuilding = 0
        #rebuildAtt = 2
        
    if options.ensemble == 0 :
        sys.exit()


    rtkmodel = "ensemble.pdb"  
    modelIn = "pre.ensemble.pdb"

    XUcaRad = [1.5] ;     XUscRad = [3.]
    
    xrayRestGen = []
    xrayRestGen.append( prepareChain.XrayRestraintsGenerator(options.map, "map", "FC", "PHIC", "F1", esmin, esmax, esmean, ["SCLbuilder","ChiBuilder","CBbuilder"],0) )

    ml = Multiloop(modelIn,None, options.mconly, options.caRad ,options.scRad, options.scReduction, guidedsampling, popsize, options.backtrack, options.nmodels*10,rtkmodel, xrayRestGen, multiPrepC,options.buildN2C,XUcaRad,XUscRad)

    ml.ranker = XrayRanker(options.map, "map", "FC", "PHIC", "F1", esmin, esmax)
    ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
    if options.ranker == 1 :
        ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble =  0
        
    ml.cellsym = None
    nb,notbuilt = ml.run('CA',ssList,0)
    
        




def restoreHelices(helices,rtkmodel):
    ## this routine compares supplied SS with SS retained in rtk models,
    ## those SS that are not retaiend are returned fro rebuilding
    from procrun import proc_run_exitOnError as execCmd
    cmd1 = "/home/anjum/downloads/promotif/promotif.scr %s"%(rtkmodel)
    cmd0 = "/home/anjum/downloads/promotif/promotif.scr initial.pdb"
    execCmd(cmd1)
    execCmd(cmd0)
    import re

    p = re.compile('model\d+')
    m = p.findall(rtkmodel)
    hlxfile1 = m[0] + ".hlx"
    hlxfile0 = "initial" + ".hlx"
    from stumpSS import getHees2
    newhelices = getHees2(hlxfile1,rtkmodel)
    oldhelices = getHees2(hlxfile0,rtkmodel)

    
    




def calcCC(crdlist,pts,mapfile,esmin,esmax):
    crdlist1 = VecVecFloat(crdlist)
    pis = VecInt( range(crdlist1.size()) )
    from restraints import EDrestraint
    er  = EDrestraint.makeEDrestraintFromMap(pis, "Xranker", mapfile, esmin, esmax, esmax)
    xcorr = EDrestraint.findBadfit1(mtzfn, f1label, pts, ptinds)
    print "xcor",xcorr,resids[ri]
    return xcorr



def makeFCEMmap(pdbfile,resolution,imap,scaling=1.00):
    from procrun import proc_run_exitOnError as execCmd
    
    
    import re
    p = re.compile('Original')
    mrcfile = pdbfile+".FC.mrc"
    cmrcfile = pdbfile+"origin.FC.mrc"
    mapfile = pdbfile+".FC.map"
    
    cmd = "proc3d %s > scratch"%imap
    execCmd(cmd)
    fbox,box = None, None
    
    for l in open("scratch", 'r').readlines() :
        m = p.match(l)
        if m!=None:
            print l
            p2 = re.compile('\d+')
            box =  p2.findall(l)[0]
            fbox = float(box)


######
    p = re.compile("\s+Sampling ")
    

    sampling = None
    for l in open("scratch", 'r').readlines() :
        m = p.match(l)
    if m!=None:
        p2 = re.compile('\d+.\d+')
        sline =  p2.findall(l)[0]
        sampling  = float(sline)
        print sampling
    else :
        print l
        
    assert(fbox!=None)
    cmd= "pdb2mrc %s %s apix=%f res=%f box=%f center"%(pdbfile,mrcfile,sampling,resolution,fbox)
    execCmd(cmd)
    cmd= "proc3d %s %s origin=0,0,0"%(mrcfile,cmrcfile)
    execCmd(cmd)
    os.rename(cmrcfile, mapfile)
    os.remove(mrcfile)

    return mapfile







def makeHelix(start,end,res,resids):
    rlist = []
    if start < end :
        for ind in range(start,end+1) :
            if ind <= end-3 :
                rlist += makeHbondRestraints(ind+3, ind, res, resids)
            if ind <= end-2 :

                rlist.append( DistanceRestraint(VecInt([res[ind][' CA '],res[ind+2][' CA ']]), "helix-CA-i-i+2 dist-restraint", 5.0, 6.0) ) # 5.1-5.8

            desc = "Dihedral Phi restraint on [%s]" %  resids[ind]

            #C(i-1) -N(i) - CA(i) - C(i)

            ai1 = res[ind-1][" C  "]
            ai2 = res[ind][' N  ']
            ai3 = res[ind][' CA ']
            ai4 = res[ind][' C  ']

            alist = []
            alist.append(ai1)
            alist.append(ai2)
            alist.append(ai3)
            alist.append(ai4)

            #rlist.append(DihedralRestraint(VecInt(alist),desc,-110.0,-40.0))
            

            desc = "Dihedral Psi restraint on [%s]" %  resids[ind]
            #N(i) - CA(i) - C(i)- N(i+1)
            ai1 = res[ind][' N  ']
            ai2 = res[ind][' CA ']
            ai3 = res[ind][' C  ']
            ai4 = res[ind+1][' N  ']

            alist = []
            alist.append(ai1)
            alist.append(ai2)
            alist.append(ai3)
            alist.append(ai4)
            
            
            #rlist.append(DihedralRestraint(VecInt(alist),desc,-75.0,0.0))

            

    return rlist


def randomize(doRandomize) :
    import misc, random
    if doRandomize :
        random.seed(doRandomize)
        misc.RanGen.instance().seedme( doRandomize )
    else : 
        random.seed(1973) # for Paul!
        misc.RanGen.instance().seedme(1942) # for Tom!



def getNOpartners(hbdfile):
    ip = open(hbdfile, 'r')
    fileList = ip.readlines()
    corr = {}
    for k in fileList :
        p1 = int(k[5:10])
        p2 = int(k[23:27])
        
        if p1 in corr.keys():
            val = corr[p1]
            val.append(p2)
            corr[p1] = val
        else : 
            corr[p1] = [p2]

            
        
    return corr
        




def fixGMXop(pdbfile) :
    from pdbr import isPdbAAline, isPdbAtomLine, line2atomname, line2resn, changeResn, changeAtomname
    lines = []
    for l in open(pdbfile, 'r').readlines() :
        if not isPdbAtomLine(l) : lines.append(l)
        elif line2resn(l) == "CU " : lines.append( changeResn(l, " CU") )
        elif line2resn(l) == "MSE" and line2atomname(l) == " SE " : lines.append( changeAtomname(l, "SE  ") )
        elif line2resn(l) == "ILE" and line2atomname(l) == " CD " : lines.append(re.sub(" CD ", " CD1", l))
        elif line2atomname(l) == " OT1" : lines.append(re.sub(" OT1", " O  ", l))
        elif line2atomname(l) == " OXT" : continue
        else : lines.append(l)
    op = open(pdbfile, 'w')
    for l in lines : op.write(l)
    op.close()


if __name__ == "__main__" :
    main()
    import sys ; sys.exit(0)
    #from scplacement import SCplacement
    import prepareChain
    scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 1, 1
    badresids = ["VAL   85 ", "ASP   86 ", "TYR   68 ", "TYR   90 ",],
    #SCplacement("premodel2.pdb", 0.5, "mmm.pdb", "dotfile", useDEE, "phased1.mtz2fofc.map", "FP", "FC", "PHIC", "2F1-F2", 0, 5, None, useGivenRot,
    #    badresids, scPrepC).run()
    import sys ; sys.exit(0)
    replaceWaters("model1.pdb", "rtk0.map")
