import math, os, sys, string , shutil , sys
from xray import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS, fft, omitmap, mapman
from builders import VecInt, VecFloat, VecVecFloat
from geometry import MapIntIntFloat, CAtraceGH, calcDist
from pdbr import protein, line2resid, line2resn, line2crd, line2atomname, line2resnum, line2resic, line2chid, makeResid
import prot2res
import data
from restraints import EDrestraint
from prepareChain import  incompleteSCcorrection , incompleteMCSCcorrection ,  incompleteMCcorrection
import misc
from prefRapperV6 import randomize
from data import sgtable , long2shortHM
import prepareChainV5
from protmodel import removeMODEL, adjustBfac2 , copyNonprotein
from copyheader import copyheader
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
            if ' CA ' not in res[ind2].keys() or ' CA ' not in res[ind1].keys() :
                continue
            if calcDist( VecFloat(pts[res[ind1][' CA ']]), VecFloat(pts[res[ind2][' CA ']]) ) > (cutoff + 7.) : continue
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


def main(pdbfile,pdbout,dir_xyzout,mapfn,mtzfn,caRes,scRes,caRad,scRad,scReduction,popsize,verbose,backtrack,rotlib,nmodels,mconly,sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freerlabel,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt,mapformat , chid12, closeCut) : 


    from checkfornone import fixdirnames
    misc.setVerbosity(verbose) ; randomize(1)
    if dir_xyzout == None : print "Directory for output files --dir_xyzout needs to be set. ",dir_xyzout; import sys ; sys.exit()
    if  os.path.isfile(dir_xyzout) : print dir_xyzout,"Is a file, Please rename output directory" ;  import sys ; sys.exit()

    dir_xyzout = fixdirnames(dir_xyzout)
    head,tail  = os.path.split(dir_xyzout)

    if head and os.path.isdir(head):
        print "Checking for directory ............",head
    else:
        print "%s does not exists"%head
        import sys ; sys.exit()

    if not os.path.isdir(dir_xyzout) :
        os.mkdir(dir_xyzout)

    os.chdir(dir_xyzout)


    if pdbout == None : print "Output filename needs to be set."; import sys ; sys.exit()    
    if (os.path.isdir(pdbout) == True )  : print "Output Coordinate file %s is a directory"%pdbout; import sys ; sys.exit()
    pdbout = fixdirnames(pdbout)
    head,tail  = os.path.split(pdbout)
    if head != "" : pdbout = tail

    if (os.path.isfile("%s/%s" % (dir_xyzout,"modelinit.pdb"))==True) :
        os.remove("%s/%s" % (dir_xyzout,"modelinit.pdb"))
    if  pdbfile == None : print "Please specify input coordinates file %s " % pdbfile ; import sys ;        sys.exit()
    if (os.path.isdir(pdbfile) == True )  : print "Input Coordinate file %s is a directory"%pdbfile; import sys ; sys.exit()
    pdbfile = fixdirnames(pdbfile)
    if (os.path.isfile(pdbfile)==False)  :        print "Cannot find file %s " %pdbfile ; import sys ;        sys.exit()



    shutil.copyfile(pdbfile, "%s/%s" % (dir_xyzout,"modelinit.pdb"))
    modelIn =  "modelinit.pdb" ;    rtkmodel = pdbout 


    #### Initializations  ####
    opsaxoff = None ; guidedsampling = None;    scvdwr = scReduction ;   
    xrayRestGen = [] ;    badresids = [] ;         nullres = [] ; missingpts = {}
    misc.setVerbosity(verbose) 
    randomize(1)
    str2bool = {"False": 0 , "True":1}
        

    ###########     CHECKS  FOR INTEGRITY OF INPUT VALUES                    #############################

    mconly =  str2bool[mconly]
    sconly =  str2bool[sconly]
    opsax =  str2bool[opsax]
    usefreer =  str2bool[usefreer]
    nativeBfac =  str2bool[nativeBfac]

    caRes =  str2bool[caRes]
    scRes =  str2bool[scRes]
    edOpt =  str2bool[edOpt]
    allOpt =  str2bool[allOpt]
    

    if mconly == 0 :        mconly = None
    if allOpt == 1 :        edOpt = 1
    if caRes == 0 :         scRes = 0
    

    ########### CHECK FOR AVAIALABLE CELL PARAMETERS, SPACE GROUP and RESOLUTION

    if (mapfn !=  None or mtzfn !=None ):
        from stump import getCRYST , getRESO
        useomitmap,usecnsmap ,useccp4map  = None,None, None
        if (a == None or b == None or c == None or alpha == None or beta == None or gamma == None) :
            print "Getting cell paramters from coordinate file....."
            a,b,c,alpha , beta , gamma,d1  = getCRYST(pdbfile)
            if (a == None or b == None or c == None or alpha == None or beta == None or gamma == None ):
                print "CRYST card cannot be read from coordinate file. Please input cell paramaters --a, --b , --c , --alpha, --beta , --gamma = ",a , b , c , alpha , beta  , gamma 
                import sys ; sys.exit()
            print "Found a b c alpha beta gamma  ", a , b , c , alpha , beta  , gamma 
        acell = a ; bcell = b ; ccell = c
        if sg == None : 
            print "Getting space group from coordinate file....."
            d1,d2,d3,d4 , d5 , d6, sg  = getCRYST(pdbfile)
            if sg == None : 
                print "Please input space group using --sg" , sg ; import sys ; sys.exit()
        ss = ""
        for sg1 in sg:
            if sg1 in ["\n","\t","\s"]:
                continue
            else :
                ss = ss+sg1
        sg = ss
        if sg  in long2shortHM.keys(): shortsg = long2shortHM[sg];  sg = shortsg
        if sg not in sgtable.keys(): print "Check --sg , Not recognised [%s]"%sg ;            import sys ; sys.exit()
        print "Setting Space Group to",sg
        
        if resolution == None : 
            print "Getting resolution limit from coordinate file........"
            resolution = getRESO(pdbfile)
            if (resolution == None):
                print "Please input resolution using --resolution value, currently" , resolution ; import sys ; sys.exit()
            print "Resolution = [ " , resolution, " ] "



    ############### CHECK RESTRAINT INFORMATION  ######################


    if caRes == 0 : print "\n********** WARNING!! In loop building mode, No C-alpha positional restraints will be used";caRad = None; 
    if scRes == 0 : print "********** WARNING!! No sidechain centroid positional restraints will be used"; scRad = None


    from protmodel import getChains
    from xray import mapdump
    
    if  mapfn != None and mtzfn != None :
        print "********** WARNING!! both mtz and map file given, please only give either MTZ file or MAP file  *************"
        import sys ; sys.exit()

    if  mapfn != None:
        if (os.path.isfile("%s/%s" % (dir_xyzout,"init.map"))==True) :
            os.remove("%s/%s" % (dir_xyzout,"init.map"))        

        esmin, esmax, esmean, rcmult  = minXSig, maxXSig, .0, 5 
        mapfilepath = mapfn
        if (os.path.isfile(mapfilepath)==False) :
            print "Cannot find file %s " %mapfn ; import sys ;            sys.exit()
        if  "map" not in mapfn : 
            print "Cannot understand i/p map format,the file should be  (*.map) format " ; import sys ; sys.exit()
        if  mapdump(mapfn)  == None :
            print "Check format of map file, needs to be CCP4 format";
            import sys ;              sys.exit()            

        shutil.copyfile(mapfilepath, "%s/init.map" % (dir_xyzout))
        ccp4map = mapfn

    elif mtzfn != None :
        if (os.path.isfile("%s/%s" % (dir_xyzout,"init.mtz"))==True) :
            os.remove("%s/%s" % (dir_xyzout,"init.mtz"))            
        mtzfilepath= mtzfn
        esmin, esmax, esmean, rcmult = minXSig, maxXSig, .0, 5 
        
        if (os.path.isfile(mtzfilepath)==False) :
            print "Cannot find file %s " %mtzfn ; import sys ;sys.exit()
        if  "mtz" not in mtzfn : 
            print "Cannot understand i/p HKLIN format,the file should be in mtz (*.mtz) format " ; import sys ; sys.exit()


        shutil.copyfile(mtzfilepath, "%s/init.mtz" % (dir_xyzout))
        mtzfilepath = "%s/init.mtz" % (dir_xyzout)


        ################   COPY MtZ COLOUMN LABELS  #########################################
        from xray import mtzdump

        if (f1label == None):
            print "Please specify FP label  " , f1label ; import sys ; sys.exit()
        if (sigf1label == None):
            print "Please specify SIGFP label  " , sigf1label ; import sys ; sys.exit()
            
        if usefreer == 1 : 
            if freerlabel == None :
                print "If you would like to use the freeR set, please specify column name for freeR"
                import sys ; sys.exit()
            
        
        if  mtzdump(mtzfn,f1label,sigf1label,freerlabel,f2label,phiclabel) == None :
            print "Check structure factor file, format needs to be MTZ"
            import sys ;              sys.exit()            
            

        if (f2label == None or phiclabel == None ):
            print "Column labels for  FC and PHI  not specified, will use input coordinate structure to obtain FC and PHI"

            
            if sfall(pdbfile, mtzfn, "phased.mtz",resolution,usefreer,f1label,sigf1label,freerlabel) == None :
                
                print "Structure factor file cannot be phased , please enter column labels for FC and PHI using keywords --FC and --PHI"
                import sys ; sys.exit()
            
            mtzfn  = "phased.mtz"
            f2label  = "FC" ; phiclabel = "PHIC"



        ### MAKE UNWEIGHTED DIFF MAP : todo sigmaA weighted maps #########################

        from xray import fftnew

        if mapformat == "cns": usecnsmap = 1
        elif mapformat == "omit": useomitmap = 1
        elif mapformat == "2fofc": useccp4map = 1
        else : print "Unrecognised mapformat ,exiting" ; import sys; sys.exit()
            
        if useccp4map == 1 or useomitmap   == 1 :
            fofcmap = "%s%dFo-%dFC.map"%(mtzfn,2,1) ; fcmap = "%sFC.map"%mtzfn
            mapfn = fofcmap
            if usefreer == 1 : 

                fftnew(mtzfn, fofcmap, pdbfile,f1label,f2label,phiclabel,sigf1label,2,1,freerlabel)
                fftnew(mtzfn, fcmap, pdbfile ,f1label,f2label,phiclabel,sigf1label,0,1,freerlabel,1)
            else :
                fftnew(mtzfn, fofcmap, pdbfile , f1label,f2label,phiclabel,sigf1label,2,1,None)
                fftnew(mtzfn, fcmap, pdbfile ,f1label,f2label,phiclabel,sigf1label,0,1,None,1)



        if useomitmap == 1 :
            omit_mapfile = "%s2Fo-1FC.omit.map"%(mtzfn)
            mapfn = omit_mapfile 
            if usefreer == 1 : 
                omitsucc = omitmap(mtzfn, omit_mapfile, f1label,f2label,phiclabel,sigf1label,2,1,freerlabel)
            else :
                omitsucc = omitmap(mtzfn, omit_mapfile, f1label,f2label,phiclabel,sigf1label,2,1,None)

            if omitsucc == None :
                useomitmap = 0 ; useccp4map = 1 ; mapfn = fofcmap

        if usecnsmap == 1 :
            fofcmap = mtzfn+"2fofc.map"
            fcmap = mtzfn+"fc.map"
            mapfn =  mtzfn+"2fofc.map"
            
        

    ############### SET UP   ######################

    if (mapfn !=  None or mtzfn !=None ):
        esmin, esmax, esmean, rcmult = minXSig, maxXSig, .0, 5 
        mapcoeff = "2F1-1F2"
        xrayRestGen.append( prepareChainV5.XrayRestraintsGenerator(mapfn, "map", f2label, phiclabel, mapcoeff, esmin, esmax, esmean, [], edOpt ) )
    



    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    #if mconly != 1 :
    #    scMissInds = incompleteSCcorrection(res, resns, pts)
    #for i in range(len(res.keys())):
    #    incompleteMCcorrection(res, resns, pts , i) 
    #if mconly == 1 : res, pts = removeSC(res, pts, res.keys())


    chainsFound  = getChains(pdbfile)    
    if len (chid12) < 2 :
        print "Atleast two chains define an interface, check value of --chids"
        print "By default interface across all chains will be sampled"
        chid12 = []
        for x in chainsFound:
            if x not in chid12:
                chid12.append(x)

        if len(chid12) < 2:
            print "Only one chain detected in PDB file"
            print "Atleast two chains define an interface, check PDB file"
            import sys ; sys.exit()
        else:
            chi = ""
            for c in chid12 :
                chi = chi+c
            chid12 = chi

    for cc in chid12 : 
        if cc not in chainsFound : 
            print "Chain id requested [ %s ] not found in input PDB file :: "%cc 
            import sys ; sys.exit()


    segs = [] 
    # which interface loops need rebuilding?
    for chi1 in chid12 :
        for chi2 in chid12 :
            if chi1 != chi2 :
                start1, end1 = findStartEnd(chids, chi1)
                start2, end2 = findStartEnd(chids, chi2)
                if start1 == None or end1 == None or start2 == None or end2  == None :
                    print "Error in reading chains..exiting now.."
                    import sys ; sys.exit()

                close1 = findCloseCA(start1,end1, start2,end2,closeCut, res,pts)
                close2 = findCloseCA(start2,end2, start1,end1,closeCut, res,pts)
                segs = segs + makeCloseSegs(close1)
                segs = segs + makeCloseSegs(close2)

    badresids = []
    print "Residues in Interface -----------------------"
    for startindex,endindex in segs :
        for rr in range(startindex, endindex+1) : 
            badresids.append(resids[rr])
            print resids[rr]
    print "------------------------------------------------\n\n\n"

    if len(badresids) == 0:
        print "No residues in contact within specified cutoff %f Angstroms"% closeCut
        import sys; sys.exit()


    from loopbuild12 import Multiloop
    prepC  = prepareChainV5.PrepareChain(rotlib)
    
    nb = 0

    multiPrepC = prepareChainV5.PrepareChain(rotlib)
    for  i in range(nmodels):
        nb = 0 ; cp1 = 0

        if nmodels == 1 > 1 :
            outpdb = "pre."+str(i)+"."+rtkmodel
            outpdbsc = str(i)+"."+rtkmodel
        else :
            outpdb = "pre."+rtkmodel
            outpdbsc = rtkmodel

        if (caRad != None and caRad < 0.49) :
            print "Minimum Ca-radii needs to be 0.5 ANGSTROM."
            if nmodels > 1 : 
                print "Copying", pdbfile ,  " to ","pre."+rtkmodel 
                shutil.copyfile(pdbfile, "pre."+rtkmodel)
            else :
                print "Copying", pdbfile ,  " to ","pre."+str(i)+"."+rtkmodel 
                shutil.copyfile(pdbfile, "pre."+str(i)+"."+rtkmodel)
            nb = 1 ; cp1 = 1 
        
        elif ( len(badresids) > 0 and sconly == 0 ):
            print "\nGenerating MODEL %d of %d asked for ..........."%(i+1,nmodels)

            ml = Multiloop(modelIn, badresids, mconly, caRad, scRad, scvdwr, guidedsampling, popsize, backtrack, 1 , outpdb, xrayRestGen, multiPrepC,sconly, cacaCutoff,[],None,natt,{},1,allOpt,edOpt,None,None) 
            ml.ranker = None 

            #if mtzfn != None or mapfn !=None:
            #    ml.ranker = XrayRanker(mapfn, "map", f2label, phiclabel, mapcoeff, esmin, esmax)
            #    ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
            #    ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble = None

            if mtzfn != None or mapfn !=None:
                ml.cellsym = [ acell, bcell, ccell, alpha, beta, gamma, sgtable[sg][0] ]
                


            nb,built,unbuilt = ml.run()

            print "\n\n==================================== BUILD REPORT ====================================\n\n"
            for ub in unbuilt[0]:
                for bands in ub :
                    if bands != [] :
                        print 
                        print "BAND : [%s] to [%s] could not be built (%s) "%(bands[0],bands[1],bands[2])
                        print 

            for ub in built[0]:
                for bands in ub :
                    if bands != [] :
                        print 
                        print "Successfully generated coordinates for BAND : [%s] to [%s] (%s) "%(bands[0],bands[1],bands[2])                
                        print 



            print "\n\n==================================== BUILD REPORT ENDS ================================\n\n"
            removeMODEL(outpdb)
            
                
        else :
            cp1 = 1
            shutil.copy(pdbfile, outpdb)

        from copyheader import copyheader




        from copyheader import copyheader
        if  sconly == 1 and opsax == 0 and mapfn != None :
            print "** In order to rebuild sidechains set --opsax True ** "
            if cp1 != 1 : 
                removeMODEL(outpdb) ## automatic water addition can go here
                copyheader(outpdb,"modelinit.pdb")
                copyNonprotein(pdbfile,outpdb)
                adjustBfac2(outpdb, pdbfile,nativeBfac, mcBfac, scBfac)
            import sys ; sys.exit()


        if  sconly == 1 and opsax == 1 and mapfn == None :
            print "** In order to rebuild sidechains using opsax MTZ/MAP file required ** "
            if cp1 != 1 : 
                removeMODEL(outpdb) ## automatic water addition can go here
                copyheader(outpdb,"modelinit.pdb")
                copyNonprotein(pdbfile,outpdb)
                adjustBfac2(outpdb, pdbfile,nativeBfac, mcBfac, scBfac)
            import sys ; sys.exit()

        if sconly == 1 and opsax == 0 and mapfn == None :
            print "** In order to rebuild sidechains for given backbone :  "
            print "** (1) Give MTZ/MAP file  ** "
            print "                 or          "
            print "** (2) Generate all atom model using --ca-restraint-radius 0.5 --use-ca-restraints True **  "

            if cp1 != 1 :
                removeMODEL(outpdb) ## automatic water addition can go here
                copyheader(outpdb,"modelinit.pdb")
                copyNonprotein(pdbfile,outpdb)
                adjustBfac2(outpdb, pdbfile,nativeBfac, mcBfac, scBfac)
                
            import sys ; sys.exit()            



        if len(badresids) > 250:
            print "Length of section > 250 residues, can not perform OPSAX"
            opsaxoff = 1 

        if sconly == 1 and opsax == 1 and mapfn != None and opsaxoff == 1 :
            print "Cannot use opsax for proteins > 200 aminoacids, try  rebuilding all atoms :"
            print "--sconly False --ca-restraint-radius 0.5 --use-ca-restraints True "
            if cp1 != 1 :
                removeMODEL(outpdb) ## automatic water addition can go here
                copyheader(outpdb,"modelinit.pdb")
                copyNonprotein(pdbfile,outpdb)
                adjustBfac2(outpdb, pdbfile,nativeBfac, mcBfac, scBfac)
                
            import sys ; sys.exit()

            
            
        if opsax == 0 or mapfn == None or mconly == 1 or opsaxoff == 1:
            shutil.copyfile(outpdb, outpdbsc  )
            print "Renaming", outpdb, " to ", outpdbsc

            
            
        elif  nb != 1 and sconly !=1 : 
            shutil.copyfile(outpdb, outpdbsc  )
            print  "Renaming", outpdb, " to ", outpdbsc




        else :
        
            print "\n\nStarting OPSAX ........"
            scPrepC, useGivenRot, useDEE = prepareChainV5.PrepareChain("PRL"), 0 , 1
            
            
            nb2 = SCplacement(outpdb, scReduction, outpdbsc, "dotfile", useDEE, mapfn, f1label, f2label , phiclabel, mapcoeff, esmin, esmax, addsc, useGivenRot, badresids, scPrepC).run()
            if nb2 != 1:
                shutil.copyfile(outpdb , outpdbsc)


        adjustBfac2(outpdb, pdbfile,nativeBfac, mcBfac, scBfac)
        adjustBfac2(outpdbsc, pdbfile,nativeBfac, mcBfac, scBfac)

        removeMODEL(outpdbsc) ## automatic water addition can go here
        removeMODEL(outpdb) ## automatic water addition can go here
        
        copyNonprotein(pdbfile,outpdbsc) ## automatic water addition can go here
        copyNonprotein(pdbfile,outpdb)

        copyheader(outpdb,"modelinit.pdb")
        copyheader(outpdbsc,"modelinit.pdb")
        




    ### orig starts ##

    #if sconly != 1 :
    #    ml = Multiloop(pdbfile, badresids, mconly, caRad, scRad, scReduction, None, popsize,backtrack, nmodels, "pre"+outpdb, restrGen , prepC)
    #    ml.ranker = None 
    #    ml.cacaCutoff = cacaCutoff
    #    ml.cellsym = cellsym
    #    nb = ml.run()
    #if opsax == 0 or mapfn == None or mconly == 1 or len(badresids) > 300 :
    #    shutil.copyfile("pre%s"%outpdb, outpdb)
    #else :
    #    if  options.sconly == 1 or nb !=1 :
    #        scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 0, 1
    #    else :
    #        scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 1, 1
    #    nb2  = 0
    #    nb2 = SCplacement("pre"+outpdb, scReduction, outfile, "dotfile", useDEE, mapfn, None, None , None, None, esmin, esmax, None, useGivenRot, badresids, scPrepC).run()
    #    if nb2 != 1:
    #        shutil.copyfile("pre"+outpdb, outpdb)
    #        
    #    adjustBfac2(outpdb, pdbfile ,nativeBfac, mcBfac, scBfac)
    #    copyNonprotein(pdbfile, outpdb)
    #    removeMODEL(outpdb) ## automatic water addition can go here
    #    from copyheader import copyheader
    #    copyheader(outpdb,pdbfile)


def callmain() :

    from commonOptions import makeParser , parseOptions
    import optparse ; parser = optparse.OptionParser()
    #    parser = makeParser()
    import sys

    parser.add_option("--xyzin", action='store', type='string', dest='pdbfile', help='Name of input coordinates file')
    parser.add_option("--xyzout", action='store', type='string', dest='pdbout', help='Name of output coordinates file',default="modelout.pdb")
    parser.add_option("--dir-xyzout", action='store', type='string', dest='dir_xyzout', help='Name of directory in which the output files will be created. Complete PATH needed')
    parser.add_option("--mapin", action='store', type='string', dest='mapfn', help='Name of input MAP file (complete PATH needed)',default=None)
    parser.add_option("--hklin", action='store', type='string', dest='mtzfn', help='Name of input (phased) MTZ file (complete PATH needed)',default=None)


    ### Assignment of restraints and radii #############################################
    
    parser.add_option("--use-ca-restraints", action='store', dest='caRes', type= 'string', help='[True/False], Apply positional restraints on the C-alpha atoms',default="True")
    parser.add_option("--use-sc-restraints", action='store', dest='scRes',type= 'string', help='[True/False],  Apply positional restraints on the centroid of the sidechain atoms',default="True",)
    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='Radius of spherical restraint ( Angstrom ) on the C-alpha atom position', default=1.0)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='Radius of spherical restraint ( Angstrom ) on the centroid of the sidechain atoms', default=2.0)

    ### Interface definition #########################################################
    
    parser.add_option("--chainids", action='store', type='string', dest='chids', help='chains to be traced', default=' ')
    parser.add_option("--close-cutoff", action='store', type='float', dest='closeCut', help='how close is close for CAs across interface', default=5)

    ############ General parameters #####################################################
    
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='Factor by which to reduce effective Van der Waals distance for sidechain atoms', default=0.75)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='Population size', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='Level of printed output [1-10], 1:Concise output log, 10: Detailed output log', default=6)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='[numsteps]X[stepsize] e.g. when set to 4X5, In case of failure to build at a residue the program will  backtrack 4 times by 5 residues. For detailed help see Rappertk wiki/webpage', default="4X5")
    parser.add_option("--rotamerlib", action='store', type='string', dest='rotLib', help='[PRL/SCL1.0/SCL0.5/SCL0.2] Name of rotamer library to use when building side chains ', default='PRL')        
    parser.add_option("--num-models", action='store', type='int', dest='nmodels', help='Number of models wanted ', default=1)


    #################### Build parameters ################################################
    
    parser.add_option("--mconly", action='store', type='string', dest='mconly', help='[True/False] Build mainchain only', default="False")
    parser.add_option("--sconly", action='store', type='string', dest='sconly', help='[True/False] Build side chains only, can only be used when MAP/MTZ file is given. See web page for further details', default="False")

    parser.add_option("--opsax", action='store', type='string', dest='opsax', help='[True/False] Reassign side chains with OPSAX, will only be used when MTZ or MAP file is given', default="True")

    parser.add_option("--attempts", action='store', type='int', dest='natt', help='Number of attempts made to build section', default=5)

    parser.add_option("--cacaCutoff", action='store', type='float', dest='cacaCutoff', help='Minimum distance ( angstrom ) between adjacent Calpha atoms in order to detect a chain-break', default=5.)


    #################    Electron density parameters ####################################
    
    parser.add_option("--a", action='store', type='float', dest='a', help='Cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='Cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='Cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='Cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='Cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='Cell spacegroup, in CCP4 notation')
    parser.add_option("--resolution", action='store', type='float', dest='resolution', help='Resolution of the data')
    
    parser.add_option("--FP", action='store', type='string', dest='f1label', help='Column label for FP in MTZ file', default="FP")
    parser.add_option("--SIGFP", action='store', type='string', dest='sigf1label', help='Column label for sigFP in MTZ file', default="SIGFP")
    parser.add_option("--FC", action='store', type='string', dest='f2label', help='Column label for FC in MTZ file', default=None)
    parser.add_option("--PHI", action='store', type='string', dest='phiclabel', help='Column label for PHI in MTZ file', default=None)
    parser.add_option("--use-FREER", action='store', type='string', dest='usefreer', help='[True/False] Use FreeR set ? ', default="False")
    parser.add_option("--FREER", action='store', type='string', dest='freerlabel', help='Column label for FreeR in MTZ file', default=None)




    parser.add_option("--models-get-native-bfactors", action='store', type='string', dest='nativeBfac', help='[True/False] Assign B-factors of remodelled atoms to original values', default="False")
    parser.add_option("--default-mainchain-b-factor", action='store', type='float', dest='mcBfac', help='The value of B-factor assigned to the newly built main chain atoms', default=20.)
    parser.add_option("--default-sidechain-b-factor", action='store', type='float', dest='scBfac', help='The value of B-factor assigned to the newly built side chain atoms', default=30.)



    ### Electron density parametets #########################################

    parser.add_option("--minimum-sig", action='store', type='float', dest='minXSig', help='Minimum sigma ', default=0.25)
    parser.add_option("--maximum-sig", action='store', type='float', dest='maxXSig', help='Maximum sigma ', default=2.0)



    ########## Optional restraints ##########################################
    
    parser.add_option("--make-ed-optional", action='store', type='string', dest='edOpt', help='[True/False]  If False, then the mainchain will be unconditionally forced to lie in positive density. If True then positive density restraint on the mainchain will be made optional.This is useful when tracing through a structure with regions in very poor (non-existent) density', default= "False")

    
    parser.add_option("--make-all-restraints-optional", action='store', type='string', dest='allOpt', help='[True / False ]  If True, then all  restraints will be made optional', default="False")    
    parser.add_option("--add-sidechains", action='store', type='string', dest='addsc', help='Build missing side chains [True/False] Note: If building all atoms then missing sidechain atoms in the region to be built will be added and if building only sidechains then missing sidechain atoms can be optionally added', default='False')

    parser.add_option("--use-given-rotamer", action='store', type='string', dest='userot', help='Use given rotamer [True/False]', default='False')

    parser.add_option("--maptype", action='store', type='string', dest='mapformat', help='Options are : [2fofc/omit] ', default="2fofc")



    options = parseOptions(parser)

    from checkfornone import checkfornone_ppi
    
    options.modelN2C  = 'False'
    pdbfile,pdbout,dir_xyzout,mapfn,mtzfn,caRes,scRes,caRad,scRad,scReduction,popsize,verbose,backtrack,rotLib,nmodels,mconly,sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freeRlabel,poordum,ccdum,loopdum,startdum,startcdum,stopdum,stopcdum,chaindum,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt,closuredum,addscdum,userotdum =    checkfornone_ppi(options.pdbfile,options.pdbout,options.dir_xyzout,options.mapfn,options.mtzfn,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,options.nmodels,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freerlabel,"False",0.9,None,None,None,None,None,None,options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt,"True","False","False")


    print "	--close-cutoff ".ljust(30), 
    print "%10f "%options.closeCut

    print "	--chids ".ljust(30), 
    print "%10s "%options.chids
       
    main(options.pdbfile,options.pdbout,options.dir_xyzout,options.mapfn,options.mtzfn,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,options.nmodels,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freerlabel,options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt,options.mapformat , options.chids, options.closeCut) 
    


if __name__ == "__main__" : callmain()
