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

def main(pdbfile, chid12, mconly, sconly, closeCut, caRad, scRad, scReduction, restrGen, outpdb, popsize, backtrack, nmodels  , sclib ,cellsym , cacaCutoff, opsax , mapfn , nativeBfac, mcBfac , scBfac , verbose) :

    misc.setVerbosity(verbose) ; randomize(1)


    
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    if mconly != 1 :
        scMissInds = incompleteSCcorrection(res, resns, pts)

    for i in range(len(res.keys())):
        incompleteMCcorrection(res, resns, pts , i) 


        
    if mconly == 1 : res, pts = removeSC(res, pts, res.keys())
    
    
    print "Residues and atom-numbers-----------------------"
    keys = res.keys() ; keys.sort()
    for k in keys :
        print resids[k] 
    print "------------------------------------------------\n\n\n"

    if len (chid12) < 2 :
        print "Atleast two chains define an interface, check value of --chids"
        import sys ; 
        #sys.exit()
        print "By default interface across all chains will be sampled"
        chid12 = []
        for x in chids.values():
            if x not in chid12:
                chid12.append(x)
        if len(chid12) < 2:
            print "Only one chain detected in PDB file"
            print "Atleast two chains define an interface, check PDB file"
            sys.exit()
        else:
            chi = ""
            for c in chid12 :
                chi = chi+c
            chid12 = chi

    for cc in chid12 : 
        if cc not in chids.values():
            print "Chain id requested [ %s ] not found in input PDB file :: "%cc 
            import sys ; 
            sys.exit()

    



    segs = [] 
    # which interface loops need rebuilding?
    for chi1 in chid12 :
        for chi2 in chid12 :
            if chi1 != chi2 :
                start1, end1 = findStartEnd(chids, chi1)
                start2, end2 = findStartEnd(chids, chi2)
                assert start1 and end1 and start2 and end2
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

    import prepareChainV3
    prepC  = prepareChainV3.PrepareChain(sclib)
    
    if len(badresids) == 0:
        print "No residues in contact within specified cutoff %f Angstroms"% closeCut
        import sys; sys.exit()

    from loopbuildV5 import Multiloop
    nb = 0
    if sconly != 1 :
        ml = Multiloop(pdbfile, badresids, mconly, caRad, scRad, scReduction, None, popsize,backtrack, nmodels, "pre"+outpdb, restrGen , prepC)
        ml.ranker = None 
        ml.cacaCutoff = cacaCutoff
        ml.cellsym = cellsym
        nb = ml.run()


    if opsax == 0 or mapfn == None or mconly == 1 or len(badresids) > 300 :
        shutil.copyfile("pre%s"%outpdb, outpdb)
    else :
        if  options.sconly == 1 or nb !=1 :
            scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 0, 1
        else :
            scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 1, 1

        nb2  = 0
        nb2 = SCplacement("pre"+outpdb, scReduction, outfile, "dotfile", useDEE, mapfn, None, None , None, None, esmin, esmax, None, useGivenRot, badresids, scPrepC).run()
        if nb2 != 1:
            shutil.copyfile("pre"+outpdb, outpdb)
            
        adjustBfac2(outpdb, pdbfile ,nativeBfac, mcBfac, scBfac)
        copyNonprotein(pdbfile, outpdb)
        removeMODEL(outpdb) ## automatic water addition can go here
        from copyheader import copyheader
        copyheader(outpdb,pdbfile)


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
    
    parser.add_option("--chids", action='store', type='string', dest='chids', help='chains to be traced', default=' ')
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
    parser.add_option("--PHIC", action='store', type='string', dest='phiclabel', help='Column label for PHIC in MTZ file', default=None)
    parser.add_option("--use-FREER", action='store', type='string', dest='usefreer', help='[True/False] Use FreeR set ? ', default="False")
    parser.add_option("--FREER", action='store', type='string', dest='freeRlabel', help='Column label for FreeR in MTZ file', default=None)




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

    from checkfornone import checkfornone
    
    options.modelN2C  = 'False'
    

    options.pdbfile,options.pdbout,options.dir_xyzout,options.mapfn,options.mtzfn,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,options.nmodels,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freeRlabel,ndum,poordum,ccdum,loopdum,startdum,startcdum,stopdum,stopcdum,chaindum,n2cdum,options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt,closuredum,addscdum,options.userot,options.mapformat = checkfornone(options.pdbfile,options.pdbout,options.dir_xyzout,options.mapfn,options.mtzfn,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,options.nmodels,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freeRlabel,2,"False",0.9,None,None,None,None,None,None,"False",options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt,"True","False",options.userot,options.mapformat)

    






    if not os.path.isdir(options.dir_xyzout) : os.mkdir(options.dir_xyzout)
    os.chdir(options.dir_xyzout)
    fp = open('parameters.txt', 'w')
    fp.close()

    if (options.mconly == 1 and  options.sconly == 1) :
        print "Mainchain only and sidechain only modes cannot be used together"
        import sys ; sys.exit()
        

    if (os.path.isfile(options.pdbfile)==False) :
        print "Cannot find file %s "%options.pdbfile
        print "No file in directory ", os.getcwd() , options.pdbfile
        import sys ; 
        sys.exit()

    shutil.copyfile(options.pdbfile, "%s/init.pdb" % (options.dir_xyzout))

    cellsym = []
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





        from prepareChain import XrayRestraintsGenerator , PrepareChain
        xrayRestGen = []
        esmin, esmax, esmean = .000, 5., .0, 

        
        #ml.ranker = XrayRanker(ccp4map, "map", options.f2label, options.phiclabel, mapcoeff, esmin, esmax)
        #ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
        #ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble = None
        from data import sgtable , long2shortHM
        cellsym = []
        if options.sg  in long2shortHM.keys():
            shortsg = long2shortHM[options.sg]
            options.sg = shortsg
        if options.sg not in sgtable.keys():
            print "Check --sg , Not recognised [%s][%d]"%( options.sg, len(options.sg))
            print "Check --sg , Not recognised [%s][%d]"%( options.sg, len(options.sg))
            sys.exit()

        cellsym = [ options.a, options.b, options.c, options.alpha, options.beta, options.gamma, sgtable[options.sg][0] ]

    if options.mapfn != None  :
        if (os.path.isfile(options.mapfn)==False) :
            print "Cannot find file %s "%options.mapfn
            print "No file in directory ", os.getcwd()
            import sys ; 
            sys.exit()

        shutil.copyfile(options.mapfn, "%s/init.map" % (options.dir_xyzout))
        xrayRestGen.append( XrayRestraintsGenerator(options.mapfn, "map", options.f2label, options.phiclabel, "", esmin, esmax, esmean,[], options.edOpt ) )
            
    
    elif options.mtzfn != None : 
        if (os.path.isfile(options.mtzfn)==False) :
            print "Cannot find file %s "%options.mtzfn
            print "No file in directory ", os.getcwd()
            import sys ; 
            sys.exit()

        shutil.copyfile(options.mtzfn, "%s/init.mtz" % (options.dir_xyzout))
        mtzfilepath  = "%s/init.mtz" % options.dir_xyzout
        ccp4map = "%s%dFo-%dFC.map"%(mtzfilepath,options.m,options.n)
        mapcoeff = "%dF1-%dF2"%(options.m, options.n)

        if (options.f2label == None or options.phiclabel == None ):
            print "Column labels for  FC and PHIC  not specified, will use input PDB structure to obtain FC and PHIC"
            print options.pdbfile, mtzfilepath, "phased.mtz",options.resolution, options.usefreer, options.f1label,options.sigf1label,options.freerlabel 

            sfall(options.pdbfile, mtzfilepath, "phased.mtz",options.resolution, options.usefreer, options.f1label,options.sigf1label,options.freerlabel)

            mtzfilepath  = "phased.mtz"
            options.f2label  = "FC"
            options.phiclabel = "PHIC"

        if options.usefreer == 1 : 
            if options.freerlabel == None :
                print "If you would like to use the freer set, please specify column name for freeR"
                sys.exit()
            fft(mtzfilepath, ccp4map, options.f1label,options.f2label,options.phiclabel,options.sigf1label,options.m,options.n,options.freerlabel)
        else :
            fft(mtzfilepath, ccp4map, options.f1label,options.f2label,options.phiclabel,options.sigf1label,options.m,options.n,None)

            
        xrayRestGen.append( XrayRestraintsGenerator(ccp4map, "map", options.f2label, options.phiclabel, mapcoeff, esmin, esmax, esmean,[], options.edOpt ) )
        options.mapfn = ccp4map

    else :
        print "WARNING !! No input mtz or map file , No electron density restraints will be used"

    
    if options.detectBreaks == 1:
        caca = options.cacaCutoff
    else :
        caca = None

    main("init.pdb", options.chids, options.mconly, options.sconly, options.closeCut,
         options.caRad, options.scRad, options.scReduction, xrayRestGen,
         options.pdbout, options.popsize, options.backtrack, options.nmodels , options.rotLib,cellsym , caca, options.opsax , options.mapfn, options.nativeBfac, options.mcBfac, options.scBfac,options.verbose   )

if __name__ == "__main__" : callmain()
