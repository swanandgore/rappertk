
from xray import sfcheck
import os, shutil, re , sys , string
import geometry
from xray import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS, fft, omitmap, mapman
from procrun import proc_run_exitOnError as execCmd
from xcheck import XrayScorer, XrayRanker
from data import sgtable , long2shortHM
from evalCAtrace import comparePhiPsiOmegaChi
from pdbr import protein, isAAres  , line2atomname
import prot2res
from scplacement import SCplacement
from loopbuild12 import Multiloop, incompleteSCcorrection2, locateRegionsRandomize3


import prepareChainV5
from copyheader import copyheader
from multProtref import joinPDBs, splitPDBs, restoreChainids
import checkProtChains
import checkProtChainsV4
from peptidebuild import ModelRenderer
import misc
from pdbr import isPdbAAline, isPdbAtomLine, line2resn , line2resnum , line2chid , line2atomnumber
from pdbr import makeResid

## Last modified 18-02-2010 :  Using fft for mtz to ccp4map
### Todo : If band building fails, eliminate OPSAX for those residues if present in missing/ missing sidechains




ccp4args = {
    0: [{"reftype":"restrained",    "wa":0.20,    "breftype":"ISOTROPIC",   "ncyc":20}], #on native

    1: [{"reftype":"unrestrained",  "wa":0.75,    "breftype":"OVER",        "ncyc":20, "assignBfac":[20,30]}, #on catrace
        {"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":40}], #on catrace

    2: [{"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[60,90]}],
    3: [{"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[60,90]}],
    4: [{"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[60,90]}],
    5: [{"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[60,90]}],
    6: [{"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[60,90]}],

    7: [{"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[15,25]}],
    8: [{"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[15,25]}],
    9: [{"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[15,25]}],
   10: [{"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[15,25]}],
   11: [{"reftype":"restrained",    "wa":0.75,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[15,25]}],
   12: [{"reftype":"restrained",    "wa":0.50,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[15,25]}],
   13: [{"reftype":"restrained",    "wa":0.20,    "breftype":"ISOTROPIC",   "ncyc":20, "assignBfac":[15,25]}],
   14: [{"reftype":"restrained",    "wa":0.20,    "breftype":"ISOTROPIC",   "ncyc":40, "assignBfac":[ 5, 6]}],
}




def getSCcentroids(filename) :
    from pdbr import protein, isAAres
    import prot2res
    resid2scen = {}
    
    
    prot = protein(filename, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    centroids = {}
    for ri in res.keys() :
        if not isAAres(resns[ri]) : continue
        if resns[ri] in ["GLY","ALA"] : continue
        cen = [0.,0.,0.] ; cnt = 0

        for an,ai in res[ri].items() :
            if an in [" N  "," CA "," C  "," O  "," CB "] : continue
            for xi in range(3) :
                cen[xi] += pts[ai][xi]
            cnt += 1
        if cnt > 0 :
            for xi in range(3) :
                cen[xi] /= cnt
            resid2scen[resids[ri]] = cen
    return resid2scen

## if a sidechain is changed, the centroid of sidechain will move significantly, at least 0.1 A
def findChangedSC(origfile, newfile) :
    from math import sqrt
    changed = []
    origcen = getSCcentroids(origfile)
    newcen = getSCcentroids(newfile)
    for k,v in origcen.items() :
        if k in newcen.keys():
            dist = (newcen[k][0]-v[0])*(newcen[k][0]-v[0])+(newcen[k][1]-v[1])*(newcen[k][1]-v[1])+(newcen[k][2]-v[2])*(newcen[k][2]-v[2])
            if sqrt(dist) > 0.1 : changed.append(k)
    return changed




def findAllSC(origfile) :


    from math import sqrt
    changed = []
    origcen = getSCcentroids(origfile)
    for k,v in origcen.items() :
        changed.append(k)
    return changed


## change ILE:CD to ILE:CD1 if necessary
# OT1 to O
def parseLoopres(pdbfile,loopres):
    loopresids = []
    
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    if (os.path.isfile(loopres)==False) :
        print "Cannot find file %s "%loopres
        print "No file in directory ", os.getcwd()
        import sys ; sys.exit()
    lines = open(loopres, 'r').readlines()
    import re
    for li in range(len(lines)) :
        l = re.sub("\n", "", lines[li])
        if l == '' or l[0] == '#' : continue
        if l[0:5] == "loop " :
            rids,chain = parseResnums(l)
            if chain not in chids.values():
                print "Chain id in loopres does not match any in Input PDB file" , chain 
                import sys ; sys.exit()
            for i in range(int(rids[0]),int(rids[1])+1) :
                for k , v in resnums.items():
                    if int(v) == i and chids[k] == chain :
                        loopresids.append(resids[k])
    return loopresids


def parseResnums(l) :
    start, resnums = None, []
    for i in range(len(l)) :
        if l[i] == '[' : start = i
        elif l[i] == ']' :
            resnums.append( l[start+1:i] )
    for i in range(len(l)) :
        if l[i] == '\'' and l[i+2] == '\'' :
            chain = l[i+1]
        
    return resnums,chain







## assert that there is only 1 model, and remove model, endmdl line
def removeMODEL(pdbfile) :
    newlines, nm = [], 0
    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        import sys ; sys.exit()

    for l in open(pdbfile,'r').readlines() :
        if l[0:5] == 'MODEL' : nm += 1
        elif l[0:6] == 'ENDMDL' : continue
        else : newlines.append(l)
    #assert nm == 1
    op = open(pdbfile, 'w')
    for l in newlines : op.write(l)
    op.close()

def removeZeroLines(filename) :
    if (os.path.isfile(filename)==False) :
        print "Cannot find file %s "%filename
        print "No file in directory ", os.getcwd()
        import sys ; sys.exit()

        

    lines = open(filename, 'r').readlines()
    fp = open(filename, 'w')
    for l in lines :
        if re.compile('ATOM.* 0.000.* 0.000.* 0.000').search(l) : continue
        if re.compile('HETATM.* 0.000.* 0.000.* 0.000').search(l) : continue
        if re.compile('HETATM.*9999.*9999.*9999').search(l) : continue
        if re.compile('ATOM.*9999.*9999.*9999').search(l) : continue
        fp.write(l)
    fp.close()

def moleman(pdbfilename) :
    input = ["READ"]
    input.append(pdbfilename)
    input.append("AUTO")
    input.append("write")
    input.append(pdbfilename+".moleman")
    for i in range(9) : input.append("")
    input.append("QUIT")


    execCmd("lx_moleman", input)
    os.rename(pdbfilename+".moleman", pdbfilename)

## verify that last line in topdb in END, change that to TER
## then copy all non-protein lines
def copyNonprotein(frompdb, topdb, waterAlso=1) :
    from pdbr import isPdbAAline, isPdbAtomLine, changeSegid
    if (os.path.isfile(topdb)==False) :
        print "Cannot find file %s "%topdb
        print "No file in directory ", os.getcwd()
        import sys ;        sys.exit()
    
    lines = open(topdb, 'r').readlines()
    #assert lines[ len(lines)-1 ][0:3] == "END"

    print frompdb , topdb , waterAlso, len(lines)
    if lines[ len(lines)-1 ] in [ "END\n", "END " ] : lines[ len(lines)-1 ] = ""
    #else : lines.append("TER\n")
    print "copying nonprotein atoms from", frompdb, "to", topdb


    if (os.path.isfile(frompdb)==False) :
        print "Cannot find file %s "%frompdb
        print "No file in directory ", os.getcwd()
        import sys ;        sys.exit()

    for l in open(frompdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if isPdbAAline(l) : continue
        if not waterAlso and line2resn(l) == "HOH" : continue
        lines.append( changeSegid(l,"    ") )
    op = open(topdb, 'w')
    for l in lines : op.write(l)
    print >> op, "END"
    op.close()


def copyWater(frompdb, topdb, waterAlso=1) :
    from pdbr import isPdbAAline, isPdbAtomLine, changeSegid
    if (os.path.isfile(topdb)==False) :
        print "Cannot find file %s "%topdb
        print "No file in directory ", os.getcwd()
        import sys ; sys.exit()
    
    lines = open(topdb, 'r').readlines()
    #assert lines[ len(lines)-1 ][0:3] == "END"
    if lines[ len(lines)-1 ] in [ "END\n", "END " ] : lines[ len(lines)-1 ] = ""
    #else : lines.append("TER\n")
    print "copying HOH atoms from", frompdb, "to", topdb


    if (os.path.isfile(frompdb)==False) :
        print "Cannot find file %s "%frompdb
        print "No file in directory ", os.getcwd()
        import sys ;        sys.exit()

    for l in open(frompdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if isPdbAAline(l) : continue
        if line2resn(l) == "HOH" : 
            lines.append( changeSegid(l,"    ") )
    op = open(topdb, 'w')
    for l in lines : op.write(l)
    print >> op, "END"
    op.close()    

def adjustBfac(pdbfile, refpdb) :
    print "\nCopy Bfactors from", refpdb, "to", pdbfile, "if coordinates are same, else set to 30"
    from pdbr import line2bfac, isPdbAtomLine, line2bfac, line2atomid, changeBfactor, line2crdstr
    id2bo = {} # read in crd, bfac, for each atom-id

    
    if (os.path.isfile(refpdb)==False) :
        print "Cannot find file %s "%refpdb
        print "No file in directory ", os.getcwd()
        import sys ;        sys.exit()


    for l in open(refpdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if id2bo.has_key(line2atomid(l)) : continue
        id2bo[line2atomid(l)] = [line2bfac(l), line2crdstr(l)]
    newlines = []


    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        import sys ;        sys.exit()



    for l in open(pdbfile, 'r').readlines() :
        if not isPdbAtomLine(l) : newlines.append(l) ; continue
        bfac, crdstr = id2bo[ line2atomid(l) ]
        if crdstr == line2crdstr(l) : l = changeBfactor(l, bfac)
        else : l = changeBfactor(l, 30.0)
        newlines.append(l)
    op = open(pdbfile, 'w')
    for l in newlines : op.write(l)
    op.close()


def adjustBfac2(pdbfile, refpdb,nativeBfac, mcBfac, scBfac) :
    
    from pdbr import line2bfac, isPdbAtomLine, line2bfac, line2atomid, changeBfactor, line2crdstr
    id2bo = {} # read in crd, bfac, for each atom-id

    
    if (os.path.isfile(refpdb)==False) :
        print "Cannot find file %s "%refpdb
        print "No file in directory ", os.getcwd()
        import sys ;        sys.exit()

    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        import sys ;        sys.exit()



    for l in open(refpdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if id2bo.has_key(line2atomid(l)) : continue

        id2bo[line2atomid(l)] = [line2bfac(l), line2crdstr(l)]
    newlines = []



    if nativeBfac == 1:
        print "\nCopy Bfactors from", refpdb, "to", pdbfile
        for l in open(pdbfile, 'r').readlines() :
            if not isPdbAtomLine(l) :
                newlines.append(l) ; continue
            if line2atomid(l)  in id2bo.keys() :
                bfac, crdstr = id2bo[ line2atomid(l) ]
                l = changeBfactor(l, bfac)
            else :
                if line2atomname(l) in [' N  ',' CA ',' C  ',' O  ']  : 
                    l = changeBfactor(l, mcBfac)
                else :
                    l = changeBfactor(l, scBfac)

            newlines.append(l)

    else :
        for l in open(pdbfile, 'r').readlines() :
            if not isPdbAtomLine(l) : newlines.append(l) ; continue
            if line2atomid(l)  in id2bo.keys() :
                bfac, crdstr = id2bo[ line2atomid(l) ]
                if crdstr == line2crdstr(l) :
                    l = changeBfactor(l, bfac)
                else :
                    if line2atomname(l) in [' N  ',' CA ',' C  ',' O  ']  : 

                        l = changeBfactor(l, mcBfac)
                    else :
                        l = changeBfactor(l, scBfac)
            else:
                if line2atomname(l) in [' N  ',' CA ',' C  ',' O  ']  : 
                    
                    l = changeBfactor(l, mcBfac)
                else :
                    l = changeBfactor(l, scBfac)

            newlines.append(l)


    op = open(pdbfile, 'w')
    for l in newlines : op.write(l)
    op.close()


# minimal copying of pdb file. ie copy onlt the coordinate records. only the first instance of a resid + atomid is allowed in.
# Bfactors changed to 20 mcB for mainchain, scB for all others.
# occupancy set to 1.


def findWorseFits(XcorrZ, cutoff=0.9) :
    badkeys = []
    #vals = list(XcorrZ.values())
    #vals.sort()
    #cutoff = vals[ (len(vals)-1)/4 ]
    for k in XcorrZ.keys() :
        if XcorrZ[k] < cutoff and k[0:3] != "HOH" : badkeys.append(k)
    #print "BADKEYS", cutoff, len(badkeys), badkeys
    return badkeys

def changeBfacs(filename, bfacs) :
    from pdbr import isPdbAAline, line2atomname, line2occu, changeBfactorOccupancy
    lines = []

    if (os.path.isfile(filename)==False) :
        print "Cannot find file %s "%filename
        print "No file in directory ", os.getcwd()
        import sys ;        sys.exit()

    for l in open(filename, 'r').readlines() :
        if isPdbAAline(l) :
            if line2atomname(l) in [" N  "," CA "," C  "," O  "] : l = changeBfactorOccupancy(l, bfacs[0], line2occu(l))
            else : l = changeBfactorOccupancy(l, bfacs[1], line2occu(l))
        lines.append(l)
    fp = open(filename, 'w')
    for l in lines : fp.write(l)
    fp.close()

def refmacRefinement(mtzin,pdbin, mtzout,pdbout, a,b,c,alpha,beta,gamma,sg,reso, ccp4Args,cycle) :
    from xray import refmac
    if os.path.isfile(pdbout) : return
    inpdb, bufpdb = pdbin, "buf.pdb"
    for ccpa in ccp4Args[cycle] :
        if "assignBfac" in ccpa.keys() : changeBfacs(inpdb, ccpa[assignBfac])
        refmac(mtzin, mtzout, inpdb, bufpdb, reso, ccpa["reftype"], ccpa["wa"], ccpa["breftype"], ccpa["ncyc"])
        inpdb = bufpdb ; bufpdb = bufpdb + "1"
    print "Renaming", inpdb, pdbout
    os.rename(inpdb, pdbout)



def getChains(pdb):
    allchains = []

    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    
    for k,v in chids.items():
        if v not in allchains:
            allchains.append(v)
    
    return allchains





def main(pdbfile,         pdbout,         dir_xyzout,         mapfn,         mtzfn,         caRes,         scRes,         caRad,         scRad,         scReduction,         popsize,         verbose,         backtrack,         rotLib,         nmodels,         mconly,         sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freeRlabel,n,poorOnly,poorThreshold,loopres,start,startcode,stop,stopcode,chainid,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt,closure="True",addsc="False",userot="True",mapformat="2fofc",snappy="False"):    




    
    import sys
    from checkfornone import fixdirnames

    
    if poorOnly == "True" and mapfn != None and mtzfn == None:
        print "If you want to rebuild poor fitting regions please provide structure factor (.mtz) file instead of the density map"
        sys.exit()
    if poorOnly == "True" and mapfn == None and mtzfn == None:
        print "If you want to rebuild poor fitting regions please provide structure factor (.mtz) "
        sys.exit()        

                                                                       


    ######### Check for / make out directory #############

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
    poorOnly =  str2bool[poorOnly]
    nativeBfac =  str2bool[nativeBfac]
    modelN2C =  "False"
    caRes =  str2bool[caRes]
    scRes =  str2bool[scRes]
    edOpt =  str2bool[edOpt]
    allOpt =  str2bool[allOpt]
    closure = str2bool[closure]
    addsc = str2bool[addsc]
    userot = str2bool[userot]
    
    if mconly == 0 :        mconly = None
    if addsc == 0 :         addsc = None
    if userot  == 0 :       userot = None
    if allOpt == 1 :        edOpt = 1
    if caRes == 0 :         scRes = 0
    
    
    ########### CHECK FOR AVAIALABLE CELL PARAMETERS, SPACE GROUP and RESOLUTION

    if (mapfn !=  None or mtzfn !=None ):
        from stump import getCRYST , getRESO

        if (a == None or b == None or c == None or alpha == None or beta == None or gamma == None) :
            print "Getting cell paramters from coordinate file....."
            a,b,c,alpha , beta , gamma,d1  = getCRYST(pdbfile)
            if (a == None or b == None or c == None or alpha == None or beta == None or gamma == None ):
                print "CRYST card cannot be read from coordinate file. Please input cell paramater a, b , c , alpha, beta , gamma = ",a , b , c , alpha , beta  , gamma 
                import sys ; sys.exit()
            print "Found a b c alpha beta gamma  ", a , b , c , alpha , beta  , gamma 

        if sg == None : 
            print "Getting space group from coordinate file....."
            d1,d2,d3,d4 , d5 , d6, sg  = getCRYST(pdbfile)
            if sg == None : 
                print "Please input space group " , sg ; import sys ; sys.exit()
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

    if (caRes == 0 and (start == None or stop == None) and sconly !=1):
        print "Need to set --use-ca-restraints True or specify [--loopseq --chainid --start --stop ] if you require loop building"
       	import sys ; sys.exit(0)

        
    ## somesensechecking

    chainsFound  = getChains(pdbfile)    
    if chainid != None and len(chainid) != 1 :
        print "Only one chain can be specified at a time" ;     import sys ; sys.exit(1)


    if  chainid != None and  chainid not in chainsFound : 
        print "Chain id not found in pdb file", chainid , chainsFound ; import sys ; sys.exit()

    if loopres != None :
        if start == None : print "Enter start residue for loop %s to be modelled"%loopres ; import sys ; sys.exit()
        if stop == None : print "Enter stop residue for loop %s to be modelled"%loopres ; sys.exit()

        if chainid == None and len(chainsFound) != 1 :
            print "More than one chain detected , enter chain using --chainid " ; import sys ;sys.exit()
            
        elif chainid == None and len(chainsFound) == 1 :
            chainid = chainsFound[0] ; print "Setting chain id to ",chainsFound[0]

    ################## Checks for range of residues #######################################

    ## All chains
    if start == None and stop == None and loopres == None and chainid == None  :
        badresids = None ; 	print "Rebuilding all chains"
	if caRes == 0 and sconly == 0 : print "--use-ca-restraints to be set to True" ; import sys ; sys.exit()
            

    ## Single chains
    if start == None and stop == None  and chainid != None :
        print "Rebuilding chain [%s] "%chainid
        if caRes == 0 and sconly == 0 : print "--use-ca-restraints to be set to True" ; import sys ; sys.exit()
        prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        for k , v in chids.items():
            if v == chainid :
                badresids.append(resids[k])

    ### if only start or stop is given irrespective of chain id then exits: ambiguos input
    if  ((start !=None and stop== None) or (start == None and stop != None))  :
        if start == None : print "Enter start residue"  ;import sys ;            sys.exit()
        if stop == None : print "Enter stop residue"    ;  import sys ;            sys.exit()

    ### If start and stop are given and no chain id then check coordinate file for options
    if chainid == None and len(chainsFound) != 1 and ( start != None and stop !=None ):
        print "More than one chain detected , enter chain using --chainid " ; import sys ;        sys.exit()
            
    if chainid == None and len(chainsFound) == 1  and start != None and stop !=None :
        chainid = chainsFound[0] ; print "Setting chain id to ",chainsFound[0]
        
    
    if  start != None and stop != None and chainid != None  and loopres == None :
        prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        startfound = None ; stopfound = None
        for k,v in resnums.items():
            if start == int(v) and chainid == chids[k] and startcode  == inscodes[k] : startfound = 1
            if stop == int(v) and chainid == chids[k] and stopcode  == inscodes[k] : stopfound = 1
        if startfound == None or stopfound == None :
            print "Residue range %d%s %d%s of chain ['%s'] not found in coordinate file %s"%(start,startcode,stop,stopcode,chainid,pdbfile)
            sys.exit()

    #### Non ab initio loop modelling
    if ( (start != None and stop != None and chainid != None) and caRes == 1):
        loopresids = parseLoopSequenceSC(modelIn,loopres,start,stop,chainid,mconly)
        badresids = loopresids



    ### Get SEQUENCE FROM PDB FILE IF NOT SPECIFIED, (ab-initio mode)
    if start != None and stop !=None and  caRes == 0  and chainid != None and  loopres == None :
        from data import three2one
        if verbose > 7 :
            print "Reading loop seqeunce from coordinates file..."
            print "Loop to be modelled...",start,stop, chainid
        
        prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        loopseq = "" ; loopcounter = start

        loopstartid = None ; loopstopid = None
        for k,v in resnums.items():
            if chids[k] != chainid : continue
            if int(v) == start and inscodes[k] == startcode : loopstartid = k
            if int(v) == stop and inscodes[k] == stopcode :  loopstopid = k 
        if loopstartid == None or loopstopid == None:
            print "Residue range %d%s %d%s of chain ['%s'] not found in coordinate file %s"%(start,startcode,stop,stopcode,chainid,pdbfile)
            import sys ; sys.exit()
            
        for k in range(loopstartid,loopstopid+1):
            if resns[k] not in three2one.keys() :
                print "Unrecognised amino acid type %s for residue %d" %(resns[k],resids[k])
                continue
            loopseq = loopseq + three2one[resns[k]]
        loopres = loopseq
        if verbose > 7 :            print "Loop sequence read from file", loopres

    if    start != None and stop !=None and  caRes == 0  and chainid != None and  loopres != None:
	loopresids,opsaxoff,nullres,missingpts = parseLoopSequence(modelIn,loopres,start,stop,chainid,mconly)
	badresids = loopresids 
        if missingpts != {} :
            modelIn = modelIn+".fixed"




   ################   COPY MAP AND MTZ FILES  #########################################


    if  mapfn == None and mtzfn == None :

        if poorOnly == 1 :
            print "In order to rebuild poor fitting regions --hklin or --mapin  needs to be specified"
            import sys ; sys.exit()
        print "********** WARNING!! no mtz or map file given, no electron density restraints will be used *************"

    if  mapfn != None and mtzfn != None :
        
        print "********** WARNING!! both mtz and map file given, please only give either MTZ file or MAP file  *************"
        import sys ; sys.exit()

    if  mapfn != None:
        if (os.path.isfile("%s/%s" % (dir_xyzout,"init.map"))==True) :
            os.remove("%s/%s" % (dir_xyzout,"init.map"))        
        useomitmap,usecnsmap,useccp4map  = None,None,None
        esmin, esmax, esmean, rcmult, xscoreCutoff = minXSig, maxXSig, .0, 5, poorThreshold
        mapfilepath= mapfn
        if (os.path.isfile(mapfilepath)==False) :
            print "Cannot find file %s " %mapfn ; import sys ;            sys.exit()
        if  "map" not in mapfn : 
            print "Cannot understand i/p map format,the file should be  (*.map) format " ; import sys ; sys.exit()
        from xray import mapdump
        if  mapdump(mapfn)  == None :
            print "Check format of map file, needs to be CCP4 format";
            import sys ;              sys.exit()            


        shutil.copyfile(mapfilepath, "%s/init.map" % (dir_xyzout))
        ccp4map = mapfn

    elif mtzfn != None :
        if (os.path.isfile("%s/%s" % (dir_xyzout,"init.mtz"))==True) :
            os.remove("%s/%s" % (dir_xyzout,"init.mtz"))            
        mtzfilepath= "%s"%(mtzfn)
        esmin, esmax, esmean, rcmult, xscoreCutoff = minXSig, maxXSig, .0, 5, poorThreshold
        
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
            if freeRlabel == None :
                print "If you would like to use the freeR set, please specify column name for freeR"
                import sys ; sys.exit()
            
        
        if  mtzdump(mtzfn,f1label,sigf1label,freeRlabel,f2label,phiclabel)  == None :
            print "Check structure factor file"
            import sys ;              sys.exit()            
            

        if (f2label == None or phiclabel == None ):
            print "Column labels for  FC and PHI  not specified, will use input coordinate structure to obtain FC and PHI"

            
            if sfall(pdbfile, mtzfn, "phased.mtz",resolution,usefreer,f1label,sigf1label,freeRlabel) == None :
                
                print "Structure factor file cannot bephased , please enter column labels for  FC and PHI"
                import sys ; sys.exit()
            
            mtzfn  = "phased.mtz"
            f2label  = "FC" ; phiclabel = "PHIC"



        ### MAKE UNWEIGHTED DIFF MAP : todo sigmaA weighted maps #########################

        from xray import fftnew
        useomitmap,usecnsmap,useccp4map  = None,None,None
        if mapformat == "cns": usecnsmap = 1
        elif mapformat == "omit": useomitmap = 1
        elif mapformat == "2fofc": useccp4map = 1
        else : print "Unrecognised mapformat ,exiting" ; import sys; sys.exit()
            
        if useccp4map == 1 or useomitmap   == 1 :
            fofcmap = "%s%dFo-%dFC.map"%(mtzfn,n,((n)-1))
            fcmap = "%sFC.map"%mtzfn
            mapfn = fofcmap
            if usefreer == 1 : 
                fftnew(mtzfn, fofcmap, pdbfile,f1label,f2label,phiclabel,sigf1label,n,(n)-1,freeRlabel)
                fftnew(mtzfn, fcmap, pdbfile ,f1label,f2label,phiclabel,sigf1label,0,1,freeRlabel,1)
            else :
                fftnew(mtzfn, fofcmap, pdbfile , f1label,f2label,phiclabel,sigf1label,n,(n)-1,None)
                fftnew(mtzfn, fcmap, pdbfile ,f1label,f2label,phiclabel,sigf1label,0,1,freeRlabel,1)



        if useomitmap == 1 :
            omit_mapfile = "%s%dFo-%dFC.omit.map"%(mtzfn,n,((n)-1))
            mapfn = omit_mapfile 
            if usefreer == 1 : 
                omitsucc = omitmap(mtzfn, omit_mapfile, f1label,f2label,phiclabel,sigf1label,n,(n)-1,freeRlabel)
            else :
                omitsucc = omitmap(mtzfn, omit_mapfile, f1label,f2label,phiclabel,sigf1label,n,(n)-1)
            if omitsucc == None :
                useomitmap = 0
                useccp4map = 1
                mapfn = fofcmap
        if usecnsmap == 1 :
            fofcmap = mtzfn+"2fofc.map"
            fcmap = mtzfn+"fc.map"
            mapfn = fofcmap
            
        
        if poorOnly == 1 :
            
           
            xscorer = XrayScorer(None, xscoreCutoff)
            mapcoeff = "%dF1-%dF2"%(n, (n)-1)
            badmcsc,id1 = xscorer.score(pdbfile, fofcmap, fcmap, f2label, phiclabel, mapcoeff, "mcsc",verbose)
            badmc,id2 = xscorer.score(pdbfile, fofcmap, fcmap, f2label, phiclabel, mapcoeff, "mc",verbose )
            badpept,id3 = xscorer.score(pdbfile, fofcmap, fcmap, f2label, phiclabel, mapcoeff, "pept",verbose )
            badsides,id4 = xscorer.score(pdbfile, fofcmap, fcmap, f2label, phiclabel, mapcoeff, "sc",verbose)
            

            bad = list ( set(badmcsc+badmc+badpept) )
            badid = list ( set(id1+id2+id3) )

            prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
            res, resids, resnums, resns, chids, inscodes, pts , atomids , bfac = prot2res.readProtRes2(prot)

            ## IF START AND STOP ARE GIVEN, THEN GET SUBSET OF BAD RESIDUES
            if start != None and stop != None and chainid != None :
                badresids = []
                for bd in badid :
                    if int(resnums[bd]) in range(start,stop+1) and chids[bd] == chainid :
                        badresids.append(resids[bd])

                badscs = []
                for bd in id4 :
                    if int(resnums[bd]) in range(start,stop+1) and chids[bd] == chainid :
                        badscs.append(resids[bd])
            
            else:
                badresids = bad
                badscs = badsides

                
            if len(badresids) == 0 and len(badscs) == 0  :
                print "No poor region identified, try higher cutoff if required" ;
                import sys ; sys.exit(0)
            if len(badresids) == 0:
                print "No poor fitting mainchain atoms identified, only sidechains will be rebuilt"
                print "Sidechains to be rebuilt are ", badscs
        

    ############### SET UP   ######################
    if (mapfn !=  None or mtzfn !=None ):
        esmin, esmax, esmean, rcmult, xscoreCutoff = minXSig, maxXSig, .0, 5, poorThreshold
        mapcoeff = "%dF1-%dF2"%(2, 1)
        xrayRestGen.append( prepareChainV5.XrayRestraintsGenerator(mapfn, "map", f2label, phiclabel, mapcoeff, esmin, esmax, esmean, [], edOpt ) )

    

    ############### GENERATE MODELS  ######################
    multiPrepC = prepareChainV5.PrepareChain(rotLib)
    for  i in range(nmodels):
        nb = 0 ; cp1 = 0
        print
        print nmodels
        if nmodels > 1 :
            outpdb = "pre."+str(i)+"."+rtkmodel
            outpdbsc = str(i)+"."+rtkmodel
        else :
            outpdb = "pre."+rtkmodel
            outpdbsc = rtkmodel

        if (caRad != None and caRad < 0.49) :
            print "Minimum Ca-radii needs to be 0.5 ANGSTROM."
            
            if nmodels == 1 : 
                print "Copying", pdbfile ,  " to ","pre."+rtkmodel 
                shutil.copyfile(pdbfile, "pre."+rtkmodel)
            else :
                print "Copying", pdbfile ,  " to ","pre."+str(i)+"."+rtkmodel 
                shutil.copyfile(pdbfile, "pre."+str(i)+"."+rtkmodel)
            nb = 1
            
            cp1 = 1 
        
        elif ( ( (badresids != None and  len(badresids) > 0  and poorOnly ==  1 ) or (poorOnly == 0 and badresids == None ) or (poorOnly == 0 and len(badresids) > 0) ) and sconly == 0 ):
            print "\nGenerating MODEL %d of %d asked for ..........."%(i+1,nmodels)


            ml = Multiloop(modelIn, badresids, mconly, caRad, scRad, scvdwr, guidedsampling, popsize, backtrack, 1 , outpdb, xrayRestGen, multiPrepC,sconly, cacaCutoff,nullres,modelN2C,natt,missingpts,1,allOpt,edOpt,None,addsc) 
            
            ml.ranker = None 
            if poorOnly == 1 :
                ml.poor = "poor"
            if start != None and stop != None and chainid != None :
                ml.loopstartnum = start
                ml.loopstopnum = stop
                ml.loopclosure = closure


            if mtzfn != None or mapfn !=None:
                ml.ranker = XrayRanker(mapfn, "map", f2label, phiclabel, mapcoeff, esmin, esmax)
                ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
                ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble = None
                ml.cellsym = [ a, b, c, alpha, beta, gamma, sgtable[sg][0] ]

            if snappy == 'True' :
                ml.snappy = 1
                
                
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
            #removeMODEL(outpdb) ## automatic water addition can go here
            #copyheader(outpdb,"modelinit.pdb")
            #copyNonprotein(pdbfile,outpdb)
            #adjustBfac2(outpdb, pdbfile,nativeBfac, mcBfac, scBfac)
            shutil.copyfile(outpdb, outpdbsc  )
            print "Renaming", outpdb, " to ", outpdbsc

        elif  badresids != None and len(badresids)==0 and poorOnly == 0 and sconly != 1 :
            #removeMODEL(outpdb) ## automatic water addition can go here
            #copyheader(outpdb,"modelinit.pdb")

            #copyNonprotein(pdbfile,outpdb)
            #adjustBfac2(outpdb, pdbfile,nativeBfac, mcBfac, scBfac)
            shutil.copyfile(outpdb, outpdbsc  )
            print  "Renaming", outpdb, " to ", outpdbsc
            print "WARN!! Nothing to be rebuilt, renaming", outpdb, " to ", outpdbsc
            
            
        elif  nb != 1 and sconly !=1 : 
            #removeMODEL(outpdb) ## automatic water addition can go here
            #copyheader(outpdb,"modelinit.pdb")
            #copyNonprotein(pdbfile,outpdb)
            #adjustBfac2(outpdb, pdbfile,nativeBfac, mcBfac, scBfac)
            shutil.copyfile(outpdb, outpdbsc  )
            print  "Renaming", outpdb, " to ", outpdbsc
            print "WARN!! Nothing to be rebuilt, renaming", outpdb, " to ", outpdbsc



        else :
        
            print "\n\nStarting OPSAX ........"
            if badresids == None :
                prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
                res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
                badresids = list( resids.values() )
                

            if cp1 == 1 and  poorOnly == 1 :
                badresids = list(  set( badscs )  )


            elif cp1 == 1  and poorOnly == 0 :
                badresids = badresids


            elif ((len(badresids)==0 and poorOnly == 1 )  or (sconly == 1 and poorOnly == 1)):

                badresids = list(  set( badscs )  )

            elif poorOnly == 1 :
                
                if len(unbuilt[0]) == 0  :
                    badresids = list(  set(  findChangedSC(pdbfile, outpdb) + badscs )  )

                else :
                    badresids = list(  set( badscs )  )

                    
            elif   poorOnly == 0 and sconly != 1 and len(unbuilt[0]) == 0   :
                badresids = list(  set(  findChangedSC(pdbfile, outpdb) )  )

                
            elif   poorOnly == 0 and sconly != 1 and len(unbuilt[0]) != 0  :
                badresids = list(  set(  findChangedSC(pdbfile, outpdb) )  )

            else:
                badresids = badresids

                ### todo remove missing sidechains ??


            if len(badresids) == 0:
                print "No residues marked for OPSAX"
                
                shutil.copyfile(outpdb , outpdbsc)
            elif len(badresids) > 250:
                
                print "Length of section > 250 residues, can not perform OPSAX"
                shutil.copyfile(outpdb , outpdbsc)

            else :

                prot = protein(outpdb, read_hydrogens=0, read_waters=0, read_hets=0)
                res, resids, resnums, resns, chids, inscodes, pts , atomids , bfac = prot2res.readProtRes2(prot)    
                for k , v in resids.items():
                    if v in badresids:
                        missing = [' N  ', ' CA ' , ' C  ', ' O  '] 
                        for a , b in res[k].items() :
                            if a in missing :
                                missing.remove(a)
                        if len(missing) != 0 :
                            print  "In order to build sidechains, mainchain atoms ", missing, "of residue ", v ,"needs to be present"
                            import sys ; sys.exit()
                            
                            
                scPrepC, useGivenRot, useDEE = prepareChainV5.PrepareChain("PRL"), userot , 1
                #if mtzfn != None :
                 #   scp = SCplacement("pre."+str(i)+"."+rtkmodel , scReduction, str(i)+"."+rtkmodel, "dotfile", useDEE, mtzfn, f1label, f2label , phiclabel, "2F1-F2", esmin, esmax, addsc, useGivenRot, badresids, scPrepC).run()
                if mtzfn != None and useomitmap == 1 :
                    nb2 = SCplacement(outpdb, scReduction, outpdbsc, "dotfile", useDEE, omit_mapfile, f1label, f2label , phiclabel, mapcoeff, esmin, esmax, addsc, useGivenRot, badresids, scPrepC).run()

                elif mtzfn != None : 
                    nb2 = SCplacement(outpdb, scReduction, outpdbsc, "dotfile", useDEE, fofcmap, f1label, f2label , phiclabel, mapcoeff, esmin, esmax, addsc, useGivenRot, badresids, scPrepC).run()
                    
                elif mapfn != None : 
                    nb2 = SCplacement(outpdb, scReduction, outpdbsc , "dotfile", useDEE, mapfn, f1label, f2label , phiclabel, mapcoeff, esmin, esmax, addsc, useGivenRot, badresids, scPrepC).run()
                                        

                #if nb2 != 1:
                 #   shutil.copy("pre."+str(i)+"."+rtkmodel ,str(i)+"."+rtkmodel)
                #scPrepC, useGivenRot, useDEE = prepareChainV5.PrepareChain("PRL"), userot , 1

                #if mapfn != None : 
                #    nb2 = SCplacement(outpdb, scReduction, outpdbsc, "dotfile", useDEE, mapfn, f1label, f2label , phiclabel, mapcoeff, esmin, esmax, addsc, useGivenRot, badresids, scPrepC).run()
                #if nb2 != 1:
                #    shutil.copyfile(outpdb , outpdbsc)
            print "\n Finished OPSAX ........"

        adjustBfac2(outpdb, pdbfile,nativeBfac, mcBfac, scBfac)
        adjustBfac2(outpdbsc, pdbfile,nativeBfac, mcBfac, scBfac)

        removeMODEL(outpdbsc) ## automatic water addition can go here
        removeMODEL(outpdb) ## automatic water addition can go here
        
        copyNonprotein(pdbfile,outpdbsc) ## automatic water addition can go here
        copyNonprotein(pdbfile,outpdb)

        copyheader(outpdb,"modelinit.pdb")
        copyheader(outpdbsc,"modelinit.pdb")
        
        print "------------------------------------ MODEL ASSESSMENT --------------------------------"
        
        report = "False"
        if report == "False" :
            print "Normal termination"

        elif mtzfn != None :
            from prefRapperV4WithRefmac import refmacNew
            from protRefine_nick import molProbity_badres, add_cryst_card 


            sfcheck(mtzfn,pdbfile,f1label,sigf1label,freeRlabel ) 
            if (os.path.isfile("sfcheck.log")==False) :
                print "no sfcheck file ";
                
            else :
                head1,tail1  = os.path.split(pdbfile) ; sfcheckfile0 = "sfcheck.%s.log"%tail1
                shutil.copyfile("sfcheck.log",sfcheckfile0)
                

            if sfall(outpdbsc, mtzfn, "phased."+str(i)+".mtz",resolution,usefreer,f1label,sigf1label,freeRlabel ) == None :
                print "sfall failed"
                import sys ; sys.exit()
            sfcheck(dir_xyzout+"/phased."+str(i)+".mtz", outpdbsc,f1label,sigf1label,freeRlabel ) 

            if (os.path.isfile("sfcheck.log")==False) :
                print "no sfcheck file ";
                
            else :
                prefix = outpdbsc ; sfcheckfile1 = "sfcheck.%s.log"%prefix
                os.rename("sfcheck.log",sfcheckfile1)

            
            sfall(outpdb, mtzfn, "phased.pre."+str(i)+".mtz",resolution,usefreer,f1label,sigf1label,freeRlabel ) 
            sfcheck(dir_xyzout+"/phased.pre."+str(i)+".mtz",outpdb,f1label,sigf1label,freeRlabel ) 


            if (os.path.isfile("sfcheck.log")==False) :
                print "no sfcheck file ";
            else :
                prefix = outpdb; sfcheckfile2 = "sfcheck.%s.log"%prefix
                os.rename("sfcheck.log",sfcheckfile2)

            


            print
            print
            print "R-factor and correlation for initial coordinates"
            parse_sfcheck_logfile(sfcheckfile0,pdbfile)
            print

            print
            print
            print "R-factor and correlation for output (rappertk) coordinates before opsax"
            parse_sfcheck_logfile(sfcheckfile2,outpdb)

            
                       
            print
            print "R-factor and correlation for output (rappertk) coordinates"
            parse_sfcheck_logfile(sfcheckfile1,outpdbsc)






        elif mapfn != None :
            from prefRapperV4WithRefmac import refmacNew
            from protRefine_nick import molProbity_badres, add_cryst_card 
            from copyheader import copyheader


            sfcheck(mapfn,pdbfile,None,None,None)
            os.rename("sfcheck.log","sfcheck.init.log")
            
            sfcheck(mapfn, outpdbsc,None,None,None)
            os.rename("sfcheck.log","sfcheck."+str(i)+".log")

            print
            print
            print "R-factor and correlation for initial coordinates"
            parse_sfcheck_logfile("sfcheck.init.log",pdbfile)
            print
            print
            print "R-factor and correlation for output (rappertk) coordinates"
            parse_sfcheck_logfile("sfcheck."+str(i)+".log",outpdbsc)


        print "------------------------------------ MODEL ASSESSMENT ENDS ---------------------------"






def parse_sfcheck_logfile(filename,pdb):
    import re
    lines = open(filename, 'r').readlines()
#    pattern = re.compile('(\d+.\d+)')
    for l in lines :
#        print l
        if 'R-factor           ' in l :
            print pdb,l
            #params = pattern.match(l)
            #if params != None:
            #    rpar = list(params.groups())
            #    rfactor = rpar[1])

        if 'Correlation factor ' in l :
            print pdb,l
        
    
    
        if " Rfree,Rrest" in l :
            print pdb,l






def parseLoopSequence(pdbfile,loopres,start,stop,chainid,mconly):

    null = [] ; missingpts = {}
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts , atomids , bfac = prot2res.readProtRes2(prot)
    seqlen = len(loopres) ;
    sslen = 0 

    badresids = [] 
    missing = 0

    res.keys().sort()
    
    for index in range(start,stop+1):
        found = 0 

        for k in res.keys():

            if int(resnums[k]) == index and chids[k] == chainid : 
                found = 1 
                if not isAAres(resns[k]) :
                    continue
                if ' CA ' in res[k].keys():
                    badresids.append(resids[k])
                else :
                    missing = 1

                sslen= sslen + 1
                from data import three2one
                print three2one[resns[k]],

        #if found == 0 :
        #missing = 1
                
    if (seqlen == sslen):
        if missing  == 0 :
            print "Returning ... ",seqlen, sslen
            return badresids,0,null,missingpts

        else :
            #from prepareChainV4 import  incompleteREScorrection2
            #nres, nresids, nresnums, nresns, nchids, ninscodes, npts = {}, {}, {}, {}, {}, {}, []
            #nres, nresids, nresnums, nresns, nchids, ninscodes, npts , null,missingpts  = incompleteREScorrection2(res, resids, resnums, resns, chids, inscodes, pts,start,stop,chainid,loopres,mconly)
            #modelRenderer = ModelRenderer(nres, nresns, nchids, nresnums, ninscodes, [],pdbfile + ".fixed")
            #modelRenderer.render(npts)

            return badresids,0,null,missingpts

        
            #import sys ; sys.exit()

    else :
        print "Warning!! Sequence length specified by loopres does not equal  length of section %s %s %s read from PDB file "%(start,stop,chainid)
        print "Ignoring sequence from coordinate file"
        
    sslen = stop - start + 1
    print sslen
    sslen = getdiff(start,stop)
    print sslen,len(loopres)

    if sslen != seqlen :
        print "(1) Sequence length does not match start and stop points",seqlen, sslen
        import sys ; sys.exit()

    from prepareChainV5 import  incompleteREScorrection2
    nres, nresids, nresnums, nresns, nchids, ninscodes, npts = {}, {}, {}, {}, {}, {}, []
    nres, nresids, nresnums, nresns, nchids, ninscodes, npts , null,missingpts  = incompleteREScorrection2(res, resids, resnums, resns, chids, inscodes, pts,start,stop,chainid,loopres,mconly)
    modelRenderer = ModelRenderer(nres, nresns, nchids, nresnums, ninscodes, [],pdbfile + ".fixed")
    modelRenderer.render(npts)

    
    rev_resnum = {} ; badresids = [] ; chainkey = [];    missing = 0
    sslen = 0 
        
    for index in range(start,stop+1):
        found = 0 
        for k in nres.keys():
            if int(nresnums[k]) == index and nchids[k] == chainid : 
                found = 1 
                badresids.append(nresids[k])
                sslen = sslen + 1
        if found == 0 :
            missing = 1

    if (seqlen != sslen):
        print "(1) Sequence length does not match start and stop points",seqlen, sslen , badresids
        import sys ; sys.exit()

    if missing  == 0 :
        return badresids,0,null , missingpts
    
    else:
        print "Cant build.Error 201, please report problem to ak459@cam.ac.uk"
        import sys ; sys.exit(1)



def getdiff(start,stop):
    sslen = None
    if stop < start :
        print "Error,stop residuemust be greater than start residue number";
        import sys ; sys.exit()
    if start < 0 and stop <= 0:
        sslen = abs(start)-abs(stop) + 1 
    if start >=0  and stop >= 0 :     
        sslen = stop -abs(start) + 1 
    if start <= 0 and stop > 0 :                 
        sslen = abs(start) + abs(stop) + 1
    if sslen == None:
        print "Sequence error..exiting..report error to ak459@cam.ac.uk"
        import sys ; sys.exit()
    return sslen
            
def parseLoopSequenceSC(pdbfile,loopres,start,stop,chainid,mconly):
    from data import resAtoms
    from pdbr import isAAres
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    badresids = []
    for k , v in resnums.items():
        if not isAAres(resns[k]) :
            continue

        if int(v) >= start and int(v) <= stop and chids[k] == chainid :
            if ' CA ' in res[k].keys():
                badresids.append(resids[k])
            else :
                print "In order to use c-alpha restraints for building, all c-alpha atoms in the section need to be present in the input coordinates file"
                print "The C-alpha atom for residue [%s] is missing in the input coordinates "%resids[k]
                #import sys ; sys.exit()
    return badresids

    










def replaceWaters(pdbfile, ccp4map, peakCutoff=2.5) :
    from pdbr import isPdbAAline, isPdbAtomLine, line2resn
    lines = []
    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        import sys ; sys.exit()
    
    for l in open(pdbfile, 'r').readlines() :
        if isPdbAtomLine(l) and line2resn(l) == "HOH" : continue
        if l ==  "END\n" or l[0:4] == "END " : continue
        lines.append(l)
    opf = open(pdbfile, 'w')
    for l in lines : opf.write(l)
    opf.close()
    cmd = "findwaters --pdbin %s --map %s --pdbout %s --sigma %f" % (pdbfile, ccp4map, "findwat", peakCutoff)
    execCmd(cmd, [])
    if (os.path.isfile("findwat")==False) :
        print "Cannot find file findwat"
        print "No file in directory ", os.getcwd()
        import sys ; sys.exit()
    for l in open('findwat', 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        lines.append(l)
    lines.append("END\n")
    opf = open(pdbfile, 'w')
    for l in lines : opf.write(l)
    opf.close()

def randomize(doRandomize) :
    import misc, random
    if doRandomize :
        random.seed(doRandomize)
        misc.RanGen.instance().seedme( doRandomize )
    else : 
        random.seed(1973) # for Paul!
        misc.RanGen.instance().seedme(1942) # for Tom!

def adjustBfacOccu(pdbfile, refpdb) :
    print "\nCopy Bfactors from", refpdb, "to", pdbfile, "if coordinates are same, else set to 30"
    from pdbr import line2bfac, isPdbAtomLine, line2occu, line2bfac, line2atomid, changeBfactorOccupancy, line2crdstr
    id2bo = {} # read in crd, bfac, occu for each atom-id


    if (os.path.isfile(refpdb)==False) :
        print "Cannot find file %s "%refpdb
        print "No file in directory ", os.getcwd()
        import sys ; sys.exit()




    for l in open(refpdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if id2bo.has_key(line2atomid(l)) : continue
        id2bo[line2atomid(l)] = [line2bfac(l), line2occu(l), line2crdstr(l)]
    newlines = []

    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        import sys ; sys.exit()

    for l in open(pdbfile, 'r').readlines() :
        if not isPdbAtomLine(l) : newlines.append(l) ; continue
        bfac, occu, crdstr = id2bo[ line2atomid(l) ]
        if crdstr == line2crdstr(l) : l = changeBfactorOccupancy(l, bfac, occu)
        else : l = changeBfactorOccupancy(l, 30.0, occu)
        newlines.append(l)
    op = open(pdbfile, 'w')
    for l in newlines : op.write(l)
    op.close()


def callmain() :
    from commonOptions import makeParser, parseOptions
    import optparse
    parser = optparse.OptionParser()

    ################# I/O Files and directories ##########################################
    
    parser.add_option("--xyzin", action='store', type='string', dest='pdbfile', help='Name of input coordinates file: Complete PATH needed')
    parser.add_option("--xyzout", action='store', type='string', dest='pdbout', help='Name of output coordinates file',default="modelout.pdb")
    parser.add_option("--dir-xyzout", action='store', type='string', dest='dir_xyzout', help='Name of directory in which the output files will be created. Complete PATH needed')
    parser.add_option("--mapin", action='store', type='string', dest='mapfn', help='Name of input MAP file (complete PATH needed)',default=None)
    parser.add_option("--hklin", action='store', type='string', dest='mtzfn', help='Name of input (phased) MTZ file (complete PATH needed)',default=None)


    ### Assignment of restraints and radii #############################################
    
    parser.add_option("--use-ca-restraints", action='store', dest='caRes', help='[True/False], Apply positional restraints on the C-alpha atoms',default="True")
    parser.add_option("--use-sc-restraints", action='store', dest='scRes',type= 'string', help='[True/False],  Apply positional restraints on the centroid of the sidechain atoms',default="True",)
    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='Radius of spherical restraint ( Angstrom ) on the C-alpha atom position', default=1.0)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='Radius of spherical restraint ( Angstrom ) on the centroid of the sidechain atoms', default=2.0)

    ############ General parameters #####################################################
    
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='Factor by which to reduce effective Van der Waals distance for sidechain atoms', default=0.75)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='Population size', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='Level of printed output [1-10], 1:Concise output log, 10: Detailed output log', default=6)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='[numsteps]X[stepsize] e.g. when set to 4X5, In case of failure to build at a residue the program will  backtrack 4 times by 5 residues. For detailed help see Rappertk webpage', default="4X5")
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
    parser.add_option("--FREER", action='store', type='string', dest='freeRlabel', help='Column label for FreeR in MTZ file', default=None)



    ############# Residues to be modelled ####################################
    

    parser.add_option("--rebuild-poor-regions-only", action='store', type='string', dest='poorOnly', help='[True/False] Rebuild regions ofinput structure with poor fit to an electron density map. Residues to be rebuilt are identified using a real space correlation coefficientscore, the cut-off for which is set using --poor-fit-threshold.', default="False")
    parser.add_option("--poor-fit-threshold", action='store', type='float', dest='poorThreshold', help='Correlation coefficient threshold to identify poor fitting regions', default=0.9)



    parser.add_option("--loopseq", action='store', type='string', dest='loopres', help='Amino acid sequence for loop to be built', default=None)
    parser.add_option("--use-loopclosure-restraints", action='store', type='string', dest='closure', help='Use geometric restraints to ensure closure of loop with anchor residues', default= "True")
    parser.add_option("--start", action='store', type='int', dest='start', help='Residue number to start building from ', default=None)
    parser.add_option("--stop", action='store', type='int', dest='stop', help='Residue number to stop building at', default=None)

    parser.add_option("--start-inscode", action='store', type='string', dest='startcode', help='Insertion code of start residue ', default=' ')
    parser.add_option("--stop-inscode", action='store', type='string', dest='stopcode', help='Insertion code of stop residue', default=' ')

    parser.add_option("--chainid", action='store', type='string', dest='chainid', help='Chain ID of section to be built.', default=None)
#    parser.add_option("--modelN2C", action='store', type='string', dest='modelN2C', help='[True/False] Model fragment without loop closure restraints. Used in conjunction with --start, --stop, --chainid. Requires --use-ca-restraints True ', default="False")

    ######### Ouptut parameters #############################################

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


    parser.add_option("--get-snapshot", action='store', type='string', dest='snappy', help='Output the population of conformers at each residue extension step [True/False]', default='False')

    parser.add_option("--FCmap", action='store', type='string', dest='FCmap', help='FC map', default=None)
    #parser.add_option("--maptype", action='store', type='string', dest='mapformat', help='Options are : 2fofc', default="2fofc")

    options = parseOptions(parser)

    from checkfornone import checkfornone 
    

    options.maptype = "2fofc"
    options.pdbfile,options.pdbout,options.dir_xyzout,options.mapfn,options.mtzfn,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,options.nmodels,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freeRlabel,options.poorOnly,options.poorThreshold,options.loopres,options.start,options.startcode,options.stop,options.stopcode,options.chainid,options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt,options.closure,options.addsc,options.userot = checkfornone(options.pdbfile,options.pdbout,options.dir_xyzout,options.mapfn,options.mtzfn,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,options.nmodels,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freeRlabel,options.poorOnly,options.poorThreshold,options.loopres,options.start,options.startcode,options.stop,options.stopcode,options.chainid,options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt,options.closure,options.addsc,options.userot)

    main(options.pdbfile,options.pdbout,options.dir_xyzout,options.mapfn,options.mtzfn,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,options.nmodels,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freeRlabel,2,options.poorOnly,options.poorThreshold,options.loopres,options.start,options.startcode,options.stop,options.stopcode,options.chainid,options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt,options.closure,options.addsc,options.userot,options.maptype,options.snappy)

if __name__ == "__main__" :
    callmain()
