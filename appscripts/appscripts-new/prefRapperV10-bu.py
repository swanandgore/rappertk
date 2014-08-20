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
from loopbuildV7 import Multiloop, incompleteSCcorrection2
import prepareChain
import prepareChainV2
import prepareChainV3
from multProtref import joinPDBs, splitPDBs, restoreChainids
import checkProtChains
import checkProtChainsV2
from peptidebuild import ModelRenderer
import misc
from pdbr import isPdbAAline, isPdbAtomLine, line2resn , line2resnum , line2chid , line2atomnumber
from pdbr import makeResid

## Last modified 18-02-2010 :  Using fft for mtz to ccp4map





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
    #import sys; sys.exit()
    for k,v in chids.items():
        if v not in allchains:
            allchains.append(v)
    
    return allchains





def main(pdbfile,pdbout,dir_xyzout,mapfn,mtzfn,caRes,scRes,caRad,scRad,scReduction,popsize,verbose,backtrack,rotLib,nmodels,mconly,sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freeRlabel,n,poorOnly,poorThreshold,loopres,start,stop,chainid,modelN2C,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt):



    str2bool = {"False": 0 , "True":1}

    print "========================== Now  printing all KEYWORD values ==========================="




    print
    print "                           Directories and Files"
    print
    
    print "	--xyzin ".ljust(50), 
    print "%30s "%pdbfile 
    print "	--xyzout ".ljust(50), 
    print "%30s "%pdbout 
    print "	--dir-xyzout ".ljust(50), 
    print "%30s "%dir_xyzout 
    print "	--mapin ".ljust(50), 
    print "%30s "%mapfn 
    print "	--hklin ".ljust(50), 
    print "%30s "%mtzfn 
    print "	--use-ca-restraints ".ljust(50), 

    print
    print "                           Restraint paramaters "
    print

    print "%20s "%caRes 
    print "	--use-sc-restraints ".ljust(50), 
    print "%20s "%scRes 
    print "	--ca-restraint-radius ".ljust(50), 
    print "%20s "%caRad 
    print "	--sc-centroid-restraint-radius ".ljust(50), 
    print "%20s "%scRad





    print
    print "                           Build paramaters "
    print
    print "	--sidechain-vdw-reduction ".ljust(50), 
    print "%20s "%scReduction 
    print "	--population-size ".ljust(50), 
    print "%20s "%popsize 
    print "	--verbose ".ljust(50), 
    print "%20s "%verbose 
    print "	--backtrack ".ljust(50), 
    print "%20s "%backtrack 
    print "	--mconly ".ljust(50), 
    print "%20s "%mconly 
    print "	--sconly ".ljust(50), 
    print "%20s "%sconly 
    print "	--make-all-restraints-optional ".ljust(50), 
    print "%20s "%allOpt 


    print "	--loopseq ".ljust(50), 
    print "%20s "%loopres 
    print "	--start ".ljust(50), 
    print "%20s "%start 
    print "	--stop ".ljust(50), 
    print "%20s "%stop 
    print "	--chainid ".ljust(50), 
    print "%20s "%chainid 


    print
    print "                           Crystallographic paramaters "
    print
    print "	--opsax ".ljust(50), 
    print "%20s "%opsax 
    print "	--use-freer ".ljust(50), 
    print "%20s "%usefreer 
    print "	--rotamerlib ".ljust(50), 
    print "%20s "%rotLib 
    print "	--num-models ".ljust(50), 
    print "%20s "%nmodels 
    print "	--a ".ljust(50), 
    print "%20s "%a 
    print "	--b ".ljust(50), 
    print "%20s "%b 
    print "	--c ".ljust(50), 
    print "%20s "%c 
    print "	--alpha ".ljust(50), 
    print "%20s "%alpha 
    print "	--beta ".ljust(50), 
    print "%20s "%beta 
    print "	--gamma ".ljust(50), 
    print "%20s "%gamma 
    print "	--sg ".ljust(50), 
    print "%20s "%sg 
    print "	--resolution ".ljust(50), 
    print "%20s "%resolution 
    print "	--FP ".ljust(50), 
    print "%20s "%f1label 
    print "	--SIGFP ".ljust(50), 
    print "%20s "%sigf1label 
    print "	--FC ".ljust(50), 
    print "%20s "%f2label 
    print "	--PHIC ".ljust(50), 
    print "%20s "%phiclabel 
    print "	--FREER ".ljust(50), 
    print "%20s "%freeRlabel 


    print "	--n ".ljust(50), 
    print "%20s "%n 
    print "	--rebuild-poor-regions-only ".ljust(50), 
    print "%20s "%poorOnly 
    print "	--poor-fit-threshold ".ljust(50), 
    print "%20s "%poorThreshold 
    print "	--default-mainchain-b-factor ".ljust(50), 
    print "%20s "%mcBfac 
    print "	--default-sidechain-b-factor ".ljust(50), 
    print "%20s "%scBfac 
    print "	--models-get-native-bfactors ".ljust(50), 
    print "%20s "%nativeBfac 
    print "	--minimum-sig ".ljust(50), 
    print "%20s "%minXSig 
    print "	--maximum-sig ".ljust(50), 
    print "%20s "%maxXSig 
    print "	--cacaCutoff ".ljust(50), 
    print "%20s "%cacaCutoff 
    print "	--modelN2C ".ljust(50), 
    print "%20s "%modelN2C 
    print "	--make-ed-optional ".ljust(50), 
    print "%20s "%edOpt 

    
 
    print
    print
    print "========================== End KEYWORD values =========================================="
    print
    print
    

    
    randomize(1)

    
                                                                       
    ######### Check for / make out directory #############
    if dir_xyzout == None :
        print "Directory for output files needs to be set. --dir_xyzout needs to be set",dir_xyzout
        import sys ; sys.exit()

    if  os.path.isfile(dir_xyzout) :
        print dir_xyzout,"Is a file"
        print "Please rename output directory"
        import sys ; sys.exit()
    head,tail  = os.path.split(dir_xyzout)
    if head and os.path.isdir(head):
        print "Checking for directory ............",head
    else:
        print "%s does not exists"%head
        import sys ; sys.exit()
    if not os.path.isdir(dir_xyzout) :
        os.mkdir(dir_xyzout)
    os.chdir(dir_xyzout)

    
    head,tail  = os.path.split(pdbout)
    if head != "" : 
        pdbout = tail

    ######### End check for output directory #############





        
    ############### SET UP i/O Coordinate file  ######################
    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s " %pdbfile ; import sys ;        sys.exit()
    shutil.copyfile(pdbfile, "%s/%s" % (dir_xyzout,"modelinit.pdb"))

    modelIn =  "modelinit.pdb" ;    rtkmodel = pdbout 
    pdboutpath = "%s/%s"%(dir_xyzout,pdbout)




    #### Initializations  ####
    
    opsaxoff = None ; guidedsampling = None;    scvdwr = scReduction ;    popsize = popsize
    xrayRestGen = [] ;    badresids = [] ;         nullres = [] ; missingpts = []



    ###########     CHECKS  FOR INTEGRITY OF INPUT VALUES                    #############################

    if mconly not in ["True", "False"]:
        print "--mconly needs to be either True or False. Incorrect value input " , mconly
        import sys ; sys.exit()

    if sconly not in ["True", "False"]:
        print "--sconly needs to be either True or False. Incorrect value input " , sconly
        import sys ; sys.exit()

    if modelN2C not in ["True", "False"]:
        print "--modelN2C needs to be either True or False, Incorrect value input " ,modelN2C
        import sys ; sys.exit()

    if scRes not in ["True", "False"]:
        print "--use-sc-restraints needs to be either True or False, Incorrect value input " ,scRes
        import sys ; sys.exit()

    if caRes not in ["True", "False"]:
        print "--use-ca-restraints needs to be either True or False , Incorrect value input",caRes
        import sys ; sys.exit()

    if nativeBfac not in ["True","False"]:
        print "--models-get-native-bfactors needs to be either True or False. Incorrect value input",nativeBfac
        import sys ; sys.exit()

    if opsax not in ["True","False"]:
        print "--opsax needs to be either True or False.  Incorrect value input",opsax
        import sys ; sys.exit()

    if usefreer not in ["True","False"]:
        print "--use-freer needs to be either True or False. Incorrect value input",usefreer
        import sys ; sys.exit()
        
    if edOpt  not in ["True","False"]:
        print "--make-ed-optional needs to be either True or False.  Incorrect value input",edOpt
        import sys ; sys.exit()

    if allOpt not in ["True","False"] :
        print "--make-all-restraints-optional needs to be either True or False .  Incorrect value input",allOpt
        import sys ; sys.exit()
        
    if poorOnly not in ["True","False"]:
        print "--rebuild-poor-regions-only needs to be either True or False,  Incorrect value input",poorOnly
        import sys ; sys.exit()

####


    
    if modelN2C == "True" and caRes == "False" :
        print "When building as a fragment (i.e. --modelN2C True), --use-ca-restraints needs to be set to True"
        import sys ; sys.exit()
        

    if (scReduction < 0. or scReduction > 1.) and mconly != "True" :
        print "value of --sidechain-vdw-reduction needs to be >= 0.0 and <= 1." ;
        print "Incorrect value given", scReduction
        import sys ; sys.exit()
        
    if rotLib not in ["PRL", "SCL1.0","SCL0.5","SCL0.2"]:
        print "Unrecognised rotamer library, options are PRL, SCL1.0, SCL0.5 ,  SCL0.2"
        print "Incorrect value given", rotLib
        import sys ; sys.exit()


    if caRad < 0. and caRes == "True" :
        print "value of --ca-restraint-radius  needs to be positive value " ;
        print "Incorrect value given", caRad
        import sys ; sys.exit()

    if scRad < 0. and scRes == "True" :
        print "value of --sc-centroid-restraint-radius  needs to be positive value " ;
        print "Incorrect value given", scRad
        import sys ; sys.exit()


    if popsize < 100  :
        print "WARNING !! value of --population-size  should be atleast 100"
        if popsize < 0  :
            print "value of --population-size  should be atleast 100"
            print "Incorrect value given", popsize
            import sys ; sys.exit()                

    if verbose < 1  :
        print "WARNING !! value of --verbose  should be atleast 1"
        print "Setting to default value of 7"
        verbose = 7
    misc.setVerbosity(verbose) 



    if natt < 1  :
        print "WARNING !! value of --attempts  should be atleast 1"
        print "Setting to default value of 5"
        natt = 5


    if cacaCutoff < 0. :
        print "WARNING !! value of --cacaCutoff  should be positive"
        print "Setting to default value of 5.0"
        cacaCutoff = 5.0

    if backtrack == "None":
        backtrack = None
    else :
        if backtrack != None : 
            xpos = backtrack.find('X')
            if xpos <= 0:
                print "Error in specifying backtrack, it should be of the form of numstepsXstepsize, where numsteps and stepsize are integer values. For e.g. --backtrack 4X5 , in case of failure to build at a residue the program will  backtrack 4 times by 5 residues."
                print "Incorrect value given", backtrack
                import sys ; sys.exit()

            try : 
                int(backtrack[0:xpos]) ; int(backtrack[xpos+1:])

            except ValueError : 
                print "Error in specifying backtrack, it should be of the form of numstepsXstepsize, where numsteps and stepsize are integer values. For e.g. --backtrack 4X5 , in case of failure to build at a residue the program will  backtrack 4 times by 5 residues."
                print "Incorrect value given ", backtrack
                import sys ; sys.exit()



    if opsax == "True" and (mapfn == None and mtzfn == None) :
        print "WARNING!! OPSAX will  only be used when MTZ/MAP file is given"
        opsax = "False"


    mconly =  str2bool[mconly]
    sconly =  str2bool[sconly]
    opsax =  str2bool[opsax]
    usefreer =  str2bool[usefreer]
    poorOnly =  str2bool[poorOnly]
    nativeBfac =  str2bool[nativeBfac]
    modelN2C =  str2bool[modelN2C]
    caRes =  str2bool[caRes]
    scRes =  str2bool[scRes]
    edOpt =  str2bool[edOpt]
    allOpt =  str2bool[allOpt]

    if mconly == 0 :
        mconly = None

    if allOpt == 1 :
        edOpt = 1

    if caRes == 0 :
        scRes = 0

    ########### CHECK FOR AVAIALABLE CELL PARAMETERS, SPACE GROUP and RESOLUTION

    if (mapfn !=  None or mtzfn !=None ):
        from stump import getCRYST , getRESO

        if (a == None or b == None or c == None or alpha == None or beta == None or gamma == None) :
            
            print "Getting cell paramters from coordinate file....."
            a,b,c,alpha , beta , gamma,d1  = getCRYST(pdbfile)

            if (a == None or b == None or c == None or alpha== None or beta==None or gamma == None ):
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
            print "setting Space Group to",sg

        if sg  in long2shortHM.keys():
            shortsg = long2shortHM[sg]
            sg = shortsg
        if sg not in sgtable.keys():
            print "Check --sg , Not recognised [%s][%d]"%( sg, len(sg))
            print "Check --sg , Not recognised [%s][%d]"%( sg, len(sg))
            import sys ; sys.exit()

        if resolution == None : 
            print "Getting resolution limit from coordinate file........"
            resolution = getRESO(pdbfile)
            if (resolution == None):
                print "Please input resolution " , resolution
                import sys ; sys.exit()
                
            print "Resolution = [ " , resolution, " ] "


        ############### IDENTIFY POOR REGIONS  ######################

        



    ############### CHECK RESTRAINT INFORMATION  ######################


    if caRes == 0 : 
        print "\n********** WARNING!! In ab initio loop building mode, No C-alpha positional restraints will be used"
        caRad = None; 
    if scRes == 0 : print "********** WARNING!! No sidechain centroid positional restraints will be used"; scRad = None


    
    if (caRes == 0 and poorOnly == 0 and (start == None or stop == None)):
        print "Need to set --use-ca-restraints True or specify [--loopseq --chainid --start --stop ] if you require loop building"
       	import sys ; sys.exit(0)



        
################## Checks for range of residues #######################################



    chainsFound  = getChains(pdbfile)    

    ## check if only a single character is given for the chain id
    if chainid != None and len(chainid) != 1 :
        print "Only one chain can be specified at a time"
        import sys ; sys.exit(1)

    ## check if given chain id is present in coordinate file
    if  chainid != None and  chainid not in chainsFound : 
        print "Chain id not found in pdb file", chainid , chainsFound
        import sys ;        sys.exit()

    ## if building ab initio loop, start, stop needs to be given. chain id is optional.
    if loopres != None :

        if start == None :
            print "Enter start residue for loop %s to be modelled"%loopres
            import sys ;            sys.exit()

        
        if stop == None :
            print "Enter stop residue for loop %s to be modelled"%loopres
            import sys ;            sys.exit()

        if chainid == None and len(chainsFound) != 1 :
            print "More than one chain detected , enter chain using --chainid "
            import sys ;        sys.exit()
            
        elif chainid == None and len(chainsFound) == 1 :
            print "Setting chain id to ",chainsFound[0]
            chainid = chainsFound[0]





    ## All chains
    if start == None and stop == None and loopres == None and poorOnly == 0 and chainid == None  :
	print "Rebuilding all chains"
	badresids = None
	if caRes == 0:
	    print "--use-ca-restraints to be set to True" ;
            import sys ; sys.exit()
            

    ## Single chains
    if start == None and stop == None  and poorOnly == 0 and chainid != None :


        prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        if chainid not in chainsFound :
            print "%s not found in coordinate file"%chainid
            import sys ; sys.exit()
            
        print "Rebuilding chain [%s] "%chainid
        if caRes == 0:
            print "--use-ca-restraints to be set to True" ;
            import sys ; sys.exit()

        for k , v in chids.items():
            if v == chainid :
                badresids.append(resids[k])


    ### If start and stop are given and no chain id then check coordinate file for options
    if chainid == None and len(chainsFound) != 1 and start != None and stop !=None :
        print "More than one chain detected , enter chain using --chainid "
        import sys ;        sys.exit()
            
    if chainid == None and len(chainsFound) == 1  and start != None and stop !=None :
        print "Setting chain id to ",chainsFound[0]
        chainid = chainsFound[0]


    ### if only start or stop is given irrespective of chain id then exits: ambiguos input
    if  poorOnly == 0 and ((start !=None and stop== None) or (start == None and stop != None))  :

        ## todo : check for sanity of start and stop
        if start == None :
            print "Enter start residue" 
            import sys ;            sys.exit()

        
        if stop == None :
            print "Enter stop residue"
            import sys ;            sys.exit()





    ### Get SEQUENCE FROM PDB FILE IF NOT SPECIFIED, (ab-initio mode)

    if start != None and stop !=None and  caRes == 0 and poorOnly == 0 and chainid != None and  loopres == None :
        from data import three2one
        if verbose > 7 :
            print "Reading loop seqeunce from coordinates file..."
            print "Loop to be modelled...",start,stop, chainid
        
        prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        loopseq = "" ; loopcounter = start

                
        for loopcounter in range(start,stop+1):
            aafound  = 0
            for k,v in resnums.items():
                if chids[k] != chainid :
                    continue
                if int(v) == loopcounter :
                    if resns[k] not in three2one.keys() :
                        print "Unrecognised amino acid type for residue ",loopcounter
                        continue
                    loopseq = loopseq + three2one[resns[k]]
                    aafound = 1
            if aafound ==  0 :
                print "Amino acid type for residue",loopcounter, "cannot be found from input coordinate file"
                print "Please input sequence for loop using --loopseq"
                import sys ;                    sys.exit()
        loopres = loopseq
        if verbose > 7 :            print "Loop sequence read from file", loopres


    if loopres != None and caRes == 0 and poorOnly == 0 :
	loopresids,opsaxoff,nullres,missingpts = parseLoopSequence(modelIn,loopres,start,stop,chainid,mconly,pdboutpath)
	badresids = loopresids 
        

    #### Non ab initio loop modelling
        
    if (poorOnly == 0 and (start != None and stop != None and chainid != None) and caRes == 1):
        loopresids = parseLoopSequenceSC(modelIn,loopres,start,stop,chainid,mconly,pdboutpath)
        badresids = loopresids
              
            

        


   ################   COPY MAP AND MTZ FILES  #########################################
    useomitmap = None

    if  mapfn == None and mtzfn == None :

        if poorOnly == 1 :
            print "In order to rebuild poor fitting regions --hklin or --mapin  needs to be specified"
            import sys ; sys.exit()
        print "********** WARNING!! no mtz or map file given, no electron density restraints will be used *************"

    if  mapfn != None:
        mapfilepath= mapfn
        if (os.path.isfile(mapfilepath)==False) :
            print "Cannot find file %s " %mapfn ; import sys ;            sys.exit()
        if  "map" not in mapfn : 
            print "Cannot understand i/p map format,the file should be  (*.map) format " ; import sys ; sys.exit()

        shutil.copyfile(mapfilepath, "%s/init.map" % (dir_xyzout))
        ccp4map = mapfn

    elif mtzfn != None :
        mtzfilepath= "%s"%(mtzfn)
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
            import sys ;              sys.exit()            
            

        if (f2label == None or phiclabel == None ):
            print "Column labels for  FC and PHIC  not specified, will use input coordinate structure to obtain FC and PHIC"
            sfall(pdbfile, mtzfn, "phased.mtz",resolution,usefreer,f1label,sigf1label,freeRlabel)
            mtzfn  = "phased.mtz"
            f2label  = "FC" ; phiclabel = "PHIC"



        ### MAKE UNWEIGHTED DIFF MAP : todo sigmaA weighted maps #########################

        from xray import fftnew
        ccp4map = "%s%dFo-%dFC.map"%(mtzfn,n,((n)-1))
        fcmap = "%sFC.map"%mtzfn
        mapfn = ccp4map

        if usefreer == 1 : 
            fftnew(mtzfn, ccp4map, pdbfile,f1label,f2label,phiclabel,sigf1label,n,(n)-1,freeRlabel)
            fftnew(mtzfn, fcmap, pdbfile ,f1label,f2label,phiclabel,sigf1label,0,1,freeRlabel,1)
        else :
            fftnew(mtzfn, ccp4map, pdbfile , f1label,f2label,phiclabel,sigf1label,n,(n)-1,None)
            fftnew(mtzfn, fcmap, pdbfile ,f1label,f2label,phiclabel,sigf1label,0,1,freeRlabel,1)




        useomitmap = 1
        if useomitmap == 1 :
            omit_mapfile = "%s%dFo-%dFC.omit.map"%(mtzfn,n,((n)-1))
            if usefreer == 1 : 
                omitmap(mtzfn, omit_mapfile, f1label,f2label,phiclabel,sigf1label,n,(n)-1,freeRlabel)
            else :
                omitmap(mtzfn, omit_mapfile, f1label,f2label,phiclabel,sigf1label,n,(n)-1)


        
        if poorOnly == 1 :
            from restraints import EDrestraint ; EDrestraint.setPenalty(2.) ; EDrestraint.setScatRad(1.) ;
            
            esmin, esmax, esmean, rcmult, xscoreCutoff = .000, 5., .0, 5, poorThreshold
            xscorer = XrayScorer(None, xscoreCutoff)
            mapcoeff = "%dF1-%dF2"%(n, (n)-1)
            badmcsc,id1 = xscorer.score(pdbfile, ccp4map, fcmap, f2label, phiclabel, mapcoeff, "mcsc",verbose)
            badmc,id2 = xscorer.score(pdbfile, ccp4map, fcmap, f2label, phiclabel, mapcoeff, "mc",verbose )
            badpept,id3 = xscorer.score(pdbfile, ccp4map, fcmap, f2label, phiclabel, mapcoeff, "pept",verbose )
            badsides,id4 = xscorer.score(pdbfile, ccp4map, fcmap, f2label, phiclabel, mapcoeff, "sc",verbose)

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
        esmin, esmax, esmean, rcmult, xscoreCutoff = .000, 5., .0, 5, poorThreshold

    multiPrepC = prepareChainV3.PrepareChain(rotLib)
    if mtzfn != None and useomitmap == 1:
        mapcoeff = "%dF1-%dF2"%(n, n-1)
        xrayRestGen.append( prepareChain.XrayRestraintsGenerator(omit_mapfile, "map", f2label, phiclabel, mapcoeff, esmin, esmax, esmean, [], ) )

    elif mtzfn != None or mapfn !=None :
        mapcoeff = "%dF1-%dF2"%(n, n-1)
        xrayRestGen.append( prepareChain.XrayRestraintsGenerator(ccp4map, "map", f2label, phiclabel, mapcoeff, esmin, esmax, esmean, [], ) )



    ############### GENERATE MODELS  ######################

    for  i in range(nmodels):
        nb = 0 
        print
        print

        
        if ((poorOnly == 1 and len(badresids) > 0) or (poorOnly == 0 and badresids == None) or (poorOnly == 0 and len(badresids) > 0 )):
            ml = Multiloop(modelIn, badresids, mconly, caRad, scRad, scvdwr, guidedsampling, popsize, backtrack, 1 , "pre."+str(i)+"."+rtkmodel, xrayRestGen, multiPrepC,sconly, cacaCutoff,nullres,modelN2C,natt,missingpts) 
            
            ml.ranker = None 
            if poorOnly == 1 :
                ml.poor = "poor"
            if start != None and stop != None and chainid != None :
                ml.loopstartnum = start
                ml.loopstopnum = stop

            if mtzfn != None  and useomitmap == 1 :
                ml.ranker = XrayRanker(omit_mapfile, "map", f2label, phiclabel, mapcoeff, esmin, esmax)
                ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
                ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble = None
                ml.cellsym = [ a, b, c, alpha, beta, gamma, sgtable[sg][0] ]

            elif mtzfn != None or mapfn !=None:
                ml.ranker = XrayRanker(ccp4map, "map", f2label, phiclabel, mapcoeff, esmin, esmax)
                ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
                ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble = None
                ml.cellsym = [ a, b, c, alpha, beta, gamma, sgtable[sg][0] ]

                


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
            
                
        else :
            
            shutil.copy(modelIn, "pre."+str(i)+"."+rtkmodel)

        if sconly == 1 and opsax == 0 and mapfn != None :
            print "** In order to rebuild sidechains set --opsax True ** "
            import sys ; sys.exit()


        if sconly == 1 and opsax == 1 and mapfn == None :
            print "** In order to rebuild sidechains using opsax MTZ/MAP file required ** "
            import sys ; sys.exit()

        if sconly == 1 and opsax == 0 and mapfn == None :
            print "** In order to rebuild sidechains for given backbone :  "
            print "** (1) Give MTZ/MAP file  ** "
            print "                 or          "
            print "** (2) Generate all atom model using --ca-restraint-radius 0.5 --use-ca-restraints True **  "

            
            import sys ; sys.exit()            


        if sconly == 1 and opsax == 1 and mapfn != None and opsaxoff == 1 :
            print "Cannot use opsax for proteins > 200 aminoacids, try  rebuilding all atoms :"
            print "--sconly False --ca-restraint-radius 0.5 --use-ca-restraints True "
            import sys ; sys.exit()

            
            
        if opsax == 0 or mapfn == None or mconly == 1 or opsaxoff == 1:
            os.rename("pre."+str(i)+"."+rtkmodel, str(i)+"."+rtkmodel )
            print "Renaming", "pre."+str(i)+"."+rtkmodel, " to ",str(i)+"."+rtkmodel 
        else :
            print "\n\nStarting OPSAX ........"
            if poorOnly == 1 and sconly == 0 : 
                if len(badresids) == 0 :
                    badresids = list(  set( badscs )  )
                if len(unbuilt[0]) == 0  :
                    badresids = list(  set(  findChangedSC(modelIn, "pre."+str(i)+"."+rtkmodel) + badscs )  )
                else :
                    badresids = list(  set( badscs )  )
                    
            elif (poorOnly) == 1 and sconly == 1 : 
                badresids = badscs 

            elif (poorOnly) == 0 and sconly == 1 : 
                badresids = badresids 

            elif   poorOnly == 0 and len(unbuilt[0]) == 0  :
                badresids = list(  set(  findChangedSC(modelIn, "pre."+str(i)+"."+rtkmodel) + badscs )  )
            else:
                badresids = badresids
                ### todo remove missing sidechains ??
            
            
            if len(badresids) > 250:
                print "Length of section > 250 residues, can not perform OPSAX"
                os.rename("pre."+str(i)+"."+rtkmodel , str(i)+"."+rtkmodel)
            else :
                if  sconly == 1 or nb !=1 :
                    scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 0, 1
                else :
                    scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 1, 1


                if useomitmap == 1 and mtzfn != None : 
                    nb2 = SCplacement("pre."+str(i)+"."+rtkmodel, scReduction, str(i)+"."+rtkmodel, "dotfile", useDEE, omit_mapfile, f1label, f2label , phiclabel, mapcoeff, esmin, esmax, None, useGivenRot, badresids, scPrepC).run()

                else : 
                    nb2 = SCplacement("pre."+str(i)+"."+rtkmodel, scReduction, str(i)+"."+rtkmodel, "dotfile", useDEE, mapfn, f1label, f2label , phiclabel, mapcoeff, esmin, esmax, None, useGivenRot, badresids, scPrepC).run()



                if nb2 != 1:
                    shutil.copy("pre."+str(i)+"."+rtkmodel ,str(i)+"."+rtkmodel)

            print "\n Finished OPSAX ........"
        adjustBfac2(str(i)+"."+rtkmodel, modelIn,nativeBfac, mcBfac, scBfac)
        print 

        #copyNonprotein(pdbfile, str(i)+"."+rtkmodel) ## automatic water addition can go here
        print 
        #copyWater(pdbfile, str(i)+"."+rtkmodel) ## automatic water addition can go here
        print
        removeMODEL(str(i)+"."+rtkmodel) ## automatic water addition can go here
        print

        print "------------------------------------ MODEL ASSESSMENT --------------------------------"
        
        if mtzfn != None :
            from prefRapperV4WithRefmac import refmacNew
            from protRefine_nick import molProbity_badres, add_cryst_card 
            from copyheader import copyheader

            sfcheck(mtzfn,pdbfile,f1label,sigf1label,freeRlabel ) 
            os.rename("sfcheck.log","sfcheck.init.log")

            sfall(str(i)+"."+rtkmodel, mtzfn, "phased."+str(i)+".mtz",resolution,usefreer,f1label,sigf1label,freeRlabel ) 
            sfcheck("phased."+str(i)+".mtz",str(i)+"."+rtkmodel,f1label,sigf1label,freeRlabel ) 
            os.rename("sfcheck.log","sfcheck."+str(i)+".log")

            
            sfall("pre."+str(i)+"."+rtkmodel, mtzfn, "phased.pre."+str(i)+".mtz",resolution,usefreer,f1label,sigf1label,freeRlabel ) 
            sfcheck("phased.pre."+str(i)+".mtz","pre."+str(i)+"."+rtkmodel,f1label,sigf1label,freeRlabel ) 
            os.rename("sfcheck.log","sfcheck.pre."+str(i)+".log")


            print
            print
            print "R-factor and correlation for initial coordinates"
            parse_sfcheck_logfile("sfcheck.init.log")
            print

            print
            print
            print "R-factor and correlation for output (rappertk) coordinates before opsax"
            parse_sfcheck_logfile("sfcheck.pre."+str(i)+".log")
            copyheader("pre."+str(i)+"."+rtkmodel,"modelinit.pdb")
                       
            print
            print "R-factor and correlation for output (rappertk) coordinates"
            parse_sfcheck_logfile("sfcheck."+str(i)+".log")
            copyheader(str(i)+"."+rtkmodel,"modelinit.pdb")


            ############### REFMAC #############################################################################
            
            #add_cryst_card(str(i)+"."+rtkmodel, a, b, c, alpha, beta, gamma, sg)
            #refmacNew("phased."+str(i)+".mtz", "refmac."+str(i)+".mtz", str(i)+"."+rtkmodel, str(i)+".refmac.pdb" )
            # sfcheck("refmac."+str(i)+".mtz",str(i)+".refmac.pdb",f1label,sigf1label,freeRlabel )
            # os.rename("sfcheck.log","sfcheck.refmac."+str(i)+".log")
            # print
            # print "R-factor and correlation for refined (refmac) coordinates"
            # parse_sfcheck_logfile("sfcheck.refmac."+str(i)+".log")
            
            #refmacNew(mtzfn, "refmac.init.mtz", pdbfile, "refmac.init.pdb" )
            #sfcheck("refmac.init.mtz","refmac.init.pdb",f1label,sigf1label,freeRlabel )
            #os.rename("sfcheck.log","sfcheck.refmac.init.log")
            #print
            #print "R-factor and correlation for refined (refmac) coordinates"
            #parse_sfcheck_logfile("sfcheck.refmac.init.log")

            #if (os.path.isfile(pdbout)==False) :
            #    print "Refmac cycle %d failed, no %s generated"%(rcycle,pdbout)
            #    import sys ; sys.exit()
            #if (os.path.isfile(hklout)==False) :
            #    print "Refmac cycle %d failed, no %s generated"%(rcycle,hklout)
            #    import sys ; sys.exit()

            ############### REFMAC ENDS #############################################################################

        elif mapfn != None :
            from prefRapperV4WithRefmac import refmacNew
            from protRefine_nick import molProbity_badres, add_cryst_card 
            from copyheader import copyheader

            sfcheck(mapfn,pdbfile,None,None,None)
            os.rename("sfcheck.log","sfcheck.init.log")

            sfcheck(mapfn, str(i)+"."+rtkmodel,None,None,None)
            os.rename("sfcheck.log","sfcheck."+str(i)+".log")

            print
            print
            print "R-factor and correlation for initial coordinates"
            parse_sfcheck_logfile("sfcheck.init.log")
            print
            print
            print "R-factor and correlation for output (rappertk) coordinates"
            parse_sfcheck_logfile("sfcheck."+str(i)+".log")

            copyheader(str(i)+"."+rtkmodel,"modelinit.pdb")


            ############### REFMAC #############################################################################
            
            #add_cryst_card(str(i)+"."+rtkmodel, a, b, c, alpha, beta, gamma, sg)
            #refmacNew("phased."+str(i)+".mtz", "refmac."+str(i)+".mtz", str(i)+"."+rtkmodel, str(i)+".refmac.pdb" )
            # sfcheck("refmac."+str(i)+".mtz",str(i)+".refmac.pdb",f1label,sigf1label,freeRlabel )
            # os.rename("sfcheck.log","sfcheck.refmac."+str(i)+".log")
            # print
            # print "R-factor and correlation for refined (refmac) coordinates"
            # parse_sfcheck_logfile("sfcheck.refmac."+str(i)+".log")
            
            #refmacNew(mtzfn, "refmac.init.mtz", pdbfile, "refmac.init.pdb" )
            #sfcheck("refmac.init.mtz","refmac.init.pdb",f1label,sigf1label,freeRlabel )
            #os.rename("sfcheck.log","sfcheck.refmac.init.log")
            #print
            #print "R-factor and correlation for refined (refmac) coordinates"
            #parse_sfcheck_logfile("sfcheck.refmac.init.log")

            #if (os.path.isfile(pdbout)==False) :
            #    print "Refmac cycle %d failed, no %s generated"%(rcycle,pdbout)
            #    import sys ; sys.exit()
            #if (os.path.isfile(hklout)==False) :
            #    print "Refmac cycle %d failed, no %s generated"%(rcycle,hklout)
            #    import sys ; sys.exit()

            ############### REFMAC ENDS #############################################################################

        print "------------------------------------ MODEL ASSESSMENT ENDS ---------------------------"




def parse_sfcheck_logfile(filename):
    import re
    lines = open(filename, 'r').readlines()
#    pattern = re.compile('(\d+.\d+)')
    for l in lines :
#        print l
        if 'R-factor           ' in l :
            print l
            #params = pattern.match(l)
            #if params != None:
            #    rpar = list(params.groups())
            #    rfactor = rpar[1])

        if 'Correlation factor ' in l :
            print l
        
    
    
        if " Rfree,Rrest" in l :
            print l






def parseLoopSequence(pdbfile,loopres,start,stop,chainid,mconly,pdbout) : 

    null = [] ; missingpts = []
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts , atomids , bfac = prot2res.readProtRes2(prot)
    seqlen = len(loopres) ; sslen = stop - start +1
    rev_resnum = {} ; badresids = [] ; chainkey = []
    missing = 0
    for k,v in chids.items():
        if chids[k] == chainid:
            chainkey.append(k)

    chainkey.sort()
    if (seqlen != sslen):
        print "Sequence length does not match start and stop points"
        import sys ; sys.exit()

    for k in chainkey:
        rev_resnum[int(resnums[k])] = k
        
    for index in range(start,stop+1):
        if index in rev_resnum.keys():
            
            revindex = rev_resnum[index]
            if not isAAres(resns[revindex]) :
                continue
            if ' CA ' in res[revindex].keys():
                badresids.append(resids[revindex])
            else :
                missing = 1
        else :
            missing = 1
    if missing  == 0 :
        return badresids,0,null,missingpts


    from prepareChainV3 import  incompleteREScorrection2

    nres, nresids, nresnums, nresns, nchids, ninscodes, npts , null,missingpts  = incompleteREScorrection2(res, resids, resnums, resns, chids, inscodes, pts,start,stop,chainid,loopres)
    res, resids, resnums, resns, chids, inscodes, pts = {}, {}, {}, {}, {}, {}, []
    res, resids, resnums, resns, chids, inscodes, pts = nres, nresids, nresnums, nresns, nchids, ninscodes, npts
    os.remove(pdbfile)

    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [],pdbfile)
    modelRenderer.render(pts)

    
    rev_resnum = {} ; badresids = [] ; chainkey = [];    missing = 0
    for k,v in chids.items():
        if chids[k] == chainid:
            chainkey.append(k)

    chainkey.sort()

    for k in chainkey:
        rev_resnum[int(resnums[k])] = k
        
    for index in range(start,stop+1):
        if index in rev_resnum.keys():
            revindex = rev_resnum[index]
            badresids.append(resids[revindex])
        else :
            missing = 1


    if missing  == 0 :
        return badresids,0,null , missingpts
    
    else:
        print "Cant build.Error 201, please report problem to ak459@cam.ac.uk"
        import sys ; sys.exit(1)





def parseLoopSequenceSC(pdbfile,loopres,start,stop,chainid,mconly,pdbout) : 
    from data import resAtoms
    from pdbr import isAAres
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    rev_resnum = {} ; badresids = []
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
    
    

    import optparse
    from commonOptions import makeParser, parseOptions
    #parser = makeParser()
    parser = optparse.OptionParser()
    ################# I/O Files and directories ##########################################
    
    parser.add_option("--xyzin", action='store', type='string', dest='pdbfile', help='Name of input coordinate file')
    parser.add_option("--xyzout", action='store', type='string', dest='pdbout', help='Name of output coordinate file',default="modelout.pdb")
    parser.add_option("--dir-xyzout", action='store', type='string', dest='dir_xyzout', help='Name of directory in which to create output files. complete PATH needed')
    parser.add_option("--mapin", action='store', type='string', dest='mapfn', help='Name of input MAP file',default=None)
    parser.add_option("--hklin", action='store', type='string', dest='mtzfn', help='Name of input (phased) MTZ file',default=None)


    ### Assignment of restraints and radii #############################################
    
    parser.add_option("--use-ca-restraints", action='store', dest='caRes', help='[True/False], Apply positional restraints on the C-alpha atoms',default="False")
    parser.add_option("--use-sc-restraints", action='store', dest='scRes',type= 'string', help='[True/False],  Apply positional restraints on the centroid of the sidechain atoms',default="False")

    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='Radius of spherical restraint ( angstrom ) on the C-alpha atom position', default=1.0)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='Radius of spherical restraint ( angstrom ) on the centroid of the sidechain atoms', default=2.0)

    ############ General parameters #####################################################
    
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='Factor by which to reduce effective Van der Waals distance for sidechain atoms', default=0.75)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='Population size', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='Level of printed output [1-10], 1:Concise output log, 10: Detailed output log', default=7)
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
    
    parser.add_option("--FP", action='store', type='string', dest='f1label', help='Column label for FP in MTZ file', default=None)
    parser.add_option("--SIGFP", action='store', type='string', dest='sigf1label', help='Column label for sigFP in MTZ file', default=None)
    parser.add_option("--FC", action='store', type='string', dest='f2label', help='Column label for FC in MTZ file', default=None)
    parser.add_option("--PHIC", action='store', type='string', dest='phiclabel', help='Column label for PHIC in MTZ file', default=None)
    parser.add_option("--use-FreeR", action='store', type='string', dest='usefreer', help='[True/False] Use FreeR set ? ', default="False")
    parser.add_option("--FreeR", action='store', type='string', dest='freeRlabel', help='Column label for FreeR in MTZ file', default=None)
    parser.add_option("--n", action='store', type='int', dest='n', help='Value of n for difference map calculations nFo-(n-1)Fc', default=2)


    ############# Residues to be modelled ####################################
    

    parser.add_option("--rebuild-poor-regions-only", action='store', type='string', dest='poorOnly', help='[True/False] Rebuild regions ofinput structure with poor fit to an electron density map. Residues to be rebuilt are identified using a real space correlation coefficientscore, the cut-off for which is set using --poor-fit-threshold.', default="False")
    parser.add_option("--poor-fit-threshold", action='store', type='float', dest='poorThreshold', help='Correlation coefficient threshold to identify poor fitting regions', default=0.9)



    parser.add_option("--loopseq", action='store', type='string', dest='loopres', help='Amino acid sequence for loop to be built', default=None)
    parser.add_option("--start", action='store', type='int', dest='start', help='Residue number to start building from ', default=None)
    parser.add_option("--stop", action='store', type='int', dest='stop', help='Residue number to stop building at', default=None)
    parser.add_option("--chainid", action='store', type='string', dest='chainid', help='Chain ID of section to be built.', default=None)
    parser.add_option("--modelN2C", action='store', type='string', dest='modelN2C', help='[True/False] Model fragment without loop closure restraints. Used in conjunction with --start, --stop, --chainid. Requires --use-ca-restraints True ', default="False")

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


    options = parseOptions(parser)


    main(options.pdbfile,options.pdbout,options.dir_xyzout,options.mapfn,options.mtzfn,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,options.nmodels,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freeRlabel,options.n,options.poorOnly,options.poorThreshold,options.loopres,options.start,options.stop,options.chainid,options.modelN2C,options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt)

if __name__ == "__main__" :
    callmain()
