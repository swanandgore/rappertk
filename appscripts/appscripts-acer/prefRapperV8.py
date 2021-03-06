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
from loopbuildV5 import Multiloop, incompleteSCcorrection2
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


cnsArgs = {}
for cycle in range(20) : cnsArgs[cycle] = {} ; cnsArgs[cycle]["num_cycles"] = 2 ; cnsArgs[cycle]["temperature"] = 5000
#cnsArgs[0]["wa"] = -1 ; cnsArgs[0]["num_cycles"] = 1 ; cnsArgs[0]["temperature"] = 50
#cnsArgs[1]["wa"] = -1 ; cnsArgs[1]["num_cycles"] = 1 ; cnsArgs[1]["temperature"] = 50

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
## OT1 to O
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
    print "Copy Bfactors from", refpdb, "to", pdbfile, "if coordinates are same, else set to 30"
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
    print "Copy Bfactors from", refpdb, "to", pdbfile
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



def main() :

    import optparse ; parser = optparse.OptionParser()
    

    ################# I/O Files and directories ##########################################
    
    parser.add_option("--xyzin", action='store', type='string', dest='pdbfile', help='Name of coordinate input structure')
    parser.add_option("--xyzout", action='store', type='string', dest='pdbout', help='Name of PDB output file',default="modelout.pdb")
    parser.add_option("--dir-xyzout", action='store', type='string', dest='dir_xyzout', help='Name of directory to create all the files during refinement. complete PATH needed')
    parser.add_option("--mapin", action='store', type='string', dest='mapfn', help='Input map file',default=None)
    parser.add_option("--hklin", action='store', type='string', dest='mtzfn', help='Input phased MTZ file',default=None)


    ### Assignment of restraints and radii #############################################
    
    parser.add_option("--use-ca-restraints", action='store', type='string', dest='caRes', help='[True/False], Apply positional restraints on the C-alpha atom',default="True")
    parser.add_option("--use-sc-restraints", action='store', type='string', dest='scRes', help='[True/False],  Apply positional restraints on the centroid of the sidechain atoms',default="True")

    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='Radius of spherical restraint on the C-alpha atom position', default=1.0)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='Radius of spherical restraint on the centroid of the sidechain atoms', default=2.0)

    ############ General parameters #####################################################
    
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='Factor by which to reduce effective Van der Waals distance for sidechain atoms', default=0.75)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='Population size', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='Level of printed output [0-10], 0:Concise output log, 10: Detailed output log', default=0)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='[numsteps]X[stepsize] e.g. 4X5 will set backtrack to numsteps and stepsize to 4 and 5 respectively. For detailed help see Rappertk wiki/webpage', default= None)
    parser.add_option("--rotamerlib", action='store', type='string', dest='rotLib', help='Rotamer library to use when building side chains [PRL/SCL1.0/SCL0.5/SCL0.2] ', default='PRL')        
    parser.add_option("--num-models", action='store', type='int', dest='nmodels', help='Number of models wanted ', default=1)


    #################### Build parameters ################################################
    
    parser.add_option("--mconly", action='store', type='string', dest='mconly', help='[True/False] Build mainchain only', default="False")
    parser.add_option("--sconly", action='store', type='string', dest='sconly', help='[True/False] Build side chains only', default="False")
    parser.add_option("--opsax", action='store', type='string', dest='opsax', help='[True/False] Reassign side chains with OPSAX, can only be True when MTZ or MAP file is given', default="False")
    parser.add_option("--attempts", action='store', type='int', dest='natt', help='Number of attempts to build a band', default=5)
    parser.add_option("--detect-breaks", action='store', type='string', dest='detectBreaks', help='[True/False] Detect chain breaks?', default= "True")
    parser.add_option("--cacaCutoff", action='store', type='float', dest='cacaCutoff', help='Minimum distance between adjacent Calpha atoms in order to detect a chain-break', default=5.)


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
    parser.add_option("--use-freer", action='store', type='string', dest='usefreer', help='[True/False] Use Free R set ? ', default="False")
    parser.add_option("--FREER", action='store', type='string', dest='freeRlabel', help='Column label for FreeR in MTZ file', default=None)
    parser.add_option("--n", action='store', type='int', dest='n', help='Value of n for difference map calculations nFo-(n-1)Fc', default=2)


    ############# Residues to be modelled ####################################
    

    parser.add_option("--rebuild-poor-regions-only", action='store', type='string', dest='poorOnly', help='[True/False] Rebuild regions of poor fit to an electron density map. Residues to be rebuilt are identified using a real space scoring function, the cut-off for which is set using --poor-fit-threshold.', default="False")
    parser.add_option("--poor-fit-threshold", action='store', type='float', dest='poorThreshold', help='Correlation coefficient threshold to identify ill fitting regions', default=0.9)



    parser.add_option("--loopseq", action='store', type='string', dest='loopres', help='Amino acid sequence for loop to be built', default=None)
    parser.add_option("--start", action='store', type='int', dest='start', help='Residue number to start building from ', default=None)
    parser.add_option("--stop", action='store', type='int', dest='stop', help='Residue number to stop building ', default=None)
    parser.add_option("--chainid", action='store', type='string', dest='chainid', help='Chain ID of section to be built.', default=None)
    parser.add_option("--modelN2C", action='store', type='string', dest='modelN2C', help='[True/False] Model fragment without loop closure restraints. Used in conjunction with --start, --stop, --chainid. Requires --use-ca-restraints True ', default="False")

    ######### Ouptut parameters #############################################

    parser.add_option("--models-get-native-bfactors", action='store', type='string', dest='nativeBfac', help='[True/False] Assign B-factors of remodelled atoms to original values', default="False")
    parser.add_option("--default-mainchain-b-factor", action='store', type='float', dest='mcBfac', help='The B-factor assigned to the newly built main chain atoms', default=20.)
    parser.add_option("--default-sidechain-b-factor", action='store', type='float', dest='scBfac', help='The B-factor assigned to the newly built side chain atoms', default=30.)



    ### Electron density parametets #########################################

    parser.add_option("--minimum-sig", action='store', type='float', dest='minXSig', help='Minimum sigma ', default=0.25)
    parser.add_option("--maximum-sig", action='store', type='float', dest='maxXSig', help='Maximum sigma ', default=2.0)



    ########## Optional restraints ##########################################
    
    parser.add_option("--make-ed-optional", action='store', type='int', dest='edOpt', help='[0/1] 1:Yes 0:No. If 1 then positive density mainchain restraint will be made optional. If 0, then the main chain will be unconditionally forced to lie in positive density. This is primarily useful when tracing through a structure with regions in very poor (non-existent) density', default=0)
    parser.add_option("--make-all-restraints-optional", action='store', type='int', dest='allOpt', help='[0/1] 1:Yes 0:No. If 1, then all  restraints will be made optional', default=0)    



    
    

    (options, args) = parser.parse_args()
    

    print "========================== Now  printing all KEYWORD values ==========================="




    print
    print "                           Directories and Files"
    print
    
    print "	--xyzin ".ljust(50), 
    print "%30s "%options.pdbfile 
    print "	--xyzout ".ljust(50), 
    print "%30s "%options.pdbout 
    print "	--dir-xyzout ".ljust(50), 
    print "%30s "%options.dir_xyzout 
    print "	--mapin ".ljust(50), 
    print "%30s "%options.mapfn 
    print "	--hklin ".ljust(50), 
    print "%30s "%options.mtzfn 
    print "	--use-ca-restraints ".ljust(50), 

    print
    print "                           Restraint paramaters "
    print

    print "%20s "%options.caRes 
    print "	--use-sc-restraints ".ljust(50), 
    print "%20s "%options.scRes 
    print "	--ca-restraint-radius ".ljust(50), 
    print "%20s "%options.caRad 
    print "	--sc-centroid-restraint-radius ".ljust(50), 
    print "%20s "%options.scRad





    print
    print "                           Build paramaters "
    print
    print "	--sidechain-vdw-reduction ".ljust(50), 
    print "%20s "%options.scReduction 
    print "	--population-size ".ljust(50), 
    print "%20s "%options.popsize 
    print "	--verbose ".ljust(50), 
    print "%20s "%options.verbose 
    print "	--backtrack ".ljust(50), 
    print "%20s "%options.backtrack 
    print "	--mconly ".ljust(50), 
    print "%20s "%options.mconly 
    print "	--sconly ".ljust(50), 
    print "%20s "%options.sconly 
    print "	--make-all-restraints-optional ".ljust(50), 
    print "%20s "%options.allOpt 


    print "	--loopseq ".ljust(50), 
    print "%20s "%options.loopres 
    print "	--start ".ljust(50), 
    print "%20s "%options.start 
    print "	--stop ".ljust(50), 
    print "%20s "%options.stop 
    print "	--chainid ".ljust(50), 
    print "%20s "%options.chainid 


    print
    print "                           Crystallographic paramaters "
    print
    print "	--opsax ".ljust(50), 
    print "%20s "%options.opsax 
    print "	--use-freer ".ljust(50), 
    print "%20s "%options.usefreer 
    print "	--rotamerlib ".ljust(50), 
    print "%20s "%options.rotLib 
    print "	--num-models ".ljust(50), 
    print "%20s "%options.nmodels 
    print "	--a ".ljust(50), 
    print "%20s "%options.a 
    print "	--b ".ljust(50), 
    print "%20s "%options.b 
    print "	--c ".ljust(50), 
    print "%20s "%options.c 
    print "	--alpha ".ljust(50), 
    print "%20s "%options.alpha 
    print "	--beta ".ljust(50), 
    print "%20s "%options.beta 
    print "	--gamma ".ljust(50), 
    print "%20s "%options.gamma 
    print "	--sg ".ljust(50), 
    print "%20s "%options.sg 
    print "	--resolution ".ljust(50), 
    print "%20s "%options.resolution 
    print "	--FP ".ljust(50), 
    print "%20s "%options.f1label 
    print "	--SIGFP ".ljust(50), 
    print "%20s "%options.sigf1label 
    print "	--FC ".ljust(50), 
    print "%20s "%options.f2label 
    print "	--PHIC ".ljust(50), 
    print "%20s "%options.phiclabel 
    print "	--FREER ".ljust(50), 
    print "%20s "%options.freeRlabel 


    print "	--n ".ljust(50), 
    print "%20s "%options.n 
    print "	--rebuild-poor-regions-only ".ljust(50), 
    print "%20s "%options.poorOnly 
    print "	--poor-fit-threshold ".ljust(50), 
    print "%20s "%options.poorThreshold 
    print "	--default-mainchain-b-factor ".ljust(50), 
    print "%20s "%options.mcBfac 
    print "	--default-sidechain-b-factor ".ljust(50), 
    print "%20s "%options.scBfac 
    print "	--models-get-native-bfactors ".ljust(50), 
    print "%20s "%options.nativeBfac 
    print "	--minimum-sig ".ljust(50), 
    print "%20s "%options.minXSig 
    print "	--maximum-sig ".ljust(50), 
    print "%20s "%options.maxXSig 
    print "	--cacaCutoff ".ljust(50), 
    print "%20s "%options.cacaCutoff 
    print "	--detect-breaks ".ljust(50), 
    print "%20s "%options.detectBreaks 
    print "	--modelN2C ".ljust(50), 
    print "%20s "%options.modelN2C 
    print "	--make-ed-optional ".ljust(50), 
    print "%20s "%options.edOpt 

    
 
    print
    print
    print "========================== End KEYWORD values =========================================="
    print
    print
    

    misc.setVerbosity(options.verbose) ; randomize(1)


    ######### Check for / make out directory #############
    
    if  os.path.isfile(options.dir_xyzout) :
        print options.dir_xyzout,"Is a file"
        print "Please rename output directory"
        import sys ; sys.exit()

    head,tail  = os.path.split(options.dir_xyzout)
    if head and os.path.isdir(head):
        print "Checking for directory ............",head
    else:
        print "%s does not exists"%head


    if not os.path.isdir(options.dir_xyzout) :
        os.mkdir(options.dir_xyzout)

    os.chdir(options.dir_xyzout)


    #### Copy coordinate file  ####

    pdbfilepath = "%s"%(options.pdbfile)
    if (os.path.isfile(pdbfilepath)==False) :
        print "Cannot find file %s " %options.pdbfile ; import sys ;        sys.exit()
    shutil.copyfile(pdbfilepath, "%s/%s" % (options.dir_xyzout,"modelinit.pdb"))




    ###########    SANITY CHECKS                      #############################


    if options.modelN2C == 1 and options.caRes == 0 :
        print "When building as a fragment (--modelN2C 1), --use-ca-restraints needs to be set to 1"
        import sys ; sys.exit()
        
    if options.mconly not in ["True", "False"]:
        print "--mconly needs to be either True or False" , options.mconly
        import sys ; sys.exit()

    if (options.scReduction < 0. or options.scReduction > 1.) and options.mconly!= 1 :
        print "--sidechain-vdw-reduction needs to be between 0 and 1." ;
        import sys ; sys.exit()
        
    if options.rotLib not in ["PRL", "SCL0.5","SCL1.0","SCL0.2"]:
        print "Unrecognised rotamer library, options are PRL, SCL0.5, SCL1.0, SCL0.2"




    if options.scRes not in [0, 1]:
        print "--use-sc-restraints needs to be either 0 or 1, Incorrect value input" ,options.scRes
        import sys ; sys.exit()
    if options.caRes not in [0, 1]:
        print "--use-ca-restraints needs to be either 0 or 1 , Incorrect value input",options.caRes
        import sys ; sys.exit()
    if options.nativeBfac not in [0,1]:
        print "--models-get-native-bfactors needs to be either 0 or 1,  Incorrect value input",options.nativeBfac
        import sys ; sys.exit()

    #if options.mcBfac != None or options.scBfac != None :
    #    options.nativeBfac = 0

    if options.opsax not in [0,1]:
        print "--opsax needs to be either 0 or 1,  Incorrect value input",options.opsax
        import sys ; sys.exit()
    if options.edOpt  not in [0,1]:
        print "--make-ed-optional needs to be either 0 or 1.  Incorrect value input",options.edOpt
        import sys ; sys.exit()
    if options.allOpt not in [0,1] :
        print "--make-all-restraints-optional needs to be either 0 or 1.  Incorrect value input",options.allOpt
        import sys ; sys.exit()
    if options.dir_xyzout == None :
        print "Directory for output files needs to be set. --dir_xyzout needs to be set",options.dir_xyzout
        import sys ; sys.exit()

    if options.backtrack == "None":
        options.backtrack = None
    else :
        if options.backtrack != None : ## this is supposed to be aXb, numstepsXstepsize, eg 4X5
            xpos = options.backtrack.find('X')
            if xpos <= 0:
                print "Error in specifying backtrack, it should be of the form of numstepsXstepsize e.g. 4X5  where numsteps and stepsize are integer values"
                print  "will set backtrack to numsteps and stepsize to 4,5 respectively. For detailed help see Rappertk wiki/webpage"
                import sys ; sys.exit()
            #if type(1) != int(options.backtrack[0:xpos]) or type(1) != int(options.backtrack[xpos+1:]) :
            #    print "Error in specifying backtrack, it should be of the form of numstepsXstepsize e.g. 4X5 "
            #    print "Please use integer values for numsteps and stepsize"
            #    import sys ; sys.exit()
                

    if options.mconly == 0 :
        options.mconly = None

    if options.allOpt == 1:
        options.edOpt = 1


    ############### CHECK RESTRAINT INFORMATION  ######################


    if options.caRes == 0 : 
        print "\n********** WARNING!! In ab initio loop building mode, No C-alpha positional restraints will be used"
        options.caRad = None; 
    if options.scRes == 0 : print "********** WARNING!! No sidechain centroid positional restraints will be used"; options.scRad = None


    
    if (options.caRes == 0 and options.poorOnly == 0 and (options.start == None or options.stop == None or options.chainid ==None)):
        print "Need to set --use-ca-restraints 1 or specify [--loopseq --chainid --start --stop ] if you require loop building"
       	import sys ; sys.exit(0)


    pdboutpath = "%s/%s"%(options.dir_xyzout,options.pdbout)





################## Checks for range of residues #######################################

    badresids = [] ;         nullres = [] ; missingpts = []

    chainsFound  = getChains(options.pdbfile)    

    ## check if only a single character is given for the chain id
    if options.chainid != None and len(options.chainid) != 1 :
        print "Only one chain can be specified at a time"
        import sys ; sys.exit(1)

    ## check if given chain id is present in coordinate file
    if  options.chainid != None and  options.chainid not in chainsFound : 
        print "Chain id not found in pdb file", options.chainid , chainsFound
        import sys ;        sys.exit()

    ## if building ab initio loop, start, stop needs to be given. chain id is optional.
    if options.loopres != None :
        ## todo : check for sanity of start and stop
        if options.start == None :
            print "Enter start residue for loop %s to be modelled"%options.loopres
            import sys ;            sys.exit()

        
        if options.stop == None :
            print "Enter stop residue for loop %s to be modelled"%options.loopres
            import sys ;            sys.exit()


        if options.chainid == None and len(chainsFound) != 1 :
            print "More than one chain detected , enter chain using --chainid "
            import sys ;        sys.exit()
            
        elif options.chainid == None and len(chainsFound) == 1 :
            print "Setting chain id to ",chainsFound[0]
            options.chainid = chainsFound[0]

    ## All chains
    if options.start == None and options.stop == None and options.loopres == None and options.poorOnly == 0 and options.chainid == None  :
	print "Rebuilding all chains"
	badresids = None
	if options.caRes == 0:
	    print "--use-ca-restraints to be set to 1" ;
            import sys ; sys.exit()
            

    ## Single chains
    if options.start == None and options.stop == None  and options.poorOnly == 0 and options.chainid != None :


        prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        if options.chainid not in chainsFound :
            print "%s not found in coordinate file"%options.chainid
            import sys ; sys.exit()
            
            print "Rebuilding chain [%s] "%options.chainid
            if options.caRes == 0:
                print "--use-ca-restraints to be set to 1" ;
                import sys ; sys.exit()

            for k , v in chids.items():
                if v == options.chainid :
                    badresids.append(resids[k])


    ### if only start or stop is given irrespective of chain id then exits: ambiguos input
    if options.loopres == None and options.poorOnly == 0 and ((options.start !=None and options.stop== None) or (options.start == None and options.stop != None))  :

        ## todo : check for sanity of start and stop
        if options.start == None :
            print "Enter start residue for loop %s to be modelled"%options.loopres
            import sys ;            sys.exit()

        
        if options.stop == None :
            print "Enter stop residue for loop %s to be modelled"%options.loopres
            import sys ;            sys.exit()






    ### Get SEQUENCE FROM PDB FILE IF NOT SPECIFIED

    if options.start != None and options.stop !=None and  options.caRes == 0 and options.poorOnly == 0 and options.chainid != None and  options.loopres == None :

        prot = protein(pdbfilepath, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        loopseq = ""
        from data import three2one

        if options.verbose > 7 :
            print "Reading loop seqeunce from coordinates file..."
            print "Loop to be modelled...",options.start,options.stop, options.chainid
        from data import three2one
        loopseq = "" ; loopcounter = options.start
                
        for loopcounter in range(options.start,options.stop+1):
            aafound  = 0
            for k,v in resnums.items():
                if chids[k] != options.chainid :
                    continue
                if int(v) == loopcounter :
                    loopseq = loopseq + three2one[resns[k]]
                    aafound = 1
            if aafound ==  0 :
                print "Amino acid type for residue",loopcounter, "cannot be found from input coordinate file"
                print "Please sequence for loop using --loopseq"
                import sys ;                    sys.exit()
        options.loopres = loopseq
        if options.verbose > 7 :            print "Loop sequence read from file", options.loopres



    if (options.poorOnly == 0 and (options.start != None and options.stop != None and options.chainid != None) and options.caRes == 1):
        loopresids = parseLoopSequenceSC(modelIn,options.loopres,options.start,options.stop,options.chainid,options.mconly,pdboutpath)
        
        badresids = loopresids
              
            
    if options.loopres != None and options.caRes == 0 :
	loopresids,opsaxoff,nullres,missingpts = parseLoopSequence(modelIn,options.loopres,options.start,options.stop,options.chainid,options.mconly,pdboutpath)

	badresids = loopresids 

        print "nullres",nullres
        print "missingpts", missingpts 
        


    ################   COPY MAP AND MTZ FILES  #########################################

    if  options.mapfn == None and options.mtzfn == None :
        if options.poorOnly == 1 :
            print "In order to rebuild poor fitting regions --hklin or --mapin  needs to be specified"
            import sys ; sys.exit()
        print "********** WARNING!! no mtz or map file given, no electron density restraints will be used *************"

    if  options.mapfn != None:
        mapfilepath= "%s"%(options.mapfn)
        if (os.path.isfile(mapfilepath)==False) :
            print "Cannot find file %s " %options.mapfn ; import sys ;            sys.exit()
        shutil.copyfile(mapfilepath, "%s/init.map" % (options.dir_xyzout))

    elif options.mtzfn != None :
        mtzfilepath= "%s"%(options.mtzfn)
        if (os.path.isfile(mtzfilepath)==False) :
            print "Cannot find file %s " %options.mtzfn ; import sys ;sys.exit()
        if  "mtz" not in options.mtzfn : 
            print "Cannot understand i/p HKLIN format,the file should be in either cif (*.cif) or mtz (*.mtz) format " ; import sys ; sys.exit()
            
        shutil.copyfile(mtzfilepath, "%s/init.mtz" % (options.dir_xyzout))
        mtzfilepath = "%s/init.mtz" % (options.dir_xyzout)


        ################   COPY MtZ COLOUMN LABELS  #########################################

        if (options.f1label == None):
            print "Please specify FP label  " , options.f1label ; import sys ; sys.exit()
        if (options.sigf1label == None):
            print "Please specify SIGFP label  " , options.sigf1label ; import sys ; sys.exit()


        if (options.f2label == None or options.phiclabel == None ):
            print "Column labels for  FC and PHIC  not specified, will use input coordinate structure to obtain FC and PHIC"
            sfall(modelIn, options.mtzfn, "phased.mtz",options.resolution,options.usefreer,options.f1label,options.sigf1label,options.freeRlabel)
            options.mtzfn  = "phased.mtz"
            options.f2label  = "FC" ; options.phiclabel = "PHIC"
            


        ### MAKE UNWEIGHTED DIFF MAP : todo sigmaA weighted maps #########################

        from xray import fftnew
        if  type(options.n) != type(1):
            int(options.n)
            print " --n needs to be an integer"

        ccp4map = "%s%dFo-%dFC.map"%(options.mtzfn,options.n,(options.n-1))
        if options.usefreer == 1 : 
            if options.freeRlabel == None :
                print "If you would like to use the freeR set, please specify column name for freeR"
                import sys ; sys.exit()
            fftnew(options.mtzfn, ccp4map, modelIn,options.f1label,options.f2label,options.phiclabel,options.sigf1label,options.n,(options.n)-1,options.freeRlabel)
        else :
            fftnew(options.mtzfn, ccp4map, modelIn , options.f1label,options.f2label,options.phiclabel,options.sigf1label,options.n,(options.n)-1,None)
            
        fcmap = "%sFC.map"%options.mtzfn

        fftnew(options.mtzfn, fcmap, modelIn ,options.f1label,options.f2label,options.phiclabel,options.sigf1label,0,1,options.freeRlabel,1)
        options.mapfn = ccp4map

        
        
    ########### CHECK FOR AVAIALABLE CELL PARAMETERS, SPACE GROUP and RESOLUTION
    if (options.mapfn !=  None or options.mtzfn !=None ):
        from stump import getCRYST , getRESO

        if (options.a == None or options.b == None or options.c == None or options.alpha == None or options.beta == None or options.gamma == None) :
            
            print "Getting cell paramters from coordinate file....."
            options.a,options.b,options.c,options.alpha , options.beta , options.gamma,d1  = getCRYST(pdbfilepath)

            if (options.a == None or options.b == None or options.c == None or options.alpha== None or options.beta==None or options.gamma == None ):
                print "CRYST card cannot be read from coordinate file. Please input cell paramater a, b , c , alpha, beta , gamma = ",options.a , options.b , options.c , options.alpha , options.beta  , options.gamma 
                import sys ; sys.exit()

            print "Found a b c alpha beta gamma  ", options.a , options.b , options.c , options.alpha , options.beta  , options.gamma 
            
        if options.sg == None : 
            print "Getting space group from coordinate file....."
            d1,d2,d3,d4 , d5 , d6, options.sg  = getCRYST(pdbfilepath)
            if options.sg == None : 
                print "Please input space group " , options.sg ; import sys ; sys.exit()
            ss = ""
            for sg1 in options.sg:
                if sg1 in ["\n","\t","\s"]:
                    continue
                else :
                    ss = ss+sg1
            options.sg = ss
            print "Space Group ",options.sg

        if options.sg  in long2shortHM.keys():
            shortsg = long2shortHM[options.sg]
            options.sg = shortsg
        if options.sg not in sgtable.keys():
            print "Check --sg , Not recognised [%s][%d]"%( options.sg, len(options.sg))
            print "Check --sg , Not recognised [%s][%d]"%( options.sg, len(options.sg))
            import sys ; sys.exit()

        if options.resolution == None : 
            print "Getting resolution limit from coordinate file........"
            options.resolution = getRESO(pdbfilepath)
            if (options.resolution == None):
                print "Please input resolution " , options.resolution
                import sys ; sys.exit()
                
            print "Resolution = [ " , options.resolution, " ] "


        ############### IDENTIFY POOR REGIONS  ######################

        
        if options.poorOnly == 1 :
            from restraints import EDrestraint ; EDrestraint.setPenalty(2.) ; EDrestraint.setScatRad(1.) ;
            esmin, esmax, esmean, rcmult, xscoreCutoff = .000, 5., .0, 5, options.poorThreshold
            xscorer = XrayScorer(None, xscoreCutoff)

            mapcoeff = "%dF1-%dF2"%(options.n, (options.n)-1)
            badmcsc,id1 = xscorer.score(options.pdbfile, ccp4map, fcmap, options.f2label, options.phiclabel, mapcoeff, "mcsc",options.verbose)
            badmc,id2 = xscorer.score(options.pdbfile, ccp4map, fcmap, options.f2label, options.phiclabel, mapcoeff, "mc",options.verbose )
            badpept,id3 = xscorer.score(options.pdbfile, ccp4map, fcmap, options.f2label, options.phiclabel, mapcoeff, "pept",options.verbose )
            badsides,id4 = xscorer.score(options.pdbfile, ccp4map, fcmap, options.f2label, options.phiclabel, mapcoeff, "sc",options.verbose)

            bad = list ( set(badmcsc+badmc+badpept) )
            badid = list ( set(id1+id2+id3) )
            prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
            res, resids, resnums, resns, chids, inscodes, pts , atomids , bfac = prot2res.readProtRes2(prot)

            ## IF START AND STOP ARE GIVEN, THEN GET SUBSET OF BAD RESIDUES

            if options.start != None and options.stop != None and options.chainid != None :
                badresids = []
                for bd in badid :
                    if int(resnums[bd]) in range(options.start,options.stop+1) and chids[bd] == options.chainid :
                        badresids.append(resids[bd])

                badscs = []
                for bd in id4 :
                    if int(resnums[bd]) in range(options.start,options.stop+1) and chids[bd] == options.chainid :
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
    
    modelIn =  "modelinit.pdb" ;    rtkmodel = options.pdbout 
    opsaxoff = None ; guidedsampling = None;    scvdwr = options.scReduction ;    popsize = options.popsize
    xrayRestGen = [] ; badresids = None ;
    multiPrepC = prepareChainV3.PrepareChain(options.rotLib)
    if options.mtzfn != None or options.mapfn !=None :
        mapcoeff = "%dF1-%dF2"%(options.n, options.n-1)
        xrayRestGen.append( prepareChain.XrayRestraintsGenerator(ccp4map, "map", options.f2label, options.phiclabel, mapcoeff, esmin, esmax, esmean, [], ) )



    if options.verbose >= 7 : 
        print "Residues to be rebuilt are ", badresids	


    ############### GENERATE MODELS  ######################

    for  i in range(options.nmodels):
        nb = 0 
        print
        print
        print "\nGenerating MODEL %d of %d asked for ..........."%(i+1,options.nmodels)
        if ( ( badresids == None or len(badresids) != 0 or options.poorOnly == 1)): 
            
            ml = Multiloop(modelIn, badresids, options.mconly, options.caRad, options.scRad, scvdwr, guidedsampling, popsize, options.backtrack, 1 , "pre."+str(i)+"."+rtkmodel, xrayRestGen, multiPrepC,options.sconly, options.cacaCutoff,nullres,options.modelN2C,options.natt,missingpts) 
        
            ml.ranker = None 
            if options.mtzfn != None or options.mapfn !=None:
                ml.ranker = XrayRanker(ccp4map, "map", options.f2label, options.phiclabel, mapcoeff, esmin, esmax)
                ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
                ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble = None
                ml.cellsym = [ options.a, options.b, options.c, options.alpha, options.beta, options.gamma, sgtable[options.sg][0] ]


                


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

        if options.opsax == 0 or options.mapfn == None or options.mconly == 1 or opsaxoff == 1:
            os.rename("pre."+str(i)+"."+rtkmodel, str(i)+"."+rtkmodel )

        else :
            print "\n\nStarting OPSAX ........"

            if options.poorOnly == 1 and options.sconly == 0 : 
                if len(unbuilt[0]) == 0  :
                    badresids = list(  set(  findChangedSC(modelIn, "pre."+str(i)+"."+rtkmodel) + badscs )  )
                else :
                    badresids = list(  set( badscs )  )
                    
            elif (options.poorOnly) == 1 and options.sconly == 1 : 
                badresids = badscs 

            elif   options.poorOnly == 0 and len(unbuilt[0]) == 0  :
                badresids = list(  set(  findChangedSC(modelIn, "pre."+str(i)+"."+rtkmodel) + badscs )  )
            else:
                badresids = badresids
                ### todo remove missing sidechains ??
            
            if len(badresids) > 250:
                print "Length of section > 250 residues, can not perform OPSAX"
                os.rename("pre."+str(i)+"."+rtkmodel , str(i)+"."+rtkmodel)
            else :
                if  options.sconly == 1 or nb !=1 :
                    scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 0, 1
                else :
                    scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 1, 1
                
                nb2 = SCplacement("pre."+str(i)+"."+rtkmodel, options.scReduction, str(i)+"."+rtkmodel, "dotfile", useDEE, options.mapfn, options.f1label, options.f2label , options.phiclabel, mapcoeff, esmin, esmax, None, useGivenRot, badresids, scPrepC).run()

                if nb2 != 1:
                    shutil.copy("pre."+str(i)+"."+rtkmodel ,str(i)+"."+rtkmodel)

        adjustBfac2(str(i)+"."+rtkmodel, modelIn,options.nativeBfac, options.mcBfac, options.scBfac)
        print 

        #copyNonprotein(options.pdbfile, str(i)+"."+rtkmodel) ## automatic water addition can go here
        print 
        #copyWater(options.pdbfile, str(i)+"."+rtkmodel) ## automatic water addition can go here
        print
        removeMODEL(str(i)+"."+rtkmodel) ## automatic water addition can go here
        print

        print "------------------------------------ MODEL ASSESSMENT --------------------------------"
        
        if options.mtzfn != None or options.mapfn != None :
            from prefRapperV4WithRefmac import refmacNew
            from protRefine_nick import molProbity_badres, add_cryst_card

            sfcheck(options.mtzfn,options.pdbfile,options.f1label,options.sigf1label,options.freeRlabel ) 
            os.rename("sfcheck.log","sfcheck.init.log")

            sfall(str(i)+"."+rtkmodel, options.mtzfn, "phased."+str(i)+".mtz",options.resolution,options.usefreer,options.f1label,options.sigf1label,options.freeRlabel ) 
            sfcheck("phased."+str(i)+".mtz",str(i)+"."+rtkmodel,options.f1label,options.sigf1label,options.freeRlabel ) 
            os.rename("sfcheck.log","sfcheck."+str(i)+".log")
            print
            print
            print "R-factor and correlation for initial coordinates"
            parse_sfcheck_logfile("sfcheck.init.log")
            print
            print
            print "R-factor and correlation for output (rappertk) coordinates"
            parse_sfcheck_logfile("sfcheck."+str(i)+".log")

            
            add_cryst_card(str(i)+"."+rtkmodel, options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg)
            refmacNew("phased."+str(i)+".mtz", "refmac."+str(i)+".mtz", str(i)+"."+rtkmodel, str(i)+".refmac.pdb" )
            
            sfcheck("refmac."+str(i)+".mtz",str(i)+".refmac.pdb",options.f1label,options.sigf1label,options.freeRlabel )
            os.rename("sfcheck.log","sfcheck.refmac."+str(i)+".log")
            print
            print "R-factor and correlation for refined (refmac) coordinates"
            parse_sfcheck_logfile("sfcheck.refmac."+str(i)+".log")

            
            refmacNew(options.mtzfn, "refmac.init.mtz", options.pdbfile, "refmac.init.pdb" )
            sfcheck("refmac.init.mtz","refmac.init.pdb",options.f1label,options.sigf1label,options.freeRlabel )
            os.rename("sfcheck.log","sfcheck.refmac.init.log")
            print
            print "R-factor and correlation for refined (refmac) coordinates"
            parse_sfcheck_logfile("sfcheck.refmac.init.log")

            

            #if (os.path.isfile(pdbout)==False) :
            #    print "Refmac cycle %d failed, no %s generated"%(rcycle,pdbout)
            #    import sys ; sys.exit()
            #if (os.path.isfile(hklout)==False) :
            #    print "Refmac cycle %d failed, no %s generated"%(rcycle,hklout)
            #    import sys ; sys.exit()

            
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
                print "The C-alpha atom for residue [%s] is missing in the input coordinates"%resids[k]
                import sys ; sys.exit()
    return badresids

    

def removeChainId(pdbfile):

    from pdbr import line2bfac, isPdbAtomLine, line2bfac, line2atomid, changeBfactor, line2crdstr
    from pdbinfo import segi
    newlines = []

    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        import sys ; sys.exit()

    ip = open(pdbfile, 'r')
    fileList = ip.readlines()
    for l in fileList:
        
        if isPdbAtomLine(l):
            l  =  l[0:21] + ' ' + l[22:]
            newlines.append(l)
        else :
            newlines.append(l)

    op = open(pdbfile, 'w')
    for l in newlines : op.write(l)
    op.close()








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
    print "Copy Bfactors from", refpdb, "to", pdbfile, "if coordinates are same, else set to 30"
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

if __name__ == "__main__" :
    main()
    import sys ; sys.exit(0)
    from scplacement import SCplacement
    import prepareChain
    scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 1, 1
    badresids = ["VAL   85 ", "ASP   86 ", "TYR   68 ", "TYR   90 ",],
    SCplacement("premodel2.pdb", 0.5, "mmm.pdb", "dotfile", useDEE, "phased1.mtz2fofc.map", "FP", "FC", "PHIC", "2F1-F2", 0, 5, None, useGivenRot,
        badresids, scPrepC).run()
    import sys ; sys.exit(0)
    replaceWaters("model1.pdb", "rtk0.map")
