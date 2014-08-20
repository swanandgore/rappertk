import os, shutil, re , sys
import geometry
from xray import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS, fft, omitmap, mapman
from procrun import proc_run_exitOnError as execCmd
from xcheck import XrayScorer, XrayRanker
from data import sgtable
from evalCAtrace import comparePhiPsiOmegaChi
from pdbr import protein, isAAres  , line2atomname
import prot2res
from scplacement import SCplacement
from loopbuild import Multiloop
import prepareChain
from multProtref import joinPDBs, splitPDBs, restoreChainids

def removeRElines(pdbfile, regexps) :
    import re
    newlines, compRegexps = [], []
    for rege in regexps : compRegexps.append( re.compile(rege) )
    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        sys.exit()
    for l in open(pdbfile,'r').readlines() :
        include = 1
        for rege in compRegexps :
            if rege.search(l) != None : include = None ; break
        if include != None : newlines.append(l)
    op = open(pdbfile, 'w')
    for l in newlines : op.write(l)
    op.close()

def userpause() :
    print "press key and enter to proceed" ;
    import sys; sys.stdin.readline()
    print "cheers!"

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
            for xi in range(3) : cen[xi] += pts[ai][xi]
            cnt += 1
        for xi in range(3) : cen[xi] /= cnt
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
        sys.exit()
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






def fixCNSop(pdbfile) :
    from pdbr import isPdbAAline, isPdbAtomLine, line2atomname, line2resn, changeResn, changeAtomname
    lines = []
    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        sys.exit()
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

## assert that there is only 1 model, and remove model, endmdl line
def removeMODEL(pdbfile) :
    newlines, nm = [], 0
    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        sys.exit()

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
        sys.exit()

        

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
    print input

    execCmd("lx_moleman", input)
    os.rename(pdbfilename+".moleman", pdbfilename)

## verify that last line in topdb in END, change that to TER
## then copy all non-protein lines
def copyNonprotein(frompdb, topdb, waterAlso=1) :
    from pdbr import isPdbAAline, isPdbAtomLine, changeSegid
    if (os.path.isfile(topdb)==False) :
        print "Cannot find file %s "%topdb
        print "No file in directory ", os.getcwd()
        sys.exit()
    
    lines = open(topdb, 'r').readlines()
    #assert lines[ len(lines)-1 ][0:3] == "END"
    if lines[ len(lines)-1 ] in [ "END\n", "END " ] : lines[ len(lines)-1 ] = "TER\n"
    else : lines.append("TER\n")
    print "copying nonprotein atoms from", frompdb, "to", topdb


    if (os.path.isfile(frompdb)==False) :
        print "Cannot find file %s "%frompdb
        print "No file in directory ", os.getcwd()
        sys.exit()

    for l in open(frompdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if isPdbAAline(l) : continue
        if not waterAlso and line2resn(l) == "HOH" : continue
        lines.append( changeSegid(l,"    ") )
    op = open(topdb, 'w')
    for l in lines : op.write(l)
    print >> op, "TER\nEND"
    op.close()

def adjustBfac(pdbfile, refpdb) :
    print "Copy Bfactors from", refpdb, "to", pdbfile, "if coordinates are same, else set to 30"
    from pdbr import line2bfac, isPdbAtomLine, line2bfac, line2atomid, changeBfactor, line2crdstr
    id2bo = {} # read in crd, bfac, for each atom-id

    
    if (os.path.isfile(refpdb)==False) :
        print "Cannot find file %s "%refpdb
        print "No file in directory ", os.getcwd()
        sys.exit()


    for l in open(refpdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if id2bo.has_key(line2atomid(l)) : continue
        id2bo[line2atomid(l)] = [line2bfac(l), line2crdstr(l)]
    newlines = []


    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        sys.exit()



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
        sys.exit()

    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        sys.exit()



    for l in open(refpdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if id2bo.has_key(line2atomid(l)) : continue
        id2bo[line2atomid(l)] = [line2bfac(l), line2crdstr(l)]
    newlines = []



    if nativeBfac == 1:
        for l in open(pdbfile, 'r').readlines() :
            if not isPdbAtomLine(l) : newlines.append(l) ; continue
            bfac, crdstr = id2bo[ line2atomid(l) ]
            l = changeBfactor(l, bfac)
            newlines.append(l)

    else :
        for l in open(pdbfile, 'r').readlines() :
            if not isPdbAtomLine(l) : newlines.append(l) ; continue
            bfac, crdstr = id2bo[ line2atomid(l) ]
            if crdstr == line2crdstr(l) :
                l = changeBfactor(l, bfac)
            else :
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
def pdbCopy(pdbfrom, pdbto) :
    from pdbr import isPdbAtomLine, line2atomid, isPdbAAline, line2atomname, changeBfactorOccupancy
    op = open(pdbto, 'w')
    idsSofar = []
    if (os.path.isfile(pdbfrom)==False) :
        print "Cannot find file %s "%pdbfrom
        print "No file in directory ", os.getcwd()
        sys.exit()


    for l in open(pdbfrom, 'r').readlines() :
        if not isPdbAtomLine(l) : op.write(l) ; continue
        if line2atomid(l) in idsSofar : continue
        idsSofar.append(line2atomid(l))
        bfactor = 30
        if isPdbAAline(l) and line2atomname(l) in [' N  ',' CA ',' C  ',' O  '] : bfactor = 20
        l = changeBfactorOccupancy(l, bfactor, 1)
        op.write(l)
    op.close()

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
        sys.exit()

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
    print "renaming", inpdb, pdbout
    os.rename(inpdb, pdbout)

def cnsRefinement(mtzin,pdbin, mtzout,pdbout, a,b,c,alpha,beta,gamma,sg,reso, cnsArgs,cycle, extraTOPfile=None, extraPARfile=None) :
    mtz2hkl(mtzin, "cns.hkl")
    cns_generate(pdbin, "generate.mtf", "generate.pdb", extraTOPfile, extraPARfile, "generate.log")
    removeZeroLines("generate.pdb") ## ???
    wa = -1 ; harmCA = None
    if cnsArgs[cycle].has_key("harmCA") and cnsArgs[cycle]["harmCA"] != None : harmCA = 1
    cns_anneal(a, b, c, alpha, beta, gamma, sgCCP4toCNS[sg], reso,
        "cns.hkl", "generate.mtf", "generate.pdb", extraPARfile, "anneal%d.log"%cycle, wa, cnsArgs[cycle]["num_cycles"], cnsArgs[cycle]["temperature"], harmCA)
    removeZeroLines("anneal.pdb") ##  ???
    fixCNSop("anneal.pdb")
    os.rename("anneal.pdb", pdbout)
    sfall(pdbout, "rfree.mtz", mtzout, reso)
    mapman("anneal_2fofc.map", mtzout+"2fofc.map")
    mapman("anneal_fc.map", mtzout+"fc.map")
    #moleman(pdbout)

def makeRMSDcomments(refpdb, cnsout, loopres) :
    from pdbr import protein
    refprot = protein(refpdb, read_hydrogens=0, read_waters=0, read_hets=0)
    testprot = protein(cnsout, read_hydrogens=0, read_waters=0, read_hets=0)
    chi1accu, chi1diff, chi12accu, chi12diff, carmsd, mcrmsd, aarmsd, top10 = comparePhiPsiOmegaChi(refprot, testprot)
    remark = "REMARK chain-comparison cns0 ca %6.3f mc %6.3f aa %6.3f chi1 %6.3f %6.3f chi12 %6.3f %6.3f\n" % (carmsd, mcrmsd, aarmsd, chi1accu, chi1diff, chi12accu, chi12diff)
    if loopres :
        chi1accu, chi1diff, chi12accu, chi12diff, carmsd, mcrmsd, aarmsd, top10 = comparePhiPsiOmegaChi(refprot, testprot, None, loopres)
        remark += "REMARK  loop-comparison cns0 ca %6.3f mc %6.3f aa %6.3f chi1 %6.3f %6.3f chi12 %6.3f %6.3f\n" % (carmsd, mcrmsd, aarmsd, chi1accu, chi1diff, chi12accu, chi12diff)
    remark += "REMARK top10 culprits "
    for t10 in top10 : remark += "[%s] " % t10
    remark += "\n"
    if (os.path.isfile(cnsout)==False) :
        print "Cannot find file %s "%cnsout
        print "No file in directory ", os.getcwd()
        sys.exit()

    lines = [remark,] + open(cnsout, 'r').readlines()
    ofp = open(cnsout,'w')
    for l in lines : ofp.write(l)
    ofp.close()

def main() :

    import optparse ; parser = optparse.OptionParser()

    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='to create all the files during refinement. it shdnt be already present.')
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='PDB input structure to be refined ')
    parser.add_option("--outpdb", action='store', type='string', dest='pdbout', help='PDB output file',default="modelout.pdb")
    parser.add_option("--mapfn1", action='store', type='string', dest='mapfn1', help='mapfn1, example: 2Fo-Fc, Fo',default=None)
    parser.add_option("--mapfn2", action='store', type='string', dest='mapfn2', help='mapfn2, example: Fc map ',default=None)

    parser.add_option("--mtz", action='store', type='string', dest='mtzfn', help='phased mtz file',default=None)
    parser.add_option("--a", action='store', type='float', dest='a', help='cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='cell spacegroup, in CCP4 notation')
    
        
    parser.add_option("--resolution", action='store', type='float', dest='resolution', help='resolution of the data')

    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid', default=2)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=1)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)

    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='seed for randomizing', default=None)
    parser.add_option("--loopres", action='store', type='string', dest='loopres', help='filename containing resids for starting perturbation', default=None)
    parser.add_option("--framework", action='store', type='int', dest='framework', help='to be used in conjunction with loopres. it puts a framwork-ca-threshold/framework-sc-threshold pos restr on non-loopres and perturbs them too', default=None)
    parser.add_option("--framework-ca-threshold", action='store', type='float', dest='fcarad', help='to be used in conjunction with loopres.', default=1.)
    parser.add_option("--framework-sc-threshold", action='store', type='float', dest='fscrad', help='to be used in conjunction with loopres', default=3.)

    ### Added params to be RAPPER-like


    parser.add_option("--rebuild-poor-regions-only", action='store', type='int', dest='poorOnly', help='To rebuild regions of poor fit to an electron density map. Residues to be rebuilt are identified using a real space scoring function, the cut off for which is set using --poor-region-threshold.Default = False', default=0)
    parser.add_option("--poor-fit-threshold", action='store', type='float', dest='poorThreshold', help='Threshold to identify ill fitting regions', default=0.9)


    parser.add_option("--mconly", action='store', type='int', dest='mconly', help='Build mainchain only', default=None)
    parser.add_option("--sconly", action='store', type='int', dest='sconly', help='Build side chains only', default=None)
        


    parser.add_option("--default-mainchain-b-factor", action='store', type='float', dest='mcBfac', help='The b-factor assigned to the newly built main chain region', default=20.)
    parser.add_option("--default-sidechain-b-factor", action='store', type='float', dest='scBfac', help='The b-factor assigned to the newly built side chain region', default=30.)
    parser.add_option("--models-get-native-bfactors", action='store', type='int', dest='nativeBfac', help='Take the b-factors from the section that is being rebuilt.', default=0)

    parser.add_option("--dontusesc", action='store', type='int', dest='dontusesc', help='Dont reassign sc with scep, set to 1 ', default=0)

    parser.add_option("--num-models", action='store', type='int', dest='nmodels', help='number of models wanted ', default=1)
    parser.add_option("--f1label", action='store', type='string', dest='f1label', help='f1label', default="FP")
    parser.add_option("--f2label", action='store', type='string', dest='f2label', help='f2label', default="FC")
    parser.add_option("--philabel", action='store', type='string', dest='philabel', help='philabel', default="PHIC")
    parser.add_option("--maptype", action='store', type='string', dest='maptype', help='maptype can only be 2F1-F2 or F1', default="2F1-F2")

    parser.add_option("--cacaCutoff", action='store', type='float', dest='cacaCutoff', help='min dist between adjacent CA to detect a chain-break', default=5.)
    parser.add_option("--make-ed-optional", action='store', type='int', dest='edOpt', help='    If = 1, then the 0.0 (positive density) mainchain restraint will be made optional. If false, then the main chain will be unconditionally forced to lie in positive density. This is primarily useful when tracing through a structure with regions in very poor (non-existent) density', default=None)    
    #parser.add_option("--chids", action='store', type='string', dest='chids', help='chains to be traced', default=' ')
    

    #parser.add_option("--start", action='store', type='int', dest='start', help='Residue number to start building from ', default=None)
   #parser.add_option("--stop", action='store', type='int', dest='stop', help='Residue number to stop building ', default=None)
    #parser.add_option("--chain-id", action='store', type='string', dest='chainid', help='Chain ID of section to be built. If all chains to be built use *', default=' ')


    #parser.add_option("--ed-mainchain-min-sigma", action='store', type='float', dest='edmcSigma', help='    Only main chain atoms in a position with greater standard deviation than this are considered to satisfy the electron density map restraint.', default=0.5)

    #parser.add_option("--ed-sidechain-min-sigma", action='store', type='float', dest='edscSigma', help='    Only side chain atoms in a position with greater standard deviation than this are considered to satisfy the electron density map restraint.', default=0.5)

    



    (options, args) = parser.parse_args()
    
    ###########################################
    

    print "========================== Now  printing all KEYWORD values ==========================="
    print "--scratchdir = ",options.scratchdir
    print "--pdb = ",options.pdbfile
    print "--outpdb = ",options.pdbout
    print "--mapfn1 = ",options.mapfn1
    print "--mapfn2 = ",options.mapfn2
    print "--mtz = ",options.mtzfn

    print "--a = ",options.a
    print "--b = ",options.b
    print "--c = ",options.c
    print "--alpha = ",options.alpha
    print "--beta = ",options.beta
    print "--gamma = ",options.gamma
    print "--sg = ",options.sg
    print "--resolution = ",options.resolution
    print "--ca-restraint-radius = ",options.caRad
    print "--sc-centroid-restraint-radius = ",options.scRad
    print "--sidechain-vdw-reduction = ",options.scReduction
    print "--population-size = ",options.popsize
    print "--verbose = ",options.verbose
    print "--backtrack = ",options.backtrack
    print "--randomize = ",options.randomize
    print "--loopres = ",options.loopres
    print "--framework = ",options.framework
    print "--framework-ca-threshold = ",options.fcarad
    print "--framework-sc-threshold = ",options.fscrad
    print "--rebuild-poor-regions-only = ",options.poorOnly
    print "--poor-fit-threshold = ",options.poorThreshold
    print "--default-mainchain-b-factor = ",options.mcBfac
    print "--default-sidechain-b-factor = ",options.scBfac
    print "--models-get-native-bfactors = ",options.nativeBfac
    print "--dontusesc = ",options.dontusesc
    print "--num-models = ",options.nmodels
    print "--f1label = ",options.f1label
    print "--f2label = ",options.f2label
    print "--philabel = ",options.philabel
    print "--maptype = ",options.maptype
    print "--make-ed-optional = ",options.edOpt
    print "--mconly = ",options.mconly
    print "--sconly = ",options.sconly
#    print "--chids = ",options.chids
    




    
    
    print "========================== End KEYWORD values ==========================="


    ########################################
    assert options.framework in [None, 1]
    import misc
    misc.setVerbosity(options.verbose)
    randomize(options.randomize)

    if options.mtzfn !=None or options.mapfn1 !=None : 
        if (options.a == None or options.b == None or options.c == None or options.alpha== None or options.beta==None or options.gamma == None or options.sg == None or options.resolution == None ):
            print "Please check cell paramater a, b , c , alpha, beta , gamma = ",options.a , options.b , options.c , options.alpha , options.beta  , options.gamma 
            print "Please check space group " , options.sg
            sys.exit()

        if options.resolution==None :
            print "Please check value of resolution = ", options.resolution
            sys.exit()

    if options.loopres != None :
        loopresids = parseLoopres(options.pdbfile,options.loopres)
        options.loopres = loopresids 
        



    if not os.path.isdir(options.scratchdir) : os.mkdir(options.scratchdir)
    

    if (os.path.isfile(options.pdbfile)==False) :
        print "Cannot find file %s "%options.pdbfile
        print "No file in directory ", os.getcwd()
        sys.exit()

    shutil.copyfile(options.pdbfile, "%s/%s" % (options.scratchdir,options.pdbfile))




    if options.mapfn1 != None:
        if (os.path.isfile(options.mapfn1)==False) :
            print "Cannot find file %s "%options.mapfn1
            print "No file in directory ", os.getcwd()
            sys.exit()

        shutil.copyfile(options.mapfn1, "%s/%s" % (options.scratchdir,options.mapfn1))


    if options.mapfn2 != None:
        if (os.path.isfile(options.mapfn2)==False) :
            print "Cannot find file %s "%options.mapfn2
            print "No file in directory ", os.getcwd()
            sys.exit()

        shutil.copyfile(options.mapfn2, "%s/%s" % (options.scratchdir,options.mapfn2))

    if options.mtzfn != None:
        if (os.path.isfile(options.mtzfn)==False) :
            print "Cannot find file %s "%options.mtzfn
            print "No file in directory ", os.getcwd()
            sys.exit()

        shutil.copyfile(options.mtzfn, "%s/%s" % (options.scratchdir,options.mtzfn))
        

    os.chdir(options.scratchdir)
    debug = 1 
    if debug == 1 :
        uniqueify(options.mtz, "phased.mtz")
        sfall(options.pdbfile, "rfree.mtz", "phased.mtz")
        sys.exit()
        

    fp = open("parameters.txt", 'w')
    l =  "--scratchdir = %s "%str(options.scratchdir) ;      print >> fp, l 
    l =  "--pdb = %s "%str(options.pdbfile) ;      print >> fp, l 
    l =  "--outpdb = %s "%str(options.pdbout) ;      print >> fp, l 
    l =  "--mapfn1 = %s "%str(options.mapfn1) ;      print >> fp, l 
    l =  "--mapfn2 = %s "%str(options.mapfn2) ;      print >> fp, l 
    l =  "--mtz = %s "%str(options.mtzfn) ;      print >> fp, l 
    l =  "--a = %s "%str(options.a) ;      print >> fp, l 
    l =  "--b = %s "%str(options.b) ;      print >> fp, l 
    l =  "--c = %s "%str(options.c) ;      print >> fp, l 
    l =  "--alpha = %s "%str(options.alpha) ;      print >> fp, l 
    l =  "--beta = %s "%str(options.beta) ;      print >> fp, l 
    l =  "--gamma = %s "%str(options.gamma) ;      print >> fp, l 
    l =  "--sg = %s "%str(options.sg) ;      print >> fp, l 
    l =  "--resolution = %s "%str(options.resolution) ;      print >> fp, l 
    l =  "--ca-restraint-radius = %s "%str(options.caRad) ;      print >> fp, l 
    l =  "--sc-centroid-restraint-radius = %s "%str(options.scRad) ;      print >> fp, l 
    l =  "--sidechain-vdw-reduction = %s "%str(options.scReduction) ;      print >> fp, l 
    l =  "--population-size = %s "%str(options.popsize) ;      print >> fp, l 
    l =  "--verbose = %s "%str(options.verbose) ;      print >> fp, l 
    l =  "--backtrack = %s "%str(options.backtrack) ;      print >> fp, l 
    l =  "--randomize = %s "%str(options.randomize) ;      print >> fp, l 
    l =  "--loopres = %s "%str(options.loopres) ;      print >> fp, l 
    l =  "--framework = %s "%str(options.framework) ;      print >> fp, l 
    l =  "--framework-ca-threshold = %s "%str(options.fcarad) ;      print >> fp, l 
    l =  "--framework-sc-threshold = %s "%str(options.fscrad) ;      print >> fp, l 
    l =  "--rebuild-poor-regions-only = %s "%str(options.poorOnly) ;      print >> fp, l 
    l =  "--poor-fit-threshold = %s "%str(options.poorThreshold) ;      print >> fp, l 
    l =  "--mconly = %s "%str(options.mconly) ;      print >> fp, l 
    l =  "--sconly = %s "%str(options.sconly) ;      print >> fp, l 
    l =  "--default-mainchain-b-factor = %s "%str(options.mcBfac) ;      print >> fp, l 
    l =  "--default-sidechain-b-factor = %s "%str(options.scBfac) ;      print >> fp, l 
    l =  "--models-get-native-bfactors = %s "%str(options.nativeBfac) ;      print >> fp, l 
    l =  "--dontusesc = %s "%str(options.dontusesc) ;      print >> fp, l 
    l =  "--num-models = %s "%str(options.nmodels) ;      print >> fp, l 
    l =  "--f1label = %s "%str(options.f1label) ;      print >> fp, l 
    l =  "--f2label = %s "%str(options.f2label) ;      print >> fp, l 
    l =  "--philabel = %s "%str(options.philabel) ;      print >> fp, l 
    l =  "--maptype = %s "%str(options.maptype) ;      print >> fp, l 
    l =  "--cacaCutoff = %s "%str(options.cacaCutoff) ;      print >> fp, l 
    l =  "--make-ed-optional = %s "%str(options.edOpt) ;      print >> fp, l
#    l =  "--chids = %s "%str(options.chids) ;      print >> fp, l 
    fp.close()
    
    guidedsampling = None
    from restraints import EDrestraint ; EDrestraint.setPenalty(2.) ; EDrestraint.setScatRad(1.) ;
    esmin, esmax, esmean, rcmult, xscoreCutoff = .000, 5., .0, 5, options.poorThreshold
    multiPrepC = prepareChain.PrepareChain("SCL1.0")
    xscorer = XrayScorer(None, xscoreCutoff)

    
    modelIn =  options.pdbfile
    rtkmodel = options.pdbout # rappertk model to be generated in this cycle


    if options.mapfn1 == None and options.mtzfn == None:
        print "Warning no mtz or map input, no electron density restraints will be used"
    
    if options.poorOnly == 1 :
        if  options.mtzfn == None and options.mapfn1 == None and options.mapfn2 == None:
            print "No mapfn1,mapfn2 or mtzfn  set"
            sys.exit()
        if options.mapfn2 == None and options.mapfn1 !=None and options.mtzfn == None:
            print "Both mapfn1 and mapfn2 need to be  set"
            sys.exit()            

        if options.mapfn1 == None and options.mapfn2 !=None and options.mtzfn == None :
            print "Both mapfn1 and mapfn2 need to be  set"
            sys.exit()            


        if options.mtzfn != None :
            badmcsc = xscorer.score(modelIn, options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, "mcsc")
            badmc = xscorer.score(modelIn, options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, "mc")
            badpept = xscorer.score(modelIn, options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype,"pept")
            badscs = xscorer.score(modelIn, options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype,"sc")
            
        elif options.mapfn1 !=None and options.mapfn2 != None  :
            
            badmcsc = xscorer.score(modelIn, options.mapfn1 , options.mapfn2, options.f2label, options.philabel, options.maptype, "mcsc")
            badmc = xscorer.score(modelIn, options.mapfn1, options.mapfn2, options.f2label, options.philabel, options.maptype, "mc")
            badpept = xscorer.score(modelIn, options.mapfn1, options.mapfn2, options.f2label, options.philabel, options.maptype,"pept")
            badscs = xscorer.score(modelIn, options.mapfn1, options.mapfn2, options.f2label, options.philabel, options.maptype,"sc")
        
            
        else :
                print "Need  phased mtz file or map file , see help pages"
                sys.exit()

        
        badresids = list ( set(badmcsc+badmc+badpept) )
        if len(badresids) < 1 and len(badscs) < 0 :
            print "No poor region identified, try higher cutoff if required" ; sys.exit(0)




    
    xrayRestGen = []
    if options.loopres and options.framework :
        xrayRestGen.append( prepareChain.CArestraintGenerator(options.loopres, options.fcarad, options.fscrad) )
        
    if options.mtzfn != None : 
        xrayRestGen.append( prepareChain.XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, esmin, esmax, esmean, ["SCLbuilder","ChiBuilder","CBbuilder"], options.edOpt ) )
    elif options.mapfn1 != None  :
        xrayRestGen.append( prepareChain.XrayRestraintsGenerator(options.mapfn1, "map", options.f2label, options.philabel, options.maptype, esmin, esmax, esmean, ["SCLbuilder","ChiBuilder","CBbuilder"], options.edOpt ) )

    if options.loopres and not options.framework :
        badresids = options.loopres
    
    elif options.loopres and options.framework :
        badresids = None 
        
    elif options.poorOnly == 1 and not options.framework and not options.loopres :
        badresids = badresids

    else:
        badresids = None 


    scvdwr = options.scReduction ;
    popsize = options.popsize


    for  i in range(options.nmodels):
        ml = Multiloop(modelIn, badresids, options.mconly, options.caRad, options.scRad, scvdwr, guidedsampling, popsize,
                       options.backtrack, 1 , "pre"+rtkmodel+"."+str(i), xrayRestGen, multiPrepC)
        
        ml.ranker = None 
        ml.cacaCutoff =  options.cacaCutoff
        if options.mtzfn != None:
            ml.ranker = XrayRanker(options.mtzfn, options.f1label,  options.f2label, options.philabel, options.maptype, esmin, esmax)
            ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
            ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble = None
            
            if options.sg not in sgtable.keys():
                print "Check --sg , Not recognised", options.sg
                sys.exit()
            ml.cellsym = [ options.a, options.b, options.c, options.alpha, options.beta, options.gamma, sgtable[options.sg][0] ]


        elif options.mapfn1 != None :

            ml.ranker = XrayRanker(options.mapfn1, options.f1label,  options.f2label, options.philabel, options.maptype, esmin, esmax)
            ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
            ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble = None
            if options.sg not in sgtable.keys():
                print "Check --sg , Not recognised", options.sg
                sys.exit()
            ml.cellsym = [ options.a, options.b, options.c, options.alpha, options.beta, options.gamma, sgtable[options.sg][0] ]

        else :
            ml.ranker = None 

        if options.sconly == None : 
            nb = ml.run()

        else :
            shutil.copy(modelIn, "pre"+rtkmodel+"."+str(i))

        if options.dontusesc==1 or (options.mapfn1==None and options.mtzfn==None) or options.mconly == 1 :
            os.rename("pre"+rtkmodel+"."+str(i), rtkmodel+"."+str(i))
        else :
            if (options.poorOnly) == 1 : 
                badresids = list(  set(  findChangedSC(modelIn, "pre"+rtkmodel+"."+str(i)) + badscs )  )
            elif options.sconly == 1 :
                if options.loopres == None :
                    badresids = list(  set(  findAllSC(modelIn)  )  )
                else :
                    print badresids
            else :
                badresids = list(  set(  findChangedSC(modelIn, "pre"+rtkmodel+"."+str(i)) )  )


            scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 1, 1
            if options.mtzfn != None:
                SCplacement("pre"+rtkmodel+"."+str(i), options.scReduction, rtkmodel+"."+str(i) , "dotfile", useDEE, options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, esmin, esmax, None, useGivenRot, badresids, scPrepC).run()
            elif options.mapfn1 != None:
                SCplacement("pre"+rtkmodel+"."+str(i), options.scReduction, rtkmodel+"."+str(i) , "dotfile", useDEE, options.mapfn1, options.f1label, options.f2label, options.philabel, options.maptype, esmin, esmax, None, useGivenRot, badresids, scPrepC).run() 
            adjustBfac2(rtkmodel+"."+str(i), modelIn,options.nativeBfac, options.mcBfac, options.scBfac)
            copyNonprotein(modelIn, rtkmodel+"."+str(i)) ## automatic water addition can go here
        #replaceWaters(rtkmodel, rtkmap, 3)
        ## remove existing waters and add new waters using findwaters, use chain W


def removeChainId(pdbfile):

    from pdbr import line2bfac, isPdbAtomLine, line2bfac, line2atomid, changeBfactor, line2crdstr
    from pdbinfo import segi
    newlines = []

    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        sys.exit()

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
        sys.exit()
    
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
        sys.exit()
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
        sys.exit()




    for l in open(refpdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if id2bo.has_key(line2atomid(l)) : continue
        id2bo[line2atomid(l)] = [line2bfac(l), line2occu(l), line2crdstr(l)]
    newlines = []

    if (os.path.isfile(pdbfile)==False) :
        print "Cannot find file %s "%pdbfile
        print "No file in directory ", os.getcwd()
        sys.exit()

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
