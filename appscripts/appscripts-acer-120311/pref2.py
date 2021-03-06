import os, shutil, re
import geometry
from xray2 import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS, fft, omitmap, mapman
from procrun import proc_run_exitOnError as execCmd
from xcheck2 import XrayScorer, XrayRanker
from data import sgtable
from evalCAtrace import comparePhiPsiOmegaChi

def removeRElines(pdbfile, regexps) :
    import re
    newlines, compRegexps = [], []
    for rege in regexps : compRegexps.append( re.compile(rege) )
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

## change ILE:CD to ILE:CD1 if necessary
## OT1 to O
def fixCNSop(pdbfile) :
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

## assert that there is only 1 model, and remove model, endmdl line
def removeMODEL(pdbfile) :
    newlines, nm = [], 0
    for l in open(pdbfile,'r').readlines() :
        if l[0:5] == 'MODEL' : nm += 1
        elif l[0:6] == 'ENDMDL' : continue
        else : newlines.append(l)
    #assert nm == 1
    op = open(pdbfile, 'w')
    for l in newlines : op.write(l)
    op.close()

def removeZeroLines(filename) :
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
    lines = open(topdb, 'r').readlines()
    #assert lines[ len(lines)-1 ][0:3] == "END"
    if lines[ len(lines)-1 ] in [ "END\n", "END " ] : lines[ len(lines)-1 ] = "TER\n"
    else : lines.append("TER\n")
    print "copying nonprotein atoms from", frompdb, "to", topdb
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
    for l in open(refpdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if id2bo.has_key(line2atomid(l)) : continue
        id2bo[line2atomid(l)] = [line2bfac(l), line2crdstr(l)]
    newlines = []
    for l in open(pdbfile, 'r').readlines() :
        if not isPdbAtomLine(l) : newlines.append(l) ; continue
        bfac, crdstr = id2bo[ line2atomid(l) ]
        if crdstr == line2crdstr(l) : l = changeBfactor(l, bfac)
        else : l = changeBfactor(l, 30.0)
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
    for l in open(filename, 'r').readlines() :
        if isPdbAAline(l) :
            if line2atomname(l) in [" N  "," CA "," C  "," O  "] : l = changeBfactorOccupancy(l, bfacs[0], line2occu(l))
            else : l = changeBfactorOccupancy(l, bfacs[1], line2occu(l))
        lines.append(l)
    fp = open(filename, 'w')
    for l in lines : fp.write(l)
    fp.close()

def refmacRefinement(mtzin,pdbin, mtzout,pdbout, a,b,c,alpha,beta,gamma,sg,reso, ccp4Args,cycle) :
    from xray import refmac1
    if os.path.isfile(pdbout) : return
    inpdb, bufpdb = pdbin, "buf.pdb"
    for ccpa in ccp4Args[cycle] :
        if "assignBfac" in ccpa.keys() : changeBfacs(inpdb, ccpa[assignBfac])
        refmac1(mtzin, mtzout, inpdb, bufpdb, reso, ccpa["reftype"], ccpa["wa"], ccpa["breftype"], ccpa["ncyc"])
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





def cnsBfacRefinement(mtzin,pdbin, mtzout,pdbout, a,b,c,alpha,beta,gamma,sg,reso, cnsArgs,cycle, extraTOPfile=None, extraPARfile=None) :
    mtz2hkl(mtzin, "cns.hkl")
    cns_generate(pdbin, "generate.mtf", "generate.pdb", extraTOPfile, extraPARfile, "generate.log")
    removeZeroLines("generate.pdb") ## ???
    wa = -1 ; harmCA = None
    if cnsArgs[cycle].has_key("harmCA") and cnsArgs[cycle]["harmCA"] != None : harmCA = 1
    cns_bindividual(a, b, c, alpha, beta, gamma, sgCCP4toCNS[sg], reso,
        "cns.hkl", "generate.mtf", "generate.pdb", extraPARfile, "bindividual%d.log"%cycle, wa, cnsArgs[cycle]["num_cycles"], cnsArgs[cycle]["temperature"], harmCA)
    removeZeroLines("bindividual.pdb") ##  ???
    fixCNSop("bindividual.pdb")
    os.rename("bindividual.pdb", pdbout)

    #sfall(pdbout, "rfree.mtz", mtzout, reso)
    #mapman("anneal_2fofc.map", mtzout+"2fofc.map")
    #mapman("anneal_fc.map", mtzout+"fc.map")
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
    lines = [remark,] + open(cnsout, 'r').readlines()
    ofp = open(cnsout,'w')
    for l in lines : ofp.write(l)
    ofp.close()

def main() :

    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='to create all the files during refinement. it shdnt be already present.')
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='starting pdb containing a model of pdb-ligand complex')

    parser.add_option("--sf", action='store', type='string', dest='sf', help='structure factors file')
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
    parser.add_option("--noRTK", action='store', type='int', dest='noRTK', help='dont rebuild bad-fits with rtk', default=0)
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='seed for randomizing', default=None)
    parser.add_option("--loopres", action='store', type='string', dest='loopres', help='filename containing resids for starting perturbation', default=None)
    parser.add_option("--framework", action='store', type='int', dest='framework', help='to be used in conjunction with loopres. it puts a 2/3 ca/sc pos restr on non-loopres and perturbs them too', default=None)


    (options, args) = parser.parse_args()

    assert options.framework in [None, 1]

    import misc
    misc.setVerbosity(options.verbose)

    randomize(options.randomize)

    if options.loopres != None : ## expecting resids in [RSNcRESIa] format, one per line
        loopresids = []
        for l in open(options.loopres, 'r').readlines() :
            if l[0] == "#" : continue
            lres = l[1:len(l)-2]
            loopresids.append(lres)
        options.loopres = loopresids # [ "VAL   47 ", "GLN   48 ", "GLY   49 ", "GLU   50 ", "GLU   51 ", "SER   52 ", "ASN   53 ", "ASP   54 ", "LYS   55 ", ]

    if not os.path.isdir(options.scratchdir) : os.mkdir(options.scratchdir)
    shutil.copyfile(options.pdbfile, "%s/model0.pdb" % options.scratchdir)
    shutil.copyfile(options.sf, "%s/strfactors.mtz" % options.scratchdir)
    os.chdir(options.scratchdir)
    from stump import getCRYST , getRESO
    if options.a == None : 
        options.a,options.b,options.c,options.alpha , options.beta , options.gamma,options.sg  = getCRYST(options.pdbfile)
    if options.resolution  == None:
        options.resolution = getRESO(options.pdbfile)
    
    cif2mtz("strfactors.mtz", "base.mtz", options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg)
    uniqueify("base.mtz", "rfree.mtz")

    from scplacement import SCplacement
    from loopbuild import Multiloop
    import prepareChain

    from restraints import EDrestraint ; EDrestraint.setPenalty(2.) ; EDrestraint.setScatRad(1.) ;
    esmin, esmax, esmean, rcmult, xscoreCutoff = .000, 5., .0, 5, 0.9

    guidedsampling = None
    multiPrepC = prepareChain.PrepareChain("SCL1.0")

    xscorer = XrayScorer(None, xscoreCutoff)
    numRefCycles = 20 ; startCycle = 0
    for cycle in range(startCycle, numRefCycles) :
        print "***********************************************"
        print "********************CYCLE %d*******************" % cycle
        print "***********************************************"
        modelIn = "model%d.pdb" % cycle
        cnsout = "cns%d.pdb" % cycle
        rtkmodel = "model%d.pdb" % (cycle+1) # rappertk model to be generated in this cycle
        sfall(modelIn, "rfree.mtz", "phased.mtz")

        phasedmtz = "phased%d.mtz" % cycle # phase the str factors with current model

        if not os.path.isfile(cnsout) :
            cnsRefinement("phased.mtz", modelIn, phasedmtz, cnsout,
                options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg, options.resolution,
                cnsArgs, cycle)
            if cycle > 0 : makeRMSDcomments("cns0.pdb", cnsout, options.loopres)
            #from pdbr import protein
            #refprot = protein("cns0.pdb", read_hydrogens=0, read_waters=0, read_hets=0)
            #testprot = protein(cnsout, read_hydrogens=0, read_waters=0, read_hets=0)
            #if cycle != 0 : ## annotate with RMSD
                #chi1accu, chi1diff, chi12accu, chi12diff, carmsd, mcrmsd, aarmsd = comparePhiPsiOmegaChi(refprot, testprot)
                #remark = "REMARK chain-comparison cns0 ca %6.3f mc %6.3f aa %6.3f chi1 %6.3f %6.3f chi12 %6.3f %6.3f\n" % (carmsd, mcrmsd, aarmsd, chi1accu, chi1diff, chi12accu, chi12diff)
                #if options.loopres :
                    #chi1accu, chi1diff, chi12accu, chi12diff, carmsd, mcrmsd, aarmsd = comparePhiPsiOmegaChi(refprot, testprot, None, options.loopres)
                    #remark += "REMARK  loop-comparison cns0 ca %6.3f mc %6.3f aa %6.3f chi1 %6.3f %6.3f chi12 %6.3f %6.3f\n" % (carmsd, mcrmsd, aarmsd, chi1accu, chi1diff, chi12accu, chi12diff)
                #lines = [remark,] + open(cnsout, 'r').readlines()
                #ofp = open(cnsout,'w')
                #for l in lines : ofp.write(l)
                #ofp.close()

        if cycle == 15 : options.scRad *= 2 ### XXX

        badmcsc = xscorer.score(cnsout, phasedmtz+"2fofc.map", phasedmtz+"fc.map", "FC", "PHIC", "2F1-F2", "mcsc")
        badmc = xscorer.score(cnsout, phasedmtz+"2fofc.map", phasedmtz+"fc.map", "FC", "PHIC", "2F1-F2", "mc")
        badpept = xscorer.score(cnsout, phasedmtz+"2fofc.map", phasedmtz+"fc.map", "FC", "PHIC", "2F1-F2", "pept")
        badscs = xscorer.score(cnsout, phasedmtz+"2fofc.map", phasedmtz+"fc.map", "FC", "PHIC", "2F1-F2", "sc")
        badresids = list ( set(badmcsc+badmc+badpept) )

        #badscs = xscorer.score(cnsout, phasedmtz, "FP", "FC", "PHIC", "2F1-F2", "sc")

        if cycle > 0 and len(badresids) < 1 and len(badscs) < 0 : print "DONE refinement" ; sys.exit(0)
        #userpause()

        xrayRestGen = []
        if options.loopres and options.framework : xrayRestGen.append( prepareChain.CArestraintGenerator(options.loopres, 1., 3.) )
        #if cycle > 0 : xrayRestGen.append( prepareChain.XrayRestraintsGenerator(phasedmtz, "FP", "FC", "PHIC", "2F1-F2", esmin, esmax, esmean ) )
        if cycle > 0 : xrayRestGen.append( prepareChain.XrayRestraintsGenerator(phasedmtz+"2fofc.map", "map", "FC", "PHIC", "2F1-F2", esmin, esmax, esmean, ["SCLbuilder","ChiBuilder","CBbuilder"], ) )
        if cycle == 0 : badresids = None #options.loopres
        if options.loopres and not options.framework :
            if badresids == None : badresids = options.loopres
            else : badresids = list( set(badresids).intersection(options.loopres) )
        if cycle > 0 and options.noRTK == 1 : shutil.copy(cnsout, rtkmodel) ; continue
        if os.path.isfile(rtkmodel) : continue
        scvdwr = options.scReduction ; popsize = options.popsize
        if cycle == 0 : scvdwr = .75 ; popsize = 500
        ml = Multiloop(cnsout, badresids, None, options.caRad, options.scRad, scvdwr, guidedsampling, popsize,
            options.backtrack, 1, "pre"+rtkmodel, xrayRestGen, multiPrepC)
        if cycle > 0 :
            ml.ranker = XrayRanker(phasedmtz+"2fofc.map", "map", "FC", "PHIC", "2F1-F2", esmin, esmax)
            ml.ranker.rankChildren = rcmult ; ml.ranker.rankRecurse = 1
            ml.ranker.rankLeaderBuilderOnly = None ; ml.ranker.rankGivenEnsemble = None
            #if cycle < 5 : ml.ranker.rankLeaderBuilderOnly = 1
            #if cycle < 5 : ml.ranker.rankChildren *= 2
        ml.cellsym = [ options.a, options.b, options.c, options.alpha, options.beta, options.gamma, sgtable[options.sg][0] ]
        nb = ml.run()
        #sys.exit(0)
        #badresids = ml.badresids
        badresids = list(  set(  findChangedSC(cnsout, "pre"+rtkmodel) + badscs )  )
        print "multiloop built", nb, "model/s"
        dontusesc = 1
        if cycle == 0 or dontusesc==1 : os.rename("pre"+rtkmodel, rtkmodel)
        else :
            scPrepC, useGivenRot, useDEE = prepareChain.PrepareChain("PRL"), 1, 1
            SCplacement("pre"+rtkmodel, options.scReduction, rtkmodel, "dotfile", useDEE, phasedmtz+"2fofc.map", "FP", "FC", "PHIC", "2F1-F2", esmin, esmax, None, useGivenRot, badresids, scPrepC).run()
            adjustBfac(rtkmodel, cnsout)
        removeMODEL(rtkmodel)
        copyNonprotein(cnsout, rtkmodel) ## automatic water addition can go here
        #replaceWaters(rtkmodel, rtkmap, 3)

## remove existing waters and add new waters using findwaters, use chain W
def replaceWaters(pdbfile, ccp4map, peakCutoff=2.5) :
    from pdbr import isPdbAAline, isPdbAtomLine, line2resn
    lines = []
    for l in open(pdbfile, 'r').readlines() :
        if isPdbAtomLine(l) and line2resn(l) == "HOH" : continue
        if l ==  "END\n" or l[0:4] == "END " : continue
        lines.append(l)
    opf = open(pdbfile, 'w')
    for l in lines : opf.write(l)
    opf.close()
    cmd = "findwaters --pdbin %s --map %s --pdbout %s --sigma %f" % (pdbfile, ccp4map, "findwat", peakCutoff)
    execCmd(cmd, [])
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
    for l in open(refpdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if id2bo.has_key(line2atomid(l)) : continue
        id2bo[line2atomid(l)] = [line2bfac(l), line2occu(l), line2crdstr(l)]
    newlines = []
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
