import os, shutil, sys, random, re
from xray import cif2mtz, sfall, uniqueify, findResolutionReflections
from pref import cnsRefinement, cnsArgs, changeBfacs, removeMODEL, adjustBfacOccu, copyNonprotein, removeRElines, randomize
from xcheck import XrayScorer
from pdbr import protein, isAAres, line2resid
import prot2res
from peptidebuild import ModelRenderer
from builders import VecInt, VecVecFloat, CBbuilder
from data import consts
from prepareChain import buildCBs
from pref import parseResnums

def parseDelres(pdbfile,loopres):
    loopresids = []
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    lines = open(loopres, 'r').readlines()
    import re
    for li in range(len(lines)) :
        l = re.sub("\n", "", lines[li])
        if l == '' or l[0] == '#' : continue
        if  "delete" in l :
            rids,chain = parseResnums(l)
            for i in range(int(rids[0]),int(rids[1])+1) :
                for k , v in resnums.items():
                    if int(v) == i and chids[k] == chain :
                        loopresids.append(resids[k])
    print loopresids
    return loopresids
    

def deleteResidues(pdbfile, resdel) :
    print resdel
    lines = []
    for l in open(pdbfile, 'r').readlines() :
        if l[0:6] in ["ATOM  ","HETATM"] and line2resid(l) in resdel :
            print "Removing" , line2resid(l)
            continue
        lines.append(l)
    ofp = open(pdbfile, 'w')
    for l in lines :
        ofp.write(l)
    ofp.close()
    

def pdbNumAtoms(pdbfile) :
    na = 0
    for l in open(pdbfile, 'r').readlines() :
        if l[0:6] in ["ATOM  ","HETATM"] : na += 1
    return na

def findBadResids(pdbfile, mtzfile) :
    from xcheck import main as xcheckMain
    xfitdata = xcheckMain(pdbfile, mtzfile, "FP", "FC", "PHIC", "2F1-F2", aasc=1) ## assess bad fit
    badresids, minbadres = None, 10
    for cut in range(80,81, 5) :
        badresids = []
        for k,v in xfitdata.items() :
            if k[0:3] != "HOH" and v < cut/100. : badresids.append(k)
        if len(badresids) >= minbadres : print "CUTOFF", cut/100. ; break
    for br in badresids : print "[%s]" % br
    if len(badresids) < minbadres : return None
    return badresids

def keepMClinesOnly(pdbfn) :
    from pdbr import isPdbAAline, isPdbAtomLine, line2atomname
    lines = []
    for l in open(pdbfn, 'r').readlines() :
        if not isPdbAtomLine(l) : lines.append(l)
        elif isPdbAAline(l) and line2atomname(l) in [" N  "," CA "," C  "," O  "] : lines.append(l)
    ofp = open(pdbfn, 'w')
    for l in lines : ofp.write(l)
    ofp.close()


def keepMClinesOnly2(pdbfn,pdbout) :
    from pdbr import isPdbAAline, isPdbAtomLine, line2atomname
    lines = []
    for l in open(pdbfn, 'r').readlines() :
        if not isPdbAtomLine(l) : lines.append(l)
        elif isPdbAAline(l) and line2atomname(l) in [" N  "," CA "," C  "," O  "," CB "] : lines.append(l)
    ofp = open(pdbout, 'w')
    for l in lines : ofp.write(l)
    ofp.close()

def main() :
    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='to create all the files during refinement. it shdnt be already present.')
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='starting pdb containing a model of pdb-ligand complex')
    parser.add_option("--offset", action='store', type='int', dest='offset', help='a +/- integer that describes the sequence shift to apply', default=0)
    parser.add_option("--truncateN", action='store', type='int', dest='truncateN', help='remove residues from Nterm', default=None)
    parser.add_option("--truncateC", action='store', type='int', dest='truncateC', help='remove residues from Cterm', default=None)
    parser.add_option("--delres", action='store', type='string', dest='delres', help='filename containing residues to be deleted. truncate and this option shdnt be used together', default=None)

    parser.add_option("--sf", action='store', type='string', dest='sf', help='structure factors file')
    parser.add_option("--a", action='store', type='float', dest='a', help='cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='cell spacegroup, in CCP4 notation')
    parser.add_option("--resolution", action='store', type='float', dest='resolution', help='resolution of the data')
    parser.add_option("--extraTOP", action='store', type='string', dest='extraTOP', help='ligand CNS topology file')
    parser.add_option("--extraPAR", action='store', type='string', dest='extraPAR', help='ligand CNS params file')

    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid', default=2)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=1)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--noRTK", action='store', type='int', dest='noRTK', help='dont rebuild bad-fits with rtk', default=None)
    parser.add_option("--uselowres", action='store', type='int', dest='uselowres', help='lower the resolution artificially to have argument*numatoms reflections', default=None)
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='seed for randomizing', default=None)

    (options, args) = parser.parse_args()

    import misc
    misc.setVerbosity(options.verbose)
    from scplacement import SCplacement
    from loopbuild import Multiloop
    import prepareChain
    from restraints import EDrestraint ; EDrestraint.setPenalty(5.) ; EDrestraint.setScatRad(1.) ;
    esmin, esmax = 0.1, 1

    randomize(options.randomize)

    if options.delres != None :
        resdel = parseDelres(options.pdbfile, options.delres)
        print resdel

        ## amk sommmented 
        #for l in open(options.delres, 'r').readlines() :
        #    if l[0] == '#' : continue
        #    l = re.sub("\n", "", l)
        #    resdel.append( re.sub( "]", "", re.sub("\[", "", l) ) )

    #import sys; sys.exit(0)
    if not os.path.isdir(options.scratchdir) : os.mkdir(options.scratchdir)
    shutil.copyfile(options.sf, "%s/strfactors.mtz" % options.scratchdir)
    shutil.copyfile(options.pdbfile, "%s/%s" % (options.scratchdir,options.pdbfile))
    firstmodel = "%s/model0.pdb" % options.scratchdir
    if not os.path.isfile(firstmodel) :
        if options.caRad < 0.49 : shutil.copyfile(options.pdbfile, firstmodel)
        else :
            mlPrepC = prepareChain.PrepareChain("PRL")
            Multiloop(options.pdbfile, None, None, options.caRad, options.scRad, options.scReduction, None, options.popsize, options.backtrack, 1, firstmodel, None, mlPrepC).run()
        if options.delres != None :
            deleteResidues(firstmodel, resdel) #; sys.exit(0)

    os.chdir(options.scratchdir)

    shutil.copyfile("model0.pdb", "model0.AA.pdb")

    numRefCycles = 10
    cif2mtz("strfactors.mtz", "base.mtz", options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg)
    uniqueify("base.mtz", "rfree.mtz")

    if options.noRTK :
        removeRElines("model0.AA.pdb", ["^MODEL", "^ENDMDL"])
        copyNonprotein(options.pdbfile, "model0.AA.pdb")
        for cycle in range(numRefCycles) :
            modelIn = "model%d.AA.pdb" % cycle
            phasedmtz = "phased%d.mtz" % cycle # phase the str factors with current model
            sfall(modelIn, "rfree.mtz", phasedmtz)
            modelOut = "model%d.AA.pdb" % (cycle+1) # model to be generated in this cycle
            mtzout = "phased%d.mtz" % (cycle+1) # model to be generated in this cycle
            cnsRefinement(phasedmtz, modelIn, mtzout, modelOut,
                options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg, options.resolution,
                cnsArgs, cycle, options.extraTOP, options.extraPAR)
        sys.exit(0)

    keepMClinesOnly("model0.pdb")

    if options.offset != None : ## offset=0 by default, not None
        prot = protein("model0.pdb", read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        keys = list(resns.keys()) ; keys.sort()
        newresns = {}
        allAA = [ "GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR", "TRP", "ASP", "GLU", "HIS", "LYS", "ARG", "CYS", "MET", "SER", "THR", "ASN", "GLN", ]
        if options.offset > 0 :
            for ki in range(options.offset, len(keys)) : newresns[ keys[ki] ] = resns[ keys[ki-options.offset] ]
            for ki in range(options.offset) : newresns[ keys[ki] ] = allAA[ int(len(allAA)*random.random()) ]
        elif options.offset < 0 :
            options.offset *= -1
            for ki in range(len(keys)-options.offset) : newresns[ keys[ki] ] = resns[ keys[ki+options.offset] ]
            for ki in range(len(keys)-options.offset, len(keys)) : newresns[ keys[ki] ] = allAA[ int(len(allAA)*random.random()) ]
        else :
            for k,v in resns.items() : newresns[k] = v
        for ri in keys : print resns[ri], newresns[ri]
        ModelRenderer(res, newresns, chids, resnums, inscodes, [], "model0.pdb").render(pts)
        buildCBs("model0.pdb")
        removeRElines("model0.pdb", ["^MODEL", "^ENDMDL"])

    if options.truncateN != None or options.truncateC != None : ## cut out residues from N/C terminii
        prot = protein("model0.pdb", read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        keys = resns.keys(); keys.sort()
        if options.truncateN : keys = keys[options.truncateN:]
        if options.truncateC : keys = keys[0:len(keys)-options.truncateC]
        newres, newresids, newresnums, newresns, newchids, newinscodes, newpts = {}, {}, {}, {}, {}, {}, []
        for ri in keys :
            newri = ri - options.truncateN
            newresids[newri], newresnums[newri], newresns[newri], newchids[newri], newinscodes[newri] \
                    = resids[ri], resnums[ri], resns[ri], chids[ri], inscodes[ri]
            newres[newri] = {}
            for an,ai in res[ri].items() :
                newpts.append( pts[ res[ri][an] ] )
                newres[newri][an] = len(newpts)-1
        res, resids, resnums, resns, chids, inscodes, pts = newres, newresids, newresnums, newresns, newchids, newinscodes, newpts
        ModelRenderer(res, newresns, chids, resnums, inscodes, [], "model0.pdb").render(pts)

    if options.uselowres != None :
        numatoms = pdbNumAtoms("model0.AA.pdb")
        effres = findResolutionReflections( "rfree.mtz", options.uselowres*pdbNumAtoms("model0.AA.pdb") )
        if effres != None : options.resolution = effres
        print "Effective resolution", options.pdbfile, options.resolution, numatoms

    xscorer = XrayScorer(None, 0.9)

    for cycle in range(numRefCycles) :
        print "***********************************************"
        print "********************CYCLE %d*******************" % cycle
        print "***********************************************"
        modelIn = "model%d.pdb" % cycle
        phasedmtz = "phased%d.mtz" % cycle # phase the str factors with current model
        sfall(modelIn, "rfree.mtz", phasedmtz, options.resolution)
        rtkmodel = "rtk%d.pdb" % cycle # rappertk model to be generated in this cycle

        if cycle == 0 : badresids = None ; addSC = 1
        else :
            addSC = None
            badresids = xscorer.score(modelIn, phasedmtz, aasc="sc") #findBadResids(modelIn, phasedmtz)
            #exclres = [ "GLN   34 ", "GLU   37 ", "GLN   48 ", "GLU   50 ", "ASN   53 ", "LYS   74 ", "LYS   77 ",
            #    "LYS   88 ", "LYS   93 ", "ASN  107 ", "ASN  108 ", ]
            #for er in exclres :
            #    if er in badresids : badresids.remove(er)
            if not badresids : print "REFINEMENT DONE"; sys.exit(0)

        scPrepC = prepareChain.PrepareChain("PRL")
        useGivenRot, useDEE = None, 1
        if cycle > 5 : useGivenRot = 1
        scp = SCplacement(modelIn, options.scReduction, rtkmodel, "dotfile", useDEE, phasedmtz, "FP", "FC", "PHIC", "2F1-F2", esmin, esmax, addSC, useGivenRot, badresids, scPrepC)
        scp.makeConfidentDecisions = None
        print scp.run()
        removeMODEL(rtkmodel)
        if cycle > 0 :
            adjustBfacOccu(rtkmodel, modelIn)
            if options.uselowres == None : copyNonprotein(modelIn, rtkmodel)
        else :
            changeBfacs(rtkmodel, [20.,30.])
            if options.uselowres == None : copyNonprotein(options.pdbfile, rtkmodel)

        modelOut = "model%d.pdb" % (cycle+1) # model to be generated in this cycle
        mtzout = "phased%d.mtz" % (cycle+1) # model to be generated in this cycle
        cnsRefinement(phasedmtz, rtkmodel, mtzout, modelOut,
            options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg, options.resolution,
            cnsArgs, cycle, options.extraTOP, options.extraPAR)

if __name__ == "__main__" :
    main()
    import sys ; sys.exit(0)
    replaceWaters("model1.pdb", "rtk0.map")
