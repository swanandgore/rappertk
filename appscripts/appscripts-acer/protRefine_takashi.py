import os, shutil, re
from xray_takashi_ver001 import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_generate_multi_chain, cns_anneal, sgCCP4toCNS, fft, refmac, phenix, cns_cis_peptide, cad, cadnat 
from procrun import proc_run_exitOnError as execCmd
from xcheck import XrayScorer, XrayRanker

cnsArgs = {
    0: {"wa":-1},
    1: {"wa":-1},
    2: {"wa":-1},
    3: {"wa":-1},
    4: {"wa":-1},
    5: {"wa":-1},
    6: {"wa":-1},
    7: {"wa":-1},
    8: {"wa":-1},
    9: {"wa":-1},
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


## add the cryst card for refmac for rtk output pdb NDF 18/01/2007 ##
def add_cryst_card(pdbfile, a, b, c, alpha, beta, gamma, sg):
    cryst_card_line = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %s\n" % (a, b, c, alpha, beta, gamma, sg)   
    lines = []
    lines.append(cryst_card_line)
    for l in open(pdbfile, 'r').readlines() :
        lines.append(l)
    op = open(pdbfile, 'w')
    for l in lines : op.write(l)
    op.close()
       

## change ILE:CD to ILE:CD1 if necessary
## OT1 to O
def fixCNSop(pdbfile) :
    from pdbr import isPdbAAline, isPdbAtomLine, line2atomname, line2resn
    lines = []
    for l in open(pdbfile, 'r').readlines() :
        if not isPdbAAline(l) : lines.append(l)
        elif line2resn(l) == "ILE" and line2atomname(l) == " CD " : lines.append(re.sub(" CD ", " CD1", l))
        elif line2atomname(l) == " OT1" : lines.append(re.sub(" OT1", " O  ", l))
        elif line2atomname(l) == " OXT" : continue
        else : lines.append(l)
    op = open(pdbfile, 'w')
    for l in lines : op.write(l)
    op.close()

## take a multimodel file and split by model removing start tags and putting END tag at end of file 15/01/2007 ##
def removeMODEL_multi(pdbfile, mod_num) :
    mod_num += 1
    newlines = []
    line_collect = "off"
    mod_num_string = "%d" % mod_num
    search_string = "MODEL        "+mod_num_string
    for l in open(pdbfile+".pdb",'r').readlines() :
        if l[0:14] == search_string:
            line_collect = 'on'
        elif l[0:6] == 'ENDMDL' and line_collect == 'on':
            newlines.append("END")
            line_collect = "off"
        elif line_collect == "on"  : newlines.append(l)
    op = open(pdbfile+mod_num_string+".pdb", 'w')
    for l in newlines : op.write(l)
    op.close()

## assert that there is only 1 model, and remove model, endmdl line
def removeMODEL(pdbfile) :
    newlines = []
    line_collect = "off"
    search_string = "MODEL"
    for l in open(pdbfile,'r').readlines() :
        if l[0:5] == search_string:
            line_collect = 'on'
        elif l[0:6] == 'ENDMDL' and line_collect == 'on':
            newlines.append("END")
            line_collect = "off"
        elif line_collect == "on"  : newlines.append(l)
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
    from pdbr import isPdbAAline, isPdbAtomLine
    lines = open(topdb, 'r').readlines()
    assert lines[ len(lines)-1 ][0:3] == "END"
    lines[ len(lines)-1 ] = "TER\n"
    print "copying nonprotein atoms from", frompdb, "to", topdb
    for l in open(frompdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if isPdbAAline(l) : continue
        if not waterAlso and line2resn(l) == "HOH" : continue
        lines.append(l)
    op = open(topdb, 'w')
    for l in lines : op.write(l)
    op.close()

## assign the occupancies/bfactors in pdbfile according to refpdb.
## use resid+atomname as the key
## copy the atoms present in refpdb but absent in pdbfile
def adjustBfacOccu(pdbfile, refpdb) :
    print "Copy Bfactors from", refpdb, "to", pdbfile
    from pdbr import line2bfac, isPdbAtomLine, line2occu, line2bfac, line2atomid, changeBfactorOccupancy
    id2bo = {}
    for l in open(refpdb, 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        if id2bo.has_key(line2atomid(l)) : continue
        id2bo[line2atomid(l)] = [line2bfac(l), line2occu(l)]
    newlines = []
    for l in open(pdbfile, 'r').readlines() :
        if not isPdbAtomLine(l) : newlines.append(l) ; continue
        bfac, occu = id2bo[ line2atomid(l) ]
        l = changeBfactorOccupancy(l, bfac, occu)
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

def molProbity_badres(pdbfile, outfile, badkeys, cycle):

    from xray_takashi_ver001 import run_reduce, run_molprobity

    tmpout = "out_%d.pdb" % cycle
    run_reduce(pdbfile, tmpout)
    run_molprobity(tmpout, outfile)

    outf = open(outfile, "r")
    lines = outf.readlines()
    outf.close()
    
    for line in lines:
        #convert file to restk style
        if line[0] != "#":
            lsp = line.split(":")
            key = lsp[1][6:9]+lsp[1][0]+lsp[1][1:5]+" "
            print "--"+key+"--"
            badkeys.append(key)
        
    return badkeys

def cnsRefinement(mtzin,pdbin, mtzout,pdbout, a,b,c,alpha,beta,gamma,sg,reso, cnsArgs,cycle,multichain) :
    if os.path.isfile(pdbout) : return
    mtz2hkl(mtzin, "cns.hkl")
    if multichain == 1:
        cns_generate_multi_chain(pdbin, "generate.mtf", "generate.pdb", None, None, "generate%d.log"%cycle)
    else:
        cns_generate(pdbin, "generate.mtf", "generate.pdb", None, None, "generate%d.log"%cycle)
    removeZeroLines("generate.pdb")
    cns_cis_peptide()
    wa = -1 #cnsArgs[cycle]["wa"]
    cns_anneal(a, b, c, alpha, beta, gamma, sgCCP4toCNS[sg], reso,
        "cns.hkl", "generate.mtf", "generate.pdb", None, "anneal%d.log"%cycle, wa)
    removeZeroLines("anneal.pdb")
    fixCNSop("anneal.pdb")
    os.rename("anneal.pdb", pdbout)
    sfall(pdbout, "rfree.mtz", mtzout)
    #moleman(pdbout) 

def main() :
    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='to create all the files during refinement.')
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='starting pdb containing a model of pdb-ligand complex')
    parser.add_option("--sf", action='store', type='string', dest='sf', help='structure factors file')
    parser.add_option("--a", action='store', type='float', dest='a', help='cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='cell spacegroup, in CCP4 notation')
    parser.add_option("--resolution", action='store', type='float', dest='resolution', help='resolution in angstroms')
#
#12052011:Takashi:Changed default value to 10 %
#
    parser.add_option("--percent-r-free", action='store', type='string', dest='prfree', help='percentage R-free to take must be at least 500 reflections', default=0.1)
#
#13052011:Takashi:Added labels of amplitude and sigma of structure factors
#
    parser.add_option("--labFin", action='store', type='string', dest='fin', help='label of structure factor amplitude', default="FP")
    parser.add_option("--labSIGFin", action='store', type='string', dest='sigfin', help='label of sigma of structure factor amplitude', default="SIGFP")
    parser.add_option("--n-models", action='store', type='int', dest='nmodel', help='the number of models required to make up the ensemble')
    parser.add_option("--restart-n", action='store', type='int', dest='restart', help='the CNS round to restart from', default=0)
    parser.add_option("--restart-n-mod", action='store', type='int', dest='restart_mod', help='the model round to restart from', default=0)
    parser.add_option("--refine-prog", action='store', type='string', dest='refinetype', help='Use CNS or PHENIX', default="CNS")
#
#12052011:Takashi:Added --phenix-param
#
    parser.add_option("--phenix-param", action='store', type='string', dest='paramin', help='paramter file for phenix.refine')
    parser.add_option("--multichain", action='store', type='int', dest='multichain', help='Does the target have more than one chain if so then use 1', default=0)
    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid', default=2)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=1)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--ranseed", action='store', type='float', dest='ranseed', help='Explicitly set the random seed.', default=1)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--noRTK", action='store', type='int', dest='noRTK',help='dont rebuild bad-fits with rtk', default=0)
    (options, args) = parser.parse_args()

    import misc
    misc.setVerbosity(options.verbose)
    misc.RanGen.instance().seedme( int(-1000 * options.ranseed) )
    
    #import sys; sys.exit(0)

    if not os.path.isdir(options.scratchdir) :
        os.mkdir(options.scratchdir)
	#
        #13052011:Takashi:Added the next line to convert input labels to FP and SIGFP
        #
        cadnat(options.sf, "%s/cadnat.mtz" % options.scratchdir, options.fin, options.sigfin)
	shutil.copyfile(options.pdbfile, "%s/PDB.pdb" % options.scratchdir)
	#
	#13052011:Takashi:Commented out the next line and added the following line
	#
        #shutil.copyfile(options.sf, "%s/strfactors.mtz" % options.scratchdir)
	shutil.copyfile("%s/cadnat.mtz" % options.scratchdir, "%s/strfactors.mtz" % options.scratchdir)
    os.chdir(options.scratchdir)

    from xcheck import main as xcheckMain
    from loopbuild import Multiloop
    from scplacement import SCplacement
    import prepareChain
    from restraints import EDrestraint ; EDrestraint.setPenalty(5.) ; EDrestraint.setScatRad(1.) ;
    esmin, esmax, esmean, rcmult = 0.0001, 5.0, .1, 10
    EDrestraint.setScatRad(.5) ;
    multiPrepC = prepareChain.PrepareChain("PRL1.0") ##here we can pass sc lib option

    if options.restart == 0 and options.restart_mod == 0:
        Multiloop("PDB.pdb", None, None, options.caRad, options.scRad, options.scReduction, None , options.popsize, options.backtrack, options.nmodel, "model.pdb", None, multiPrepC).run()
        ## change to remove n models requested NDF 15/01/07 ###
        for mod_num in range(options.nmodel):
            print "Model file manipulation on %d .pdb" %mod_num 
            removeMODEL_multi("model", mod_num)
            mod_file_num = mod_num + 1
            mfn_st = "%d" % mod_file_num
            copyNonprotein("PDB.pdb", "model"+mfn_st+".pdb")
    #
    #12052011:Takashi:commented out cif2mtz and added the next line
    #
    #cif2mtz("strfactors.mtz", "base.mtz", options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg)
    shutil.copyfile("%s/strfactors.mtz" % options.scratchdir, "%s/base.mtz" % options.scratchdir)
    uniqueify("base.mtz", "rfree.mtz", options.prfree) ## add here the percentage need to generate r-free - note for low res need to be 10% otherwise 5%

    
    ## for each model generate a subdir and copy in model then run script in each sub file NDF 15/01/07 ##
    for mod_dir in range(options.restart_mod, options.nmodel):
        mod_num_st = "%d" % (mod_dir + 1)
        if not os.path.isdir("model_"+mod_num_st+"_refine"): 
            os.mkdir("model_"+mod_num_st+"_refine")
            shutil.copyfile("model"+mod_num_st+".pdb", "model_"+mod_num_st+"_refine/model"+mod_num_st+"_0.pdb")
	    #
	    #12052011:Takashi:Removed "../"
	    #13052011:Takashi:Added cadnat line
	    #
            #shutil.copyfile(options.sf, "model_"+mod_num_st+"_refine/strfactors.mtz")
            cadnat(options.sf, "%s/cadnat.mtz" % options.scratchdir, options.fin, options.sigfin)
	    shutil.copyfile("cadnat.mtz", "model_"+mod_num_st+"_refine/strfactors.mtz")
            shutil.copyfile("rfree.mtz", "model_"+mod_num_st+"_refine/rfree.mtz")
            shutil.copyfile("base.mtz", "model_"+mod_num_st+"_refine/base.mtz")
        ##if we get a cluster working here we want to paralalise the job NDF 17/01/2007 ##
        
        os.chdir("model_"+mod_num_st+"_refine")
        print "HAVE MOVED TO model_"+mod_num_st+"_refine" 
            
        justcns = None
        numRefCycles = 10
        
        for cycle in range(options.restart, numRefCycles):
	#
	#13052011:Takashi:changed CNS to Refine
	#
            print "********************************************************"
            print "**************** Refinement CYCLE %d ********************" % cycle
            print "********************************************************"
            
            xscorer = XrayScorer(None, 0.90)

            model_stem = "model"+mod_num_st+"_%d" % cycle
            modelIn = "model"+mod_num_st+"_%d.pdb" % cycle
            cnsout = "ref%d.pdb" % cycle
            rtkmodel = "model"+mod_num_st+"_%d.pdb" % (cycle+1) # rappertk model to be generated in this cycle
            cnsdiffmap = "cns%ddiff.map" % cycle
            
            if cycle >= 8:
                multiPrepC = prepareChain.PrepareChain("SCL0.5") 
            else:
                multiPrepC = prepareChain.PrepareChain("PRL1.0") 

            #
	    #13052011:Takashi:Modified sfall to be able specif labels of structure factors
	    #
            sfall(modelIn, "rfree.mtz", "phased.mtz",  "FP", "SIGFP")
            phasedmtz = "phased_"+mod_num_st+"_%d.mtz" % cycle # create a phased 2Fo-Fc map from current model and str factors
	    #
	    # choosing cns or phenix here
	    #
            if options.refinetype == "CNS":
                cnsRefinement("phased.mtz", modelIn, phasedmtz, cnsout,
                              options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg, options.resolution,
                              cnsArgs, cycle, options.multichain)
            elif options.refinetype == "PHENIX":
                phenix(modelIn, "phased.mtz", options.paramin)
                os.rename(model_stem+"_refine_001.pdb", cnsout)           
		#
		#13052011:Takashi:added the next line because output labels of phenix are different from other programs
		#
                cad("model"+mod_num_st+"_%d_refine_%03d.mtz" % (cycle, cycle+1), phasedmtz, 
		    "F-obs", "SIGF-obs", "F-model", "PHIF-model", "FP", "SIGFP")

            fft(phasedmtz, cnsdiffmap, "FP", "FC", "PHIC", "SIGFP")

            #badresids = xscorer.score(cnsout, phasedmtz, "FP", "FC", "PHIC", "2F1-F2")
	    badresids = xscorer.score(cnsout, phasedmtz)

            if cycle > 0 and len(badresids) < 1 : print "DONE refinement" ; sys.exit(0)
            xrayRestGen = []
            
            #xrayRestGen.append( prepareChain.XrayRestraintsGenerator(phasedmtz, "FP", "FC", "PHIC", "2F1-F2", 0.1, 2, 1.0 ) ) ##hard
            addSC = None
            
            status_ml = Multiloop(cnsout, badresids, addSC, options.caRad, options.scRad, options.scReduction, None , options.popsize, options.backtrack, 1, "pre"+rtkmodel, xrayRestGen,multiPrepC)
            
            status_ml.ranker = XrayRanker(phasedmtz, "FP", "FC", "PHIC", "2F1-F2", esmin, esmax)
            status_ml.ranker.rankGivenEnsemble = None
            status_ml.ranker.rankChildren = rcmult
            status_ml.ranker.rankRecurse = 1
            nb = status_ml.run()
            #badresids = status_ml.badresids
            badresids = findChangedSC(cnsout, "pre"+rtkmodel)
            print "multiloop built", nb, "model/s"
                        
            
            SCplacement("pre"+rtkmodel, 0.5 , rtkmodel, "dotfile", None, phasedmtz, "FP", "FC", "PHIC", "2F1-F2", esmin, esmax, None, 1, badresids,multiPrepC).run()
            adjustBfacOccu(rtkmodel, cnsout)
                        
            removeMODEL(rtkmodel)
            copyNonprotein(cnsout, rtkmodel)
            rtkmap = phasedmtz+"FCFPPHIC_B.map"
            cns_round_final_pdb = "model"+mod_num_st+"_%d" % (cycle+1)
            cns_final_mtz = phasedmtz

        if options.restart == 10:
            cns_round_final_pdb = "model"+mod_num_st+"_10"
            cns_final_mtz = "phased_"+mod_num_st+"_9.mtz"
        if options.restart_mod != 0:
            options.restart = 0
            
        ## add refamc *5 cycles NDF 18/01/2007 ##
        refmac_cycle = 5
        pdbin = ""
        for rcycle in range(refmac_cycle):
            #set up files
            xscorer = XrayScorer(None, 0.90)
            rtkmodel = "model"+mod_num_st+"_%d" % (rcycle+1)
            hklout = "phased_"+mod_num_st+"_9_refmac%d.mtz" % (rcycle+1)
            multiPrepC = prepareChain.PrepareChain("SCL0.2") ##here we can pass sc lib option
            if rcycle == 0:
                pdbin = cns_round_final_pdb+".pdb"
                pdbout = cns_round_final_pdb+"_refmac%d.pdb" % (rcycle+1)
                hklin = cns_final_mtz
            else:
                pdbin = cns_round_final_pdb+"_refmac%d.pdb" % (rcycle)
                pdbout = cns_round_final_pdb+"_refmac%d.pdb" % (rcycle+1)
                hklin = "phased_"+mod_num_st+"_9_refmac%d.mtz" % (rcycle)
            add_cryst_card(pdbin, options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg)
            moleman(pdbin)
            
            print "***********************************************"
            print "******************* Refmac CYCLE %d************" % rcycle
            print "***********************************************"
            refmac(hklin, hklout, pdbin , pdbout )

            badresids_xray = xscorer.score(pdbout, hklout, "FP", "FC", "PHIC", "2F1-F2")

            molp_outfile = "molp_%d.out" % rcycle
            badresids = molProbity_badres(pdbout, molp_outfile, badresids_xray, rcycle)            
            
            if rcycle > 0 and len(badresids) < 1 : print "DONE refinement" ; sys.exit(0)
            xrayRestGen = []
            xrayRestGen.append( prepareChain.XrayRestraintsGenerator(hklout, "FP", "FC", "PHIC", "2F1-F2", 0.1, 2, 1.0 ) ) ##hard

            addSC = None
            status_ml = Multiloop(pdbout, badresids, addSC, options.caRad, options.scRad, 0.75 , None , options.popsize, options.backtrack, 1, rtkmodel, xrayRestGen,multiPrepC)
            status_ml.ranker = XrayRanker(hklout, "FP", "FC", "PHIC", "2F1-F2", esmin, esmax)
            status_ml.ranker.rankGivenEnsemble = None
            status_ml.ranker.rankChildren = rcmult
            status_ml.ranker.rankRecurse = 1
            nb = status_ml.run()
            #badresids = status_ml.badresids
            #badresids = findChangedSC(pdbout, "pre"+rtkmodel)
                     
            #print "multiloop built", nb, "model/s"
            #scPrepC = prepareChain.PrepareChain("PRL1.0")
            #SCplacement("pre"+rtkmodel, options.scReduction, rtkmodel, "dotfile", None, hklout, "FP", "FC", "PHIC", "2F1-F2", esmin, esmax, None, 1, badresids, scPrepC).run()
                        
            removeMODEL(rtkmodel)
            adjustBfacOccu(rtkmodel, pdbout)
            copyNonprotein(pdbout, rtkmodel)
            rtkmap = hklout+"FCFPPHIC_B.map"
            if rcycle >= 4:
                replaceWaters(rtkmodel, rtkmap)
            rtkpdb_final = rtkmodel
            mtzfinal = hklout

        add_cryst_card(rtkpdb_final, options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg)
        moleman(rtkpdb_final)
        hklout = "RAPPERtk_final_mod"+mod_num_st+".mtz"
        pdbout = "RAPPERtk_final_mod"+mod_num_st+".pdb"
        refmac(hklin, hklout, pdbin , pdbout )    
        
        os.chdir("../")
        
## remove existing waters and add new waters using findwaters, use chain W
def replaceWaters(pdbfile, ccp4map) :
    from pdbr import isPdbAAline, isPdbAtomLine, line2resn
    lines = []
    for l in open(pdbfile, 'r').readlines() :
        if isPdbAtomLine(l) and line2resn(l) == "HOH" : continue
        if l ==  "END\n" or l[0:4] == "END " : continue
        lines.append(l)
    opf = open(pdbfile, 'w')
    for l in lines : opf.write(l)
    opf.close()
    cmd = "findwaters --pdbin %s --map %s --pdbout %s --sigma 2.5" % (pdbfile, ccp4map, "findwat")
    execCmd(cmd, [])
    for l in open('findwat', 'r').readlines() :
        if not isPdbAtomLine(l) : continue
        lines.append(l)
    lines.append("END\n")
    opf = open(pdbfile, 'w')
    for l in lines : opf.write(l)
    opf.close()

if __name__ == "__main__" :
    main()
    import sys ; sys.exit(0)
    #replaceWaters("model1.pdb", "rtk0.map")
