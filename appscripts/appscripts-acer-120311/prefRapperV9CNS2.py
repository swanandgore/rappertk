import os, shutil, re
import geometry
from xray import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS, fft, omitmap, mapman
from procrun import proc_run_exitOnError as execCmd
from xcheck import XrayScorer, XrayRanker
from data import sgtable
from evalCAtrace import comparePhiPsiOmegaChi
from pdbr import protein, isAAres
import prot2res
from pref import removeZeroLines
from pref import fixCNSop



from data import sgtable , long2shortHM
from scplacement import SCplacement
from loopbuild import Multiloop
import prepareChain
from stump import getCRYST , getRESO


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
    parser.add_option("--use-ca-restraints", action='store', dest='caRes', help='[True/False], Apply positional restraints on the C-alpha atoms',default="True")
    parser.add_option("--use-sc-restraints", action='store', dest='scRes',type= 'string', help='[True/False],  Apply positional restraints on the centroid of the sidechain atoms',default="True",)
    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='radius of spherical restraint on CA position', default=1)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='radius of spherical restraint on sidechain centroid', default=2)
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default=1)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--noRTK", action='store', type='int', dest='noRTK', help='dont rebuild bad-fits with rtk', default=0)
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='seed for randomizing', default=None)
    parser.add_option("--loopres", action='store', type='string', dest='loopres', help='filename containing resids for starting perturbation', default=None)
    parser.add_option("--framework", action='store', type='int', dest='framework', help='to be used in conjunction with loopres. it puts a 1/3 ca/sc pos restr on non-loopres and perturbs them too', default=None)


    (options, args) = parser.parse_args()






    if not os.path.isdir(options.scratchdir) : os.mkdir(options.scratchdir)
    shutil.copyfile(options.pdbfile, "%s/0.model0.pdb" % options.scratchdir)
    shutil.copyfile(options.sf, "%s/rfree.mtz" % options.scratchdir)
    os.chdir(options.scratchdir)



    if (options.a == None or options.b == None or options.c == None or options.alpha == None or options.beta == None or options.gamma == None) :
        options.a,options.b,options.c,options.alpha , options.beta , options.gamma,d1  = getCRYST(options.pdbfile)
        
        if (options.a == None or options.b == None or options.c == None or options.alpha== None or options.beta==None or options.gamma == None ):
            print "CRYST card cannot be read from coordinate file. Please input cell paramater a, b , c , alpha, beta , gamma = ",options.a , options.b , options.c , options.alpha , options.beta  , options.gamma 
            import sys ; sys.exit()



            
    if options.sg == None : 
        d1,d2,d3,d4 , d5 , d6, options.sg  = getCRYST(options.pdbfile)
        if options.sg == None : 
            print "Please input space group " , options.sg ; import sys ; sys.exit()
        ss = ""
        for sg1 in options.sg:
            if sg1 in ["\n","\t","\s"]:
                continue
            else :
                ss = ss+sg1
        options.sg = ss
                
    if options.sg  in long2shortHM.keys():
        shortsg = long2shortHM[options.sg]
        options.sg = shortsg
    if options.sg not in sgtable.keys():
        print "Check --sg , Not recognised [%s][%d]"%( options.sg, len(options.sg))
        import sys ; sys.exit()

    if options.resolution == None : 
        options.resolution = getRESO(options.pdbfile)
        if (options.resolution == None):
            print "Please input resolution " , options.resolution
            import sys ; sys.exit()




    rotlib = "SCL1.0"
    numRefCycles = 10 ; startCycle = 0




    for cycle in range(startCycle, numRefCycles) :

        if cycle < 10 :
            xscoreCutoff = 0.8
        else :
            xscoreCutoff = 0.9



        if cycle == 15 : options.scRad *= 2 



        if cycle == 0 :
            scvdwr = .75 ; popsize = 500

        modelIn = "0.model%d.pdb" % cycle
        cnsout = "cns%d.pdb" % cycle
        rtkmodel = "model%d.pdb" % (cycle+1) # rappertk model to be generated in this cycle

        sfall(modelIn, "rfree.mtz", "phased.mtz")
        phasedmtz = "phased%d.mtz" % cycle # phase the str factors with current model


        
        if not os.path.isfile(cnsout) :
            cnsRefinement("phased.mtz", modelIn, phasedmtz, cnsout,
                          options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg, options.resolution,
                          cnsArgs, cycle)



        from pref12 import main as  prefRapperMain
        ## options.pdbfile,options.pdbout,options.dir_xyzout,options.mapfn,options.mtzfn,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,options.nmodels,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freeRlabel,options.n,options.poorOnly,options.poorThreshold,options.loopres,options.start,options.stop,options.chainid,options.modelN2C,options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt)
        
        prefRapperMain(cnsout,rtkmodel,options.scratchdir,None,phasedmtz,options.caRes,options.scRes,options.caRad,options.scRad,scvdwr,popsize,options.verbose,options.backtrack,rotlib, 1 , "False", "False" , "True" , 5 , 5.0 ,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,"FP","SIGFP",None,None,"True","FreeR_flag",2,"True",xscoreCutoff,None,None,None,None,"False","False",20.0,30.0,0.25,2.0,"False","False")
        

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
