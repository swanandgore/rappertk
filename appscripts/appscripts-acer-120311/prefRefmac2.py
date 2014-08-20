import os, shutil, re , sys
import geometry
from xray import cif2mtz, uniqueify, sfall2, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS, fft, omitmap, mapman, sfall
from procrun import proc_run_exitOnError as execCmd
from xcheck import XrayScorer, XrayRanker
from data import sgtable , long2shortHM
from evalCAtrace import comparePhiPsiOmegaChi
from pdbr import protein, isAAres  , line2atomname
import prot2res
from scplacement import SCplacement
from loopbuildV5 import Multiloop, incompleteSCcorrection2
import prepareChainV3
from multProtref import joinPDBs, splitPDBs, restoreChainids
import checkProtChains
import checkProtChainsV2
from peptidebuild import ModelRenderer
import misc
from pdbr import isPdbAAline, isPdbAtomLine, line2resn , line2resnum , line2chid , line2atomnumber
from xray import refmac
from protRefine_nick import molProbity_badres, add_cryst_card

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


def changeBfacs1(filename, bfacs) :
    from pdbr import isPdbAAline, line2atomname, line2occu, changeBfactor
    lines = []

    if (os.path.isfile(filename)==False) :
        print "Cannot find file %s "%filename
        print "No file in directory ", os.getcwd()
        import sys ;        sys.exit()

    for l in open(filename, 'r').readlines() :
        if isPdbAAline(l) :
            if line2atomname(l) in [" N  "," CA "," C  "," O  "] : l = changeBfactor(l, bfacs[0])
            else : l = changeBfactor(l, bfacs[1])
        lines.append(l)
    fp = open(filename, 'w')
    for l in lines : fp.write(l)
    fp.close()

def refmacRefinement(mtzin,pdbin, mtzout,pdbout, a,b,c,alpha,beta,gamma,sg,reso, ccp4Args,cycle) :

    from xray import refmac1
    if os.path.isfile(pdbout) : return
    inpdb, bufpdb = pdbin, "buf.pdb"
    for ccpa in ccp4Args[cycle] :
        if "assignBfac" in ccpa.keys() : changeBfacs1(inpdb, ccpa["assignBfac"])
        refmac1(mtzin, mtzout, inpdb, bufpdb, reso, ccpa["reftype"], ccpa["wa"], ccpa["breftype"], ccpa["ncyc"])
        inpdb = bufpdb ; bufpdb = bufpdb + "1"
    print "renaming", inpdb, pdbout
    os.rename(inpdb, pdbout)


def main() :

    import optparse ; parser = optparse.OptionParser()
    
    
    parser.add_option("--dir-xyzout", action='store', type='string', dest='dir_xyzout', help='to create all the files during refinement. it shdnt be already present.')
    parser.add_option("--xyzin", action='store', type='string', dest='pdbfile', help='starting pdb containing a model of pdb-ligand complex')
    parser.add_option("--hklin", action='store', type='string', dest='mtzfn', help='structure factors file')

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
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='factor to reduce effective vdw dist in case of sidechains', default= 0.75)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='population size for PopulationStrategy', default=100)


    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='use backtracking version of PopulationStrategy. eg 4X5 will set backtrack numsteps and stepsize to 4,5 respectively. not used by default.', default=None)
    parser.add_option("--rotamerlib", action='store', type='string', dest='rotLib', help='[PRL/SCL1.0/SCL0.5/SCL0.2] Name of rotamer library to use when building side chains ', default='SCL1.0')
    parser.add_option("--add-sidechains", action='store', type='string', dest='addsc', help='Build missing side chains ', default='False')
    parser.add_option("--use-given-rotamer", action='store', type='string', dest='userot', help='Use given rotamer', default='False')        
    parser.add_option("--noRTK", action='store', type='int', dest='noRTK', help='dont rebuild bad-fits with rtk', default=0)
    parser.add_option("--randomize", action='store', type='int', dest='randomize', help='seed for randomizing', default=None)

    parser.add_option("--mconly", action='store', type='string', dest='mconly', help='[True/False] Build mainchain only', default="False")
    parser.add_option("--sconly", action='store', type='string', dest='sconly', help='[True/False] Build side chains only, can only be used when MAP/MTZ file is given. See web page for further details', default="False")

    parser.add_option("--opsax", action='store', type='string', dest='opsax', help='[True/False] Reassign side chains with OPSAX, will only be used when MTZ or MAP file is given', default="True")

    parser.add_option("--attempts", action='store', type='int', dest='natt', help='Number of attempts made to build section', default=5)

    parser.add_option("--cacaCutoff", action='store', type='float', dest='cacaCutoff', help='Minimum distance ( angstrom ) between adjacent Calpha atoms in order to detect a chain-break', default=5.)


    #################    Electron density parameters ####################################
    

    
    parser.add_option("--FP", action='store', type='string', dest='f1label', help='Column label for FP in MTZ file', default=None)
    parser.add_option("--SIGFP", action='store', type='string', dest='sigf1label', help='Column label for sigFP in MTZ file', default=None)
    parser.add_option("--FC", action='store', type='string', dest='f2label', help='Column label for FC in MTZ file', default=None)
    parser.add_option("--PHIC", action='store', type='string', dest='phiclabel', help='Column label for PHIC in MTZ file', default=None)
    parser.add_option("--use-FreeR", action='store', type='string', dest='usefreer', help='[True/False] Use FreeR set ? ', default="False")
    parser.add_option("--FreeR", action='store', type='string', dest='freeRlabel', help='Column label for FreeR in MTZ file', default=None)
    parser.add_option("--n", action='store', type='int', dest='n', help='Value of n for difference map calculations nFo-(n-1)Fc', default=2)
    parser.add_option("--mapformat", action='store', type='string', dest='mapformat', help='Value of n for difference map calculations nFo-(n-1)Fc', default="ccp4")

    ############# Residues to be modelled ####################################
    

    parser.add_option("--rebuild-poor-regions-only", action='store', type='string', dest='poorOnly', help='[True/False] Rebuild regions ofinput structure with poor fit to an electron density map. Residues to be rebuilt are identified using a real space correlation coefficientscore, the cut-off for which is set using --poor-fit-threshold.', default="False")
    parser.add_option("--poor-fit-threshold", action='store', type='float', dest='poorThreshold', help='Correlation coefficient threshold to identify poor fitting regions', default=0.9)



    parser.add_option("--loopseq", action='store', type='string', dest='loopres', help='Amino acid sequence for loop to be built', default=None)
    parser.add_option("--use-loopclosure-restraints", action='store', type='string', dest='closure', help='Use geometric restraints to ensure closure of loop with anchor residues', default= "True")
    parser.add_option("--start", action='store', type='int', dest='start', help='Residue number to start building from ', default=None)
    parser.add_option("--stop", action='store', type='int', dest='stop', help='Residue number to stop building at', default=None)
    parser.add_option("--chainid", action='store', type='string', dest='chainid', help='Chain ID of section to be built.', default=None)
    parser.add_option("--modelN2C", action='store', type='string', dest='modelN2C', help='[True/False] Model fragment without loop closure restraints. Used in conjunction with --start, --stop, --chainid. Requires --use-ca-restraints True ', default="False")

    ######### Ouptut parameters #############################################

    parser.add_option("--start-inscode", action='store', type='string', dest='startcode', help='Residue number to start building from ', default=' ')
    parser.add_option("--stop-inscode", action='store', type='string', dest='stopcode', help='Residue number to stop building at', default=' ')
    parser.add_option("--models-get-native-bfactors", action='store', type='string', dest='nativeBfac', help='[True/False] Assign B-factors of remodelled atoms to original values', default="False")
    parser.add_option("--default-mainchain-b-factor", action='store', type='float', dest='mcBfac', help='The value of B-factor assigned to the newly built main chain atoms', default=20.)
    parser.add_option("--default-sidechain-b-factor", action='store', type='float', dest='scBfac', help='The value of B-factor assigned to the newly built side chain atoms', default=30.)



    ### Electron density parametets #########################################

    parser.add_option("--minimum-sig", action='store', type='float', dest='minXSig', help='Minimum sigma ', default=0.25)
    parser.add_option("--maximum-sig", action='store', type='float', dest='maxXSig', help='Maximum sigma ', default=2.0)



    ########## Optional restraints ##########################################
    
    parser.add_option("--make-ed-optional", action='store', type='string', dest='edOpt', help='[True/False]  If False, then the mainchain will be unconditionally forced to lie in positive density. If True then positive density restraint on the mainchain will be made optional.This is useful when tracing through a structure with regions in very poor (non-existent) density', default= "False")

    parser.add_option("--make-all-restraints-optional", action='store', type='string', dest='allOpt', help='[True / False ]  If True, then all  restraints will be made optional', default="False")
    
    

    (options, args) = parser.parse_args()
    


    

    if not os.path.isdir(options.dir_xyzout) : os.mkdir(options.dir_xyzout)
    os.chdir(options.dir_xyzout)


    pdbfilepath = options.pdbfile

    
    shutil.copyfile(options.mtzfn, "%s/init.mtz" % (options.dir_xyzout))
    mtzfilepath = "%s/init.mtz" % (options.dir_xyzout)

    num_refmac_cycles = 15 ; startCycle = 0
    from pref14 import main as  prefRapperMain

    from stump import getCRYST , getRESO
    if options.mtzfn != None :
        if (options.a == None or options.b == None or options.c == None or options.alpha == None or options.beta == None or options.gamma == None) :
            print "Getting cell paramters from coordinate file....."
            options.a,options.b,options.c,options.alpha , options.beta , options.gamma,d1  = getCRYST(options.pdbfile)

            if (options.a == None or options.b == None or options.c == None or options.alpha == None or options.beta == None or options.gamma == None ):
                print "CRYST card cannot be read from coordinate file. Please input cell paramater a, b , c , alpha, beta , gamma = "
                import sys ; sys.exit()
            
        if options.sg == None : 
            print "Getting space group from coordinate file....."
            d1,d2,d3,d4 , d5 , d6, options.sg  = getCRYST(options.pdbfile)
            if options.sg == None : 
                print "Please input space group " ; import sys ; sys.exit()
        ss = ""
        for sg1 in options.sg:
            if sg1 in ["\n","\t","\s"]:
                continue
            else :
                ss = ss+sg1

        options.sg = ss
        if options.sg  in long2shortHM.keys(): shortsg = long2shortHM[options.sg];  options.sg = shortsg
        if options.sg not in sgtable.keys(): print "Check --sg , Not recognised [%s]"%options.sg ;            import sys ; sys.exit()
        print "Setting Space Group to",options.sg
        
        if options.resolution == None : 
            print "Getting resolution limit from coordinate file........"
            options.resolution = getRESO(options.pdbfile)
            if (options.resolution == None):
                print "Please input resolution " , options.resolution ; import sys ; sys.exit()
            print "Resolution = [ " , options.resolution, " ] "
    
    for rcycle in range(num_refmac_cycles):
        if rcycle == 0:
            pdbin = options.pdbfile
            hklin = mtzfilepath
            pdbout = "refmac%d.pdb" % (rcycle+1)
            hklout = "phased_refmac%d.mtz" % (rcycle+1)
            rtkmodel = "model%d.pdb" % (rcycle+1)
        else:
            pdbin = "0.model%d.pdb" % (rcycle)
            hklin = "phased_refmac%d.mtz" % (rcycle)
            pdbout = "refmac%d.pdb" % (rcycle+1)
            hklout = "phased_refmac%d.mtz" % (rcycle+1)
            rtkmodel = "model%d.pdb" % (rcycle+1)


        print "***********************************************"
        print "******************* Refmac CYCLE %d************" % rcycle
        print "***********************************************"
        
        
        if (os.path.isfile(pdbout)==False) : 
            refmacRefinement(hklin,pdbin, hklout,pdbout, options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution, ccp4args,rcycle) 
            #refmac(hklin, hklout, pdbin , pdbout )
            if (os.path.isfile(pdbout)==False) :
                print "Refmac cycle %d failed, no %s generated"%(rcycle,pdbout)
                import sys ; sys.exit()
            if (os.path.isfile(hklout)==False) :
                print "Refmac cycle %d failed, no %s generated"%(rcycle,hklout)
                import sys ; sys.exit()
        
                
        prefRapperMain(pdbout,rtkmodel,options.dir_xyzout,None,hklout,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,1,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freeRlabel,options.n,options.poorOnly,options.poorThreshold,options.loopres,options.start,options.startcode,options.stop,options.stopcode,options.chainid,options.modelN2C,options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt,options.closure,options.addsc,options.userot,options.mapformat)
        
        from copyheader import copyheader
        copyheader("0."+rtkmodel,options.pdbfile)

def molProbity_badres(pdbfile, outfile, badkeys, cycle):

    from xray import run_reduce, run_molprobity

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





if __name__ == "__main__" :
    main()
    import sys ; sys.exit(0)
    from scplacement import SCplacement
    import prepareChainV3
    scPrepC, useGivenRot, useDEE = prepareChainV3.PrepareChain("SCL1.0"), 1, 1
    badresids = ["VAL   85 ", "ASP   86 ", "TYR   68 ", "TYR   90 ",],
    SCplacement("premodel2.pdb", 0.5, "mmm.pdb", "dotfile", useDEE, "phased1.mtz2fofc.map", "FP", "FC", "PHIC", "2F1-F2", 0, 5, None, useGivenRot,
        badresids, scPrepC).run()
    import sys ; sys.exit(0)
    replaceWaters("model1.pdb", "rtk0.map")
