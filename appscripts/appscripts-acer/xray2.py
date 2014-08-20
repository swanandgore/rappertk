'''we assume here that ccp4 and cns environments are properly initilized'''

from procrun import proc_run_exitOnError as execCmd
import os, re, string

sgCCP4toCNS = {}
sgCCP4toCNS["P1"] = "P1"
sgCCP4toCNS["P21"] = "P2(1)"
sgCCP4toCNS["P43"] = "P4(3)"
sgCCP4toCNS["P21212"] = "P2(1)2(1)2"
sgCCP4toCNS["P212121"] = "P2(1)2(1)2(1)"
sgCCP4toCNS["P1211"] = "P2(1)"
sgCCP4toCNS["P43212"] = "P4(3)2(1)2"
sgCCP4toCNS["P41212"] = "P4(1)2(1)2"
sgCCP4toCNS["P4122"] = "P4(1)22"
sgCCP4toCNS["P6522"] = "P6(5)22"
sgCCP4toCNS["P3221"] = "P3(2)21"
sgCCP4toCNS["C2221"] = "C222(1)"
sgCCP4toCNS["P6422"] = "P6(4)22"
sgCCP4toCNS["C121"] = "C121"
sgCCP4toCNS["C2"] = "C2"

def mapman(cnsmap, ccp4map) :
    cmd, inputs = "~/downloads/rave_linux/lx_mapman ", ["read x", cnsmap, "xplor", "write x", ccp4map, "ccp4", "quit"]
    execCmd(cmd, inputs)

def phenix(pdbfile, mtzfile):
    cmd, inputs = "phenix.refine %s %s simulated_annealing=true" % (pdbfile, mtzfile), []
    execCmd(cmd, inputs)
    
## FIX ME need to add relative path ##
def run_reduce(pdbfile, outfile):
    cmd, inputs = "reduce -Quiet -build -DB /home/nick/RAPPERtk_new/newrtk/molprobity/lib %s > %s" % (pdbfile, outfile), []
    execCmd(cmd, inputs, [0, 256])

def run_molprobity(pdbfile, outfile):
    cmd, inputs = "multichart-rtk %s > %s" % (pdbfile, outfile), []
    execCmd(cmd, inputs)

def mtz2hkl(phasedmtz, cnshkl) :
    cmd, inputs = "mtz2various hklin %s hklout %s" % (phasedmtz, cnshkl), []
    inputs.append("LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag")
    #xoinputs.append("LABIN FP=FP FREE=FreeR_flag")
    inputs.append("OUTPUT CNS")
    inputs.append("END")
    execCmd(cmd, inputs)

def mtz2hkl_TEST(phasedmtz, cnshkl) :
    cmd, inputs = "mtz2various hklin %s hklout %s" % (phasedmtz, cnshkl), []
    inputs.append("LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag")
    #xoinputs.append("LABIN FP=FP FREE=FreeR_flag")
    inputs.append("OUTPUT CNS")
    inputs.append("END")
    execCmd(cmd, inputs)

def omitmap(phasedmtz, omapname) :
    cmd, inputs = "omit hklin %s mapout %s" % (phasedmtz, omapname), []
    inputs.append("scale 1.0 -1.0")
    inputs.append("labin -")
    inputs.append("    FP=FP PHI=PHIC FC=FC")
    inputs.append("end")
    execCmd(cmd, inputs)


def fft(phasedmtz, ccp4map, scale1=2., scale2=1.) :
    cmd, inputs = "fft hklin %s mapout %s" % (phasedmtz, ccp4map), []
    inputs.append("xyzlim asu")
    inputs.append("scale F1 %f" % scale1)
    inputs.append("scale F2 %f" % scale2)
    inputs.append("labin -")
    inputs.append("    F1=FP SIG1=SIGFP PHI=PHIC FREE=FreeR_flag F2=FC")
    inputs.append("end")
    execCmd(cmd, inputs)


def sfall(pdb, rfreemtz, phasedmtz, reso=None) :
    cmd, inputs = "sfall xyzin %s hklin %s hklout %s" % (pdb, rfreemtz, phasedmtz), []
    inputs.append("LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag")
    inputs.append("LABOUT -")
    inputs.append("    FC=FC PHIC=PHIC FP=FP")
    inputs.append("    MODE SFCALC -")
    inputs.append("    XYZIN -")
    inputs.append("    HKLIN")
    if reso != None : inputs.append("RESOLUTION %f 50" % reso)
    inputs.append("badd 0.0")
    inputs.append("vdwr 2.5")
    inputs.append("end")
    execCmd(cmd, inputs)

## refmac refinment - restrained refinement 10 cycles with no hydrogens NDF 18/01/2007 ##
## need to change to use refmac latest version !!!!!! ##
def refmac(hklin, hklout, pdbin, pdbout) :
    cmd, inputs = "/home/nick/Downloads/refmac_linintel hklin %s hklout %s xyzin %s xyzout %s" % (hklin, hklout, pdbin, pdbout), []
    inputs.append("bfac set 20")
    inputs.append("NCYC 40")
    inputs.append("make hydr no")
    inputs.append("end")
    execCmd(cmd, inputs)

def refmac_weights(hklin, hklout, pdbin, pdbout, weight) :
    cmd, inputs = "/home/nick/Downloads/refmac_linintel hklin %s hklout %s xyzin %s xyzout %s" % (hklin, hklout, pdbin, pdbout), []
    inputs.append("bfac set 20")
    inputs.append("NCYC 40")
    inputs.append("make hydr no")
    inputs.append("WEIG %s" % weight)
    inputs.append("end")
    execCmd(cmd, inputs)

def refmac1(hklin, hklout, pdbin, pdbout, reso, reftype, wa, breftype, ncyc) :
    cmd, inputs = "refmac5 hklin %s hklout %s xyzin %s xyzout %s" % (hklin, hklout, pdbin, pdbout), []
    inputs.append("MAKE_restraints HYDRogens No")
    inputs.append("MAKE CHECk 0")
    inputs.append("BFAC  1  5  6  6  7.5")
    inputs.append("LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag")
    inputs.append("LABO FC=FC PHIC=PHIC    FWT=2FOFCWT PHWT=PH2FOFCWT DELFWT=FOFCWT  PHDELWT=PHFOFCWT")
    inputs.append("REFI TYPE %s RESOLUTION  20 %f" % (reftype,reso) )
    inputs.append("REFI RESI MLKF")
    inputs.append("REFI BREF %s" % breftype)
    inputs.append("WEIGHT MATRIX %f" % wa)
    inputs.append("SCALe TYPE BULK")
    inputs.append("SCALe LSSCale FIXBulk 200.0")
    inputs.append("NCYC %d" % ncyc)
    inputs.append("MONI FEW")
    inputs.append("end")
    execCmd(cmd, inputs)

    
def uniqueify(basemtz, rfreemtz, percentage=None) :
    if percentage != None : cmd = "uniqueify -p %f %s %s" % (percentage, basemtz, rfreemtz)
    else : cmd = "uniqueify %s %s" % (basemtz, rfreemtz)
    execCmd(cmd, [])


def cif2mtz(strfac, outmtz, a, b, c, A, B, G, sg) :
    cmd, inputs = "cif2mtz hklin %s hklout %s" % (strfac,outmtz), []
    inputs.append( "CELL %f %f %f %f %f %f" % (a,b,c,A,B,G)  )
    inputs.append( "SYMMETRY %s" % sg )
    execCmd(cmd, inputs)
    

def pdbSF2phasedMTZ(pdb, strfac, basemtz, rfreemtz, phasedmtz, a, b, c, A, B, G, sg) :
    '''all arguments are filenames. first 2 are given, rest are generated in this function.'''
    # cif2mtz
    cif2mtz(strfac, basemtz, a, b, c, A, B, G, sg)

def findSegids(pdbfn) :
    segids = set() ; from pdbr import line2segid
    for l in open(pdbfn,'r').readlines() :
        if l[0:6] in ["HETATM","ATOM  "] and not line2segid(l) in [None,"    "] :
            segids.add( line2segid(l) )
    return segids

def makeIgroupLine(segids, atomselect) :
    lines = []
    intline = "igroup "
    for segi in segids :
        if len(intline) > 200 : lines.append( intline ) ; intline = ""
        intline += 'interaction (%s segid "%s") (%s segid "%s") ' % (atomselect,segi,atomselect,segi)
    retline = ""
    for l in lines : retline += l + " \n"
    retline += intline
    if retline[ len(retline)-1 ] == "\n" : retline = retlines[0:len(retline)-1]
    print len(retline), retline
    return retline

def cns_generate(pdbfn, mtfout, pdbout, topCNS=None, parCNS=None, logfile=None) :
    generate_inp = "%s/appscripts/generate.inp" % os.environ["RTKROOT"]
    inputs = []
    segids = findSegids(pdbfn)
    for l in open(generate_inp, 'r').readlines() :
        if re.compile("prot_coordinate_infile_1.*amy.pdb").search(l) : inputs.append( '{===>} prot_coordinate_infile_1="%s";' % pdbfn )
        elif re.compile("structure_outfile=").search(l) : inputs.append( '{===>} structure_outfile="%s";' % mtfout )
        elif re.compile("coordinate_outfile=").search(l) : inputs.append( '{===>} coordinate_outfile="%s";' % pdbout )
        elif topCNS and re.compile("lig_topology_infile=").search(l) : inputs.append( '{===>} lig_topology_infile="%s";' % topCNS )
        elif parCNS and re.compile("lig_parameter_infile=").search(l) : inputs.append( '{===>} lig_parameter_infile="%s";' % parCNS )
        elif re.compile("igroup").search(l) and len(segids) > 0 : inputs.append( makeIgroupLine(segids,"") )
        else : inputs.append(re.sub("\n", "", l))
    fp = open('CNSgenerate.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close() ;
    if logfile == None : logfile = "generate.log"
    execCmd("cns > "+logfile, inputs)

def cns_cis_peptide(logfile=None) :
    generate_inp = "%s/appscripts/cis_peptide.inp" % os.environ["RTKROOT"]
    inputs = []
    for l in open(generate_inp, 'r').readlines() :
        inputs.append(re.sub("\n", "", l))
    fp = open('CNScis_peptide.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close()
    if logfile == None : logfile = "cis_peptide.log"
    execCmd("cns > "+logfile, inputs)
    
def cns_generate_multi_chain(pdbfn, mtfout, pdbout, topCNS=None, parCNS=None, logfile=None) :
    generate_inp = "%s/appscripts/generate_easy.inp" % os.environ["RTKROOT"]
    inputs = []
    for l in open(generate_inp, 'r').readlines() :
        if re.compile("prot_coordinate_infile_1.*amy.pdb").search(l) : inputs.append( '{===>} prot_coordinate_infile_1="%s";' % pdbfn )
        elif re.compile("structure_outfile=").search(l) : inputs.append( '{===>} structure_outfile="%s";' % mtfout )
        elif re.compile("coordinate_outfile=").search(l) : inputs.append( '{===>} coordinate_outfile="%s";' % pdbout )
        elif topCNS and re.compile("lig_topology_infile=").search(l) : inputs.append( '{===>} lig_topology_infile="%s";' % topCNS )
        elif parCNS and re.compile("lig_parameter_infile=").search(l) : inputs.append( '{===>} lig_parameter_infile="%s";' % parCNS )
        else : inputs.append(re.sub("\n", "", l))
    fp = open('CNSgenerate.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close()
    if logfile == None : logfile = "generate.log"
    execCmd("cns > "+logfile, inputs)


## changed to use old school cns_refine routine NDF 17/01/2007 ##
def cns_anneal(a,b,c,A,B,G,sg,reso, cnshkl, mtffn, pdbfn, parCNS=None, logfile=None, wa=None, num_cycles=2, temp=5000, harmCA=None) :
    print "CNS Annealing for %d cycles starting at %dK" % (num_cycles, temp)
    print sg, reso
    anneal_inp = "%s/appscripts/anneal.inp" % os.environ["RTKROOT"]
    inputs = []
    segids = findSegids(pdbfn)
    for l in open(anneal_inp, 'r').readlines() :
        if re.compile("{===>} a=").search(l) : inputs.append( '{===>} a=%f;' % a )
        elif re.compile("{===>} b=").search(l) : inputs.append( '{===>} b=%f;' % b )
        elif re.compile("^{===>} temperature=").search(l) : inputs.append( '{===>} temperature=%d;' % temp )
        elif re.compile("^{===>} num_cycles=").search(l) : inputs.append( '{===>} num_cycles=%d;' % num_cycles )
        elif re.compile("{===>} c=").search(l) : inputs.append( '{===>} c=%f;' % c )
        elif re.compile("{===>} alpha=").search(l) : inputs.append( '{===>} alpha=%f;' % A )
        elif re.compile("{===>} beta=").search(l) : inputs.append( '{===>} beta=%f;' % B )
        elif re.compile("{===>} gamma=").search(l) : inputs.append( '{===>} gamma=%f;' % G )
        elif re.compile("{===>} sg=").search(l) : inputs.append( '{===>} sg="%s";' % sg )
        elif re.compile("{===>} high_res=").search(l) : inputs.append( '{===>} high_res=%f;' % reso )
        elif wa != None and re.compile("{===>} wa=").search(l) : inputs.append( '{===>} wa=%f;' % wa )
        elif re.compile("{===>} reflection_infile_1=").search(l) : inputs.append( '{===>} reflection_infile_1="%s";' % cnshkl )
        elif re.compile("{===>} structure_infile=").search(l) : inputs.append( '{===>} structure_infile="%s";' % mtffn )
        elif re.compile("{===>} coordinate_infile=").search(l) : inputs.append( '{===>} coordinate_infile="%s";' % pdbfn )
        elif re.compile("{===>} atom_harm=").search(l) and harmCA : inputs.append('{===>} atom_harm=(name ca);')
        elif re.compile("parameter_infile_extra=").search(l) :
            if parCNS : inputs.append( '{===>} parameter_infile_extra="%s";' % parCNS )
        elif re.compile("igroup").search(l) and len(segids) > 0 : inputs.append( makeIgroupLine(segids," &atom_select and ") )
        else : inputs.append(re.sub("\n", "", l))
        #elif re.compile("atom_harm=").search(l) : inputs.append("atom_harm = ( name ca );")
    fp = open('CNSanneal.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close()
    if logfile == None : logfile = "anneal.log"
    execCmd("cns > "+logfile, inputs) 




def cns_bindividual(a,b,c,A,B,G,sg,reso, cnshkl, mtffn, pdbfn, parCNS=None, logfile=None, wa=None, num_cycles=2, temp=5000, harmCA=None) :

    print sg, reso
    bindividual_inp = "%s/appscripts/bindividual.inp" % os.environ["RTKROOT"]
    inputs = []
    segids = findSegids(pdbfn)
    for l in open(bindividual_inp, 'r').readlines() :

        if re.compile("{===>} a=").search(l) : inputs.append( '{===>} a=%f;' % a )
        elif re.compile("{===>} b=").search(l) : inputs.append( '{===>} b=%f;' % b )
        elif re.compile("{===>} c=").search(l) : inputs.append( '{===>} c=%f;' % c )

        elif re.compile("{===>} alpha=").search(l) : inputs.append( '{===>} alpha=%f;' % A )
        elif re.compile("{===>} beta=").search(l) : inputs.append( '{===>} beta=%f;' % B )
        elif re.compile("{===>} gamma=").search(l) : inputs.append( '{===>} gamma=%f;' % G )
        elif re.compile("{===>} sg=").search(l) : inputs.append( '{===>} sg="%s";' % sg )
        elif re.compile("{===>} high_res=").search(l) : inputs.append( '{===>} high_res=%f;' % reso )

        elif re.compile("{===>} reflection_infile_1=").search(l) : inputs.append( '{===>} reflection_infile_1="%s";' % cnshkl )
        elif re.compile("{===>} structure_infile=").search(l) : inputs.append( '{===>} structure_infile="%s";' % mtffn )
        elif re.compile("{===>} coordinate_infile=").search(l) : inputs.append( '{===>} coordinate_infile="%s";' % pdbfn )

        elif re.compile("igroup").search(l) and len(segids) > 0 : inputs.append( makeIgroupLine(segids," &atom_select and ") )
        else : inputs.append(re.sub("\n", "", l))

    fp = open('CNSbindividual.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close()
    if logfile == None : logfile = "bindividual.log"

    execCmd("cns > "+logfile, inputs) 



def cns_minimize(a,b,c,A,B,G,sg, mtffn, pdbfn, parCNS=None, logfile=None):

    minimize_inp = "%s/appscripts/model_minimize.inp" % os.environ["RTKROOT"]
    inputs = []
    segids = findSegids(pdbfn)
    for l in open(minimize_inp, 'r').readlines() :

        if re.compile("{===>} a=").search(l) : inputs.append( '{===>} a=%f;' % a )
        elif re.compile("{===>} b=").search(l) : inputs.append( '{===>} b=%f;' % b )
        elif re.compile("{===>} c=").search(l) : inputs.append( '{===>} c=%f;' % c )

        elif re.compile("{===>} alpha=").search(l) : inputs.append( '{===>} alpha=%f;' % A )
        elif re.compile("{===>} beta=").search(l) : inputs.append( '{===>} beta=%f;' % B )
        elif re.compile("{===>} gamma=").search(l) : inputs.append( '{===>} gamma=%f;' % G )
        elif re.compile("{===>} sg=").search(l) : inputs.append( '{===>} sg="%s";' % sg )

        elif re.compile("{===>} structure_infile=").search(l) : inputs.append( '{===>} structure_infile="%s";' % mtffn )
        elif re.compile("{===>} coordinate_infile=").search(l) : inputs.append( '{===>} coordinate_infile="%s";' % pdbfn )
        else : inputs.append(re.sub("\n", "", l))

    fp = open('CNSminimize.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close()
    if logfile == None : logfile = "minimize.log"

    execCmd("cns > "+logfile, inputs) 




def findResolutionReflections(mtzfile, numref) :
    for reso in range(100, 0, -1) :
        inputs = []
        inputs.append("stats nbin 1 reso 1000 %f" % (reso/10.) )
        inputs.append("go")
        es, ol, el = execCmd("mtzdump hklin %s | grep 'No. of reflections used in FILE STATISTICS'" % mtzfile, inputs)
        nr = string.atoi(re.sub("No. of reflections used in FILE STATISTICS", "", ol[0]))
        print "resoNR %5.1f %6d" % (reso/10., nr)
        if numref < nr : return reso/10.
    return None

if __name__ == "__main__" :
    print findResolutionReflections("rfree.mtz", 21000)
    sys.exit(0)
    omitmap("phased1.mtz", "omit1.map")
    sys.exit(0)
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdb', help='phase estimates in form of PDB file')
    parser.add_option("--sf", action='store', type='string', dest='sf', help='structure factors file')
    parser.add_option("--basemtz", action='store', type='string', dest='basemtz', help='base mtz filename')
    parser.add_option("--rfreemtz", action='store', type='string', dest='rfreemtz', help='rfree mtz filename')
    parser.add_option("--phasedmtz", action='store', type='string', dest='phasedmtz', help='phased mtz filename')
    parser.add_option("--ccp4map", action='store', type='string', dest='ccp4map', help='CCP4 map filename')
    parser.add_option("--cnshkl", action='store', type='string', dest='cnshkl', help='CNS hkl filename')
    parser.add_option("--ccp4ORcns", action='store', type='string', dest='ccp4ORcns', help='run in ccp4 or cns mode', default="ccp4")
    parser.add_option("--a", action='store', type='float', dest='a', help='cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='A', help='cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='B', help='cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='G', help='cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='cell spacegroup, in CCP4 notation')
    parser.add_option("--topCNS", action='store', type='string', dest='topCNS', help='additional topology file for CNS', default=None)
    parser.add_option("--parCNS", action='store', type='string', dest='parCNS', help='additional parameter file for CNS', default=None)
    (options, args) = parser.parse_args()

    if options.ccp4ORcns == "ccp4" :
        cif2mtz(options.sf, options.basemtz, options.a, options.b, options.c, options.A, options.B, options.G, options.sg)
        uniqueify(options.basemtz, options.rfreemtz)
        sfall(options.pdb, options.rfreemtz, options.phasedmtz)
        fft(options.phasedmtz, options.ccp4map)
    elif options.ccp4ORcns == "cns" :
        mtz2hkl(options.phasedmtz, options.cnshkl)
        cns_generate(options.pdb, "generate.mtf", "generate.pdb", options.topCNS, options.parCNS)
        cns_anneal(options.a, options.b, options.c, options.A, options.B, options.G, sgCCP4toCNS[options.sg], options.cnshkl, "generate.mtf", "generate.pdb", options.parCNS)
    else : print "cnshkl options shd be ccp4 or cns" ; assert None
