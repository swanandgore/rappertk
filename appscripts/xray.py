'''we assume here that ccp4 and cns environments are properly initilized'''

from procrun import proc_run_exitOnError as execCmd
import os, re

sgCCP4toCNS = {}
sgCCP4toCNS["P21"] = "P2(1)"
sgCCP4toCNS["P212121"] = "P2(1)2(1)2(1)"


def mtz2hkl(phasedmtz, cnshkl) :
    cmd, inputs = "mtz2various hklin %s hklout %s" % (phasedmtz, cnshkl), []
    inputs.append("LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag")
    inputs.append("OUTPUT CNS")
    inputs.append("END")
    execCmd(cmd, inputs)


def fft(phasedmtz, ccp4map) :
    cmd, inputs = "fft hklin %s mapout %s" % (phasedmtz, ccp4map), []
    inputs.append("xyzlim asu")
    inputs.append("scale F1 2.0")
    inputs.append("scale F2 1.0")
    inputs.append("labin -")
    inputs.append("    F1=FP SIG1=SIGFP PHI=PHIC FREE=FreeR_flag F2=FC")
    inputs.append("end")
    execCmd(cmd, inputs)


def sfall(pdb, rfreemtz, phasedmtz) :
    cmd, inputs = "sfall xyzin %s hklin %s hklout %s" % (pdb, rfreemtz, phasedmtz), []
    inputs.append("LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag")
    inputs.append("LABOUT -")
    inputs.append("    FC=FC PHIC=PHIC")
    inputs.append("    MODE SFCALC -")
    inputs.append("    XYZIN -")
    inputs.append("    HKLIN")
    inputs.append("badd 0.0")
    inputs.append("vdwr 2.5")
    inputs.append("end")
    execCmd(cmd, inputs)


def uniqueify(basemtz, rfreemtz) :
    cmd, inputs = "uniqueify %s %s" % (basemtz, rfreemtz), []
    execCmd(cmd, inputs)


def cif2mtz(strfac, outmtz, a, b, c, A, B, G, sg) :
    cmd, inputs = "cif2mtz hklin %s hklout %s" % (strfac,outmtz), []
    inputs.append( "CELL %f %f %f %f %f %f" % (a,b,c,A,B,G)  )
    inputs.append( "SYMMETRY %s" % sg )
    execCmd(cmd, inputs)
    

def pdbSF2phasedMTZ(pdb, strfac, basemtz, rfreemtz, phasedmtz, a, b, c, A, B, G, sg) :
    '''all arguments are filenames. first 2 are given, rest are generated in this function.'''
    # cif2mtz
    cif2mtz(strfac, basemtz, a, b, c, A, B, G, sg)


def cns_generate(pdbfn, mtfout, pdbout, topCNS=None, parCNS=None, logfile=None) :
    generate_inp = "%s/appscripts/generate.inp" % os.environ["RTKROOT"]
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


def cns_anneal(a,b,c,A,B,G,sg, cnshkl, mtffn, pdbfn, parCNS=None, logfile=None) :
    print sg
    anneal_inp = "%s/appscripts/anneal.inp" % os.environ["RTKROOT"]
    inputs = []
    for l in open(anneal_inp, 'r').readlines() :
        if re.compile("{===>} a=").search(l) : inputs.append( '{===>} a=%f;' % a )
        elif re.compile("{===>} b=").search(l) : inputs.append( '{===>} b=%f;' % b )
        elif re.compile("{===>} c=").search(l) : inputs.append( '{===>} c=%f;' % c )
        elif re.compile("{===>} alpha=").search(l) : inputs.append( '{===>} alpha=%f;' % A )
        elif re.compile("{===>} beta=").search(l) : inputs.append( '{===>} beta=%f;' % B )
        elif re.compile("{===>} gamma=").search(l) : inputs.append( '{===>} gamma=%f;' % G )
        elif re.compile("{===>} sg=").search(l) : inputs.append( '{===>} sg="%s";' % sg )
        elif re.compile("{===>} reflection_infile_1=").search(l) : inputs.append( '{===>} reflection_infile_1="%s";' % cnshkl )
        elif re.compile("{===>} structure_infile=").search(l) : inputs.append( '{===>} structure_infile="%s";' % mtffn )
        elif re.compile("{===>} coordinate_infile=").search(l) : inputs.append( '{===>} coordinate_infile="%s";' % pdbfn )
        elif parCNS and re.compile("parameter_infile_extra=").search(l) : inputs.append( '{===>} parameter_infile_extra="%s";' % parCNS )
        else : inputs.append(re.sub("\n", "", l))
    fp = open('CNSanneal.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close()
    if logfile == None : logfile = "anneal.log"
    execCmd("cns > "+logfile, inputs) 


if __name__ == "__main__" :
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
