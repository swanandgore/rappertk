'''we assume here that ccp4 and cns environments are properly initilized'''

from procrun import proc_run_exitOnError as execCmd
from procrun import proc_run as execCmd2
import os, re, string


sgCCP4toCNS = {}
sgCCP4toCNS["P1"] = "P1"
sgCCP4toCNS["P6"] = "P6"
sgCCP4toCNS["P63"] = "P6(3)"
sgCCP4toCNS["P65"] = "P6(5)"
sgCCP4toCNS["P21"] = "P2(1)"
sgCCP4toCNS["P43"] = "P4(3)"
sgCCP4toCNS["P422"] = "P422"
sgCCP4toCNS["P21212"] = "P2(1)2(1)2"
sgCCP4toCNS["P212121"] = "P2(1)2(1)2(1)"
sgCCP4toCNS["P1211"] = "P2(1)"
sgCCP4toCNS["P43212"] = "P4(3)2(1)2"
sgCCP4toCNS["P41212"] = "P4(1)2(1)2"
sgCCP4toCNS["P6522"] = "P6(5)22"
sgCCP4toCNS["P3221"] = "P3(2)21"
sgCCP4toCNS["P3121"] = "P3(1)21"
sgCCP4toCNS["P4132"] = "P4(1)32"

sgCCP4toCNS["P4122"] = "P4(1)22"
sgCCP4toCNS["C2221"] = "C222(1)"
sgCCP4toCNS["P6422"] = "P6(4)22"

sgCCP4toCNS["C2"] = "C2"
sgCCP4toCNS["I4122"] = "I4(1)22"
sgCCP4toCNS["I213"] = "I213"
## new

sgCCP4toCNS["P6322"] = "P6(3)22"
sgCCP4toCNS["P622"] = "P622"
sgCCP4toCNS["P61"] = "P6(1)"
sgCCP4toCNS["P321"] = "P321"

#
#13052011:Takashi:Added CAD
#
def cadnat(mtzin, mtzout, Fobs, SIGFobs) :
    cmd, inputs = "cad hklin1 %s hklout %s" % (mtzin, mtzout), []
    inputs.append("LABIN FILE_NUMBER 1 E1=%s E2=%s" % (Fobs, SIGFobs))
    inputs.append("LABO E1=FP E2=SIGFP")
    inputs.append("CTYPEIN FILE_NUMBER 1 E1=F E2=Q")
    execCmd(cmd, inputs)

def cad(mtzin, mtzout, Fobs, SIGFobs, Fcalc, PHIcalc, FOout, SIGOout) :
    cmd, inputs = "cad hklin1 %s hklout %s" % (mtzin, mtzout), []
    inputs.append("LABIN FILE_NUMBER 1 E1=%s E2=%s E3=%s E4=%s" % (Fobs, SIGFobs, Fcalc, PHIcalc))
    inputs.append("LABO E1=%s E2=%s E3=FC E4=PHIC" % (FOout, SIGOout))
    inputs.append("CTYPEIN FILE_NUMBER 1 E1=F E2=Q E3=F E4=P")
    execCmd(cmd, inputs)

def mapman(cnsmap, ccp4map) :
    cmd, inputs = "/home/ak459/downloads/rave_linux/lx_mapman ", ["read x", cnsmap, "xplor", "write x", ccp4map, "ccp4", "quit"]
    execCmd(cmd, inputs)

# 12052011: Takashi: Added parameter file as an input file
def phenix(pdbfile, mtzfile, param):
    cmd, inputs = "phenix.refine %s %s %s" % (pdbfile, mtzfile, param), []
    execCmd(cmd, inputs)
    
## FIX ME need to add relative path ##
def run_reduce(pdbfile, outfile):
    cmd, inputs = " /home/ak459/downloads/molprobity3/bin/linux/reduce -Quiet -build  %s > %s" % (pdbfile, outfile), []
    execCmd(cmd, inputs, [0, 256])

def run_molprobity(pdbfile, outfile):
    cmd, inputs = "~/rtk152-O/molprobity/cmdline/multichart-rtk  %s > %s" % (pdbfile, outfile), []
    execCmd(cmd, inputs)

def mtz2hkl(phasedmtz, cnshkl) :
    cmd, inputs = "mtz2various hklin %s hklout %s" % (phasedmtz, cnshkl), []
    inputs.append("LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag")
    inputs.append("OUTPUT CNS")
    inputs.append("END")
    execCmd(cmd, inputs)


def mtz2hkl2(phasedmtz, cnshkl,fp,sigfp,freer) :
    cmd, inputs = "mtz2various hklin %s hklout %s" % (phasedmtz, cnshkl), []
    inputs.append("LABIN FP=%s SIGFP=%s FREE=%s"%(fp,sigfp,freer))
    inputs.append("OUTPUT CNS")
    inputs.append("END")
    execCmd(cmd, inputs)    

def omitmap(phasedmtz, omapname,fplabel,fclabel,phiclabel,siglabel,n,m,freeRlabel=None) :
    cmd, inputs = "omit hklin %s mapout %s" % (phasedmtz, omapname), []
    aline = "scale %d -%d"%(n,m)
    inputs.append(aline)
    inputs.append("labin -")
    if freeRlabel == None :
        iline = "    FP=%s PHI=%s FC=%s "%(fplabel,phiclabel,fclabel)
        inputs.append(iline)
    else :
        iline = "    FP=%s PHI=%s FC=%s "%(fplabel,phiclabel,fclabel)
        inputs.append(iline)
        

    inputs.append("end")


    fp = open('omit.sh', 'w')
    print >> fp, cmd
    for l in inputs : print >> fp, l

    fp.close() ;
    execCmd2(cmd, inputs)
    import os, sys
    if  not os.path.isfile(omapname) :
        return None
    return 1 

def mapdump(mtzfile):

    inputs = []
    #inputs.append("stats nbin 1 reso 1000 %f" % (reso/10.) )
    inputs.append("go")
    es, ol, el = execCmd2("mapdump mapin %s " % mtzfile, inputs)
    if es !=0 :
        return None
    return 1
        
def mtzdump(mtzfile,fplabel,sigfplabel,freerlabel=None,fclabel=None,phiclabel=None):

    inputs = []
    inputs.append("go")
    es, ol, el = execCmd2("mtzdump hklin %s " % mtzfile, inputs)


    if es !=0 :
        return None
    lc = 0 ;     checkline = 0 ;    typeline = 0
    
    for oline in ol :
        if '  * Column Labels :' in oline:
            checkline = lc
        if '  * Column Types :' in oline:
            typeline = lc
        else :
            lc = lc + 1
    if typeline == 0 : typeline = lc
    fp_found = 0 ; sigfp_found = 0 ; freer_found = 0 ; fc_found = 0 ; phic_found = 0
    
    for xline in range(checkline,typeline):
        if  re.compile(' H K L').search(ol[xline]) == None :
            continue
        words = ol[xline].split()
        for abc in words :

            if abc == sigfplabel  :
                sigfp_found = 1
            if abc == fplabel  :
                fp_found = 1
                
            if freerlabel != None:
                if  abc == freerlabel :
                    freer_found = 1
                
            if fclabel != None:
                if  abc == fclabel : 
                    fc_found  = 1
 
            if phiclabel != None:
                if  abc == phiclabel:
                    phic_found = 1
    if  fp_found == 0 :
        print "%s not found in mtz file %s"%(fplabel,mtzfile)
        print "exiting.."
        return None
    if  sigfp_found == 0 :
        print "%s not found in mtz file %s"%(sigfplabel,mtzfile)
        print "exiting.."
        return None

    if freerlabel != None:

        if  freer_found == 0:
            print "%s not found in mtz file %s"%(freerlabel,mtzfile)
            print "exiting.."
            return None
    if fclabel != None:
        if  fc_found == 0 :
            print "%s not found in mtz file %s"%(fclabel,mtzfile)
            print "exiting.."
            return None        
    if phiclabel != None:
        if  phic_found == 0:
            print "%s not found in mtz file %s"%(phiclabel,mtzfile)
            print "exiting.."
            return None
    return 1


    

def fft(phasedmtz, ccp4map, fplabel,fclabel,phiclabel,sigfplabel,scale1=2., scale2=1.,freerlabel=None) :
    cmd, inputs = "fft hklin %s mapout %s" % (phasedmtz, ccp4map), []
    inputs.append("xyzlim asu")
    inputs.append("scale F1 %f" % scale1)
    inputs.append("scale F2 %f" % scale2)
    inputs.append("labin -")
    if freerlabel == None :
        inputs.append("    F1=%s SIG1=%s PHI=%s F2=%s"%(fplabel,sigfplabel,phiclabel,fclabel))
    else :
        inputs.append("    F1=%s SIG1=%s PHI=%s FREE=%s F2=%s"%(fplabel,sigfplabel,phiclabel,freerlabel,fclabel))
    inputs.append("end")
    execCmd(cmd, inputs)


def rstats(phasedmtz, scaledmtz, fplabel,sigfplabel, fclabel,phiclabel,freerlabel):
    cmd, inputs = "rstats hklin %s hklout %s" % (phasedmtz, scaledmtz), []
    inputs.append("labin -")
    if freerlabel != None :
        inputs.append("    FP=%s SIGFP=%s PHIC=%s FC=%s FREE=%s"%(fplabel,sigfplabel,phiclabel,fclabel,freerlabel))
    else :
        inputs.append("    FP=%s SIGFP=%s PHIC=%s FC=%s "%(fplabel,sigfplabel,phiclabel,fclabel))
    inputs.append("end")

    es, ol, el = execCmd2(cmd, inputs)
    if es !=0 :
        return None
    return 1





def fftnew(phasedmtz, ccp4map, pdbfile,fplabel,fclabel,phiclabel,sigfplabel,scale1=2., scale2=1.,freerlabel=None,onlyFC=0) :
    #rstats_status = rstats(phasedmtz,"scaled.mtz",fplabel,sigfplabel,fclabel,phiclabel,freerlabel)

    #if rstats_status == None :
    #    print "Cannot generate difference maps..exiting"
    #    import sys ; sys.exit()

    imap =   phasedmtz # "scaled.mtz" # 
    map = "pre."+ccp4map
    cmd, inputs = "fft hklin %s mapout %s" % (imap, map), []
    inputs.append("xyzlim asu")
    if onlyFC == 1:
        inputs.append("labin -")
        if freerlabel == None :
            inputs.append("    F1=%s  PHI=%s "%(fclabel,phiclabel))
        else :
            inputs.append("    F1=%s  PHI=%s FREE=%s "%(fclabel,phiclabel,freerlabel))
    else :
        inputs.append("scale F1 %f" % scale1)
        inputs.append("scale F2 %f" % scale2)
        inputs.append("labin -")
        if freerlabel == None :
            inputs.append("    F1=%s SIG1=%s PHI=%s F2=%s"%(fplabel,sigfplabel,phiclabel,fclabel))
        else :
            inputs.append("    F1=%s SIG1=%s PHI=%s FREE=%s F2=%s"%(fplabel,sigfplabel,phiclabel,freerlabel,fclabel))

    inputs.append("end")
    es, ol, el = execCmd2(cmd, inputs)
    if es !=0 :
        print "Cannot generate difference maps..exiting"
        import sys ; sys.exit()
    mapmask(map,ccp4map,pdbfile)


def mapmask(mapin,mapout,pdbfile):
    cmd, inputs = "mapmask mapin %s mapout %s xyzin %s" % ( mapin,mapout,pdbfile), []
    inputs.append("border 4")
    inputs.append("mode mapin")
    inputs.append("end")
    #execCmd(cmd, inputs)
    es, ol, el = execCmd2(cmd, inputs)
    if es !=0 :
        print "Mapmask failed ..Cannot generate difference maps..exiting"
        import sys ; sys.exit()


def sfall2(pdb, rfreemtz, phasedmtz,FREE ="FreeR_flag" , reso=None) :
    cmd, inputs = "sfall xyzin %s hklin %s hklout %s" % (pdb, rfreemtz, phasedmtz), []
    dd = "LABIN FP=FP SIGFP=SIGFP FREE=%s"%FREE
    inputs.append(dd)
    inputs.append("LABOUT -")
    inputs.append("    FC=FC PHIC=PHIC")
    inputs.append("    MODE SFCALC -")
    inputs.append("    XYZIN -")
    inputs.append("    HKLIN")
    if reso != None : inputs.append("RESOLUTION %f 50" % reso)
    inputs.append("badd 0.0")
    inputs.append("vdwr 2.5")
    inputs.append("end")
    execCmd(cmd, inputs)
#
#13052011:Takashi:Removed the default labels for amplitude and sigma of structure factors
#
def sfall(pdb, rfreemtz, phasedmtz, fplabel, sigfplabel, freerlabel = "FreeR_flag",reso=None, usefreer=1) :
    
    cmd, inputs = "sfall xyzin %s hklin %s hklout %s" % (pdb, rfreemtz, phasedmtz), []

    if usefreer == 0 : 
        inputs.append("LABIN FP=%s SIGFP=%s "%(fplabel,sigfplabel))

    else :
        inputs.append("LABIN FP=%s SIGFP=%s FREE=%s"%(fplabel,sigfplabel,freerlabel))


    inputs.append("LABOUT -")
    inputs.append("    FC=FC PHIC=PHIC")
    inputs.append("    MODE SFCALC -")
    inputs.append("    XYZIN -")
    inputs.append("    HKLIN")
    if reso != None : inputs.append("RESOLUTION %f 50" % reso)
    inputs.append("badd 0.0")
    inputs.append("vdwr 2.5")
    inputs.append("end")

    es, ol, el = execCmd2(cmd, inputs)
    if es !=0 :
        return None
    return 1


## refmac refinment - restrained refinement 10 cycles with no hydrogens NDF 18/01/2007 ##
## need to change to use refmac latest version !!!!!! ##
def refmac(hklin, hklout, pdbin, pdbout) :
    cmd, inputs = "~/downloads/CCP4/ccp4-6.0.2/bin/refmac5 hklin %s hklout %s xyzin %s xyzout %s" % (hklin, hklout, pdbin, pdbout), []
    inputs.append("bfac set 20")
    inputs.append("NCYC 100")
    inputs.append("make hydr no")
    inputs.append("end")
    execCmd(cmd, inputs)

def sfcheck(hklin,pdbin,fplabel,sigfplabel,freerlabel) :
    cmd, inputs = "sfcheck -m %s -f  %s " % (pdbin, hklin), []
    if fplabel == None or sigfplabel  == None :
        dd  = " "
    elif freerlabel != None : 
        dd = "          -lf %s  -lsf %s -lfree %s "%(fplabel,sigfplabel,freerlabel)
    else :
        dd = "          -lf %s  -lsf %s "%(fplabel,sigfplabel)
    dd = " "
    inputs.append(dd)


    execCmd(cmd, inputs)

    

def refmacNew(hklin, hklout, pdbin, pdbout,fp="FP",sigfp="SIGFP") :
    cmd, inputs = "~/downloads/CCP4/ccp4-6.0.2/bin/refmac5 hklin %s hklout %s xyzin %s xyzout %s" % (hklin, hklout, pdbin, pdbout), []
    inputs.append("weights AUTO")
    inputs.append("LABIN FP=%s SIGFP=%s FREE=FreeR_flag"%(fp,sigfp))
    inputs.append("make check  NONE")
    inputs.append("refi type REST PHASE SCBL 1.0 BBLU 0.0 resi MLKF meth CGMAT bref ISOT")
    inputs.append("ncyc 40")
    inputs.append("scal type SIMP LSSC ANISO EXPE")
    inputs.append("solvent YES VDWProb 1.4 IONProb 0.8 RSHRink 0.8")
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
    generate_inp = "%s/appscripts-acer/generate.inp" % os.environ["RTKROOT"]
    inputs = []
    segids = findSegids(pdbfn)
    for l in open(generate_inp, 'r').readlines() :
        if re.compile("prot_coordinate_infile_1.*amy.pdb").search(l) : inputs.append( '{===>} prot_coordinate_infile_1="%s";' % pdbfn )
        elif re.compile("structure_outfile=").search(l) : inputs.append( '{===>} structure_outfile="%s";' % mtfout )
        elif re.compile("coordinate_outfile=").search(l) : inputs.append( '{===>} coordinate_outfile="%s";' % pdbout )
        elif topCNS and re.compile("lig_topology_infile=").search(l) : inputs.append( '{===>} lig_topology_infile="%s";' % topCNS )
        elif parCNS and re.compile("lig_parameter_infile=").search(l) : inputs.append( '{===>} lig_parameter_infile="%s";' % parCNS )
        #elif re.compile("igroup").search(l) and len(segids) > 0 : inputs.append( makeIgroupLine(segids,"") )
        else : inputs.append(re.sub("\n", "", l))
    fp = open('CNSgenerate.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close() ;
    if logfile == None : logfile = "generate.log"
    #execCmd("/home/anjum/downloads/cns_solve_1.1/intel-i686-linux_g77/bin/cns_solve > "+logfile, inputs)
    execCmd("/home/ak459//downloads/CNS/cns_solve_1.2/intel-x86_64bit-linux/bin/cns_solve > "+logfile, inputs)

def cns_cis_peptide(logfile=None) :
    generate_inp = "%s/appscripts-acer/cis_peptide.inp" % os.environ["RTKROOT"]
    inputs = []
    for l in open(generate_inp, 'r').readlines() :
        inputs.append(re.sub("\n", "", l))
    fp = open('CNScis_peptide.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close()
    if logfile == None : logfile = "cis_peptide.log"
    #execCmd("/home/anjum/downloads/cns_solve_1.1/intel-i686-linux_g77/bin/cns_solve > "+logfile, inputs)
    execCmd("/home/ak459//downloads/CNS/cns_solve_1.2/intel-x86_64bit-linux/bin/cns_solve  > "+logfile, inputs)
    
def cns_generate_multi_chain(pdbfn, mtfout, pdbout, topCNS=None, parCNS=None, logfile=None) :
    generate_inp = "%s/appscripts-acer/generate_easy.inp" % os.environ["RTKROOT"]
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
    execCmd("/home/ak459//downloads/CNS/cns_solve_1.2/intel-x86_64bit-linux/bin/cns_solve  > "+logfile, inputs)


## changed to use old school cns_refine routine NDF 17/01/2007 ##
def cns_anneal(a,b,c,A,B,G,sg,reso, cnshkl, mtffn, pdbfn, parCNS=None, logfile=None, wa=None, num_cycles=2, temp=5000, harmCA=None) :
    import os
    print "CNS Annealing for %d cycles starting at %dK" % (num_cycles, temp)
    print sg, reso
    anneal_inp = "%s/appscripts-acer/anneal.inp" % os.environ["RTKROOT"]
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
        #elif re.compile("igroup").search(l) and len(segids) > 0 : inputs.append( makeIgroupLine(segids," &atom_select and ") )
        else : inputs.append(re.sub("\n", "", l))
        #elif re.compile("atom_harm=").search(l) : inputs.append("atom_harm = ( name ca );")
    fp = open('CNSanneal.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close()
    if logfile == None : logfile = "anneal.log"
    #execCmd("/home/anjum/downloads/cns_solve_1.1/intel-i686-linux_g77/bin/cns_solve > "+logfile, inputs)
    import os # os.system
    #import sys ; sys.exit()
    #os.system("/home/ak459//downloads/CNS/cns_solve_1.2/intel-x86_64bit-linux/bin/cns_solve  > CNSanneal.inp %s"%logfile)
    execCmd("/home/ak459//downloads/CNS/cns_solve_1.2/intel-x86_64bit-linux/bin/cns_solve  > "+logfile, inputs)
    #import sys ; sys.exit()



def cns_anneal2(a,b,c,A,B,G,sg,reso, cnshkl, mtffn, pdbfn, parCNS=None, logfile=None, wa=None, num_cycles=2, temp=5000, harmCA=None) :
    import os
    print "CNS Annealing for %d cycles starting at %dK" % (num_cycles, temp)
    print sg, reso
    anneal_inp = "%s/appscripts-acer/anneal2.inp" % os.environ["RTKROOT"]
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
        #elif re.compile("igroup").search(l) and len(segids) > 0 : inputs.append( makeIgroupLine(segids," &atom_select and ") )
        else : inputs.append(re.sub("\n", "", l))
        #elif re.compile("atom_harm=").search(l) : inputs.append("atom_harm = ( name ca );")
    fp = open('CNSanneal.inp', 'w')
    for l in inputs : print >> fp, l
    fp.close()
    if logfile == None : logfile = "anneal.log"
    #execCmd("/home/anjum/downloads/cns_solve_1.1/intel-i686-linux_g77/bin/cns_solve > "+logfile, inputs)
    import os # os.system
    #import sys ; sys.exit()
    #os.system("/home/ak459//downloads/CNS/cns_solve_1.2/intel-x86_64bit-linux/bin/cns_solve  > CNSanneal.inp %s"%logfile)
    execCmd("/home/ak459//downloads/CNS/cns_solve_1.2/intel-x86_64bit-linux/bin/cns_solve  > "+logfile, inputs)
    #import sys ; sys.exit()






    
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
