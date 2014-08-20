

#from pref import removeMODEL,removeCRYST
from pdbr import protein, makeResid
import prot2res
from builders import VecVecFloat, VecInt
import data
from geometry import CAtraceGH
import prepareChain
from data import vdwr
from misc import verbose
from protRefine import fixCNSop
import sys
from prepareChain import removeSC
import os, shutil, re
from geom import vec_dist
from multProtref import joinPDBs
import sys, os, shutil, re
from xray import cif2mtz, uniqueify, cns_generate, cns_anneal, mtz2hkl, sgCCP4toCNS
from xcheck import XrayRanker, XrayScorer
from pref import removeRElines,copyNonprotein
import os, shutil, sys, random, re
from xray import cif2mtz, sfall, uniqueify, findResolutionReflections
from pref import cnsRefinement,  changeBfacs, removeMODEL, adjustBfacOccu, copyNonprotein, removeRElines, randomize
from xcheck import XrayScorer
from pdbr import protein, isAAres, line2resid
import prot2res
from peptidebuild import ModelRenderer
from builders import VecInt, VecVecFloat, CBbuilder
from data import consts
from prepareChain import buildCBs
from pref import findChangedSC,adjustBfac
from protinfo import *
import pdbinfo
from geom import vec_dist
import os, shutil, re,sys
import geometry
from pdbr import protein, makeResid
import prot2res
from xray import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS, fft, omitmap, mapman
from builders import VecVecFloat, VecInt
from procrun import proc_run_exitOnError as execCmd
from xcheck import XrayScorer, XrayRanker
from data import sgtable
from evalCAtrace import comparePhiPsiOmegaChi

from loopbuild import Multiloop
import prepareChain

#REMARK   2 RESOLUTION.    2.10 ANGSTROMS.                                       
def getRESO(pdbfile):
    reso = None
    import re
    x = []
    p = re.compile('\S+')
    for l in open(pdbfile,'r').readlines() :
        chara = len(l)
        if l[0:19] in "REMARK   2 RESOLUTION":
            pl = l[8:chara]
            x = p.findall(pl)
            #print x
            return float(x[2])
            #print x
        elif  "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) " in l :
            p = re.compile('\d+.\d+')
            pl = l[8:chara]
            x = p.findall(pl)
            return float(x[0])
        else : continue
        
    return None

def removespace(s):
    cleansed_s = "";
    for letter in s :
        if letter == ' ' :
            continue
        else:
            cleansed_s = cleansed_s + letter
    return cleansed_s

def getCRYST(pdbfile):
    a,b,c,alpha,beta,gamma = None,None,None,None,None,None
    import re
    x = []
    p = re.compile('\S+')
    for l in open(pdbfile,'r').readlines() :
        chara = len(l)
        if l[0:5] in "CRYST1":
            cryst = l[0:6]
            a = float(removespace(l[6:15])) ; b = float(removespace(l[15:24])) ; c = float(removespace(l[24:33])) ; alpha = float(removespace(l[33:40])) ; beta= float(removespace(l[40:47])) ; gamma = float(removespace(l[47:54]))
            spacegroup = removespace(l[55:66])
            
            #pl = l[8:chara]
            #x = p.findall(pl)
            #sg = x[6]+x[7]

            return a , b ,c , alpha , beta , gamma , spacegroup
            #print float(x[0]),float(x[1]),float(x[2]),float(x[3]),float(x[4]),float(x[5]) , x[6],x[7] , sg
            #return float(x[0]),float(x[1]),float(x[2]),float(x[3]),float(x[4]),float(x[5]) , sg
            #print x
        else : continue
        
    return None,None,None,None,None,None, None







def getRESOLUTION(pdbfile):
    res = None
    import re
    x = []
    p = re.compile('\S+')
    for l in open(pdbfile,'r').readlines() :
        chara = len(l)
        if "REMARK   2 RESOLUTION.    " in l:
            pl = float(l[23:30])
            break
    return pl



def makeMaps(mtzfile,pdbfile,sg):
    a,b,c,alpha,beta,gamma = getCRYST(pdbfile)
    cif2mtz(mtzfile, "base.mtz", a, b, c, alpha, beta, gamma, sg)
    sfallA(pdbfile,"base.mtz" , "%s.phased.mtz" %pdbfile)
    phasedmtz = "%s.phased.mtz"%pdbfile
    
    
    ## 0 : 2F1-F2 , 1: F1 , 2: F2
    from restraints import EDrestraint
    esmin, esmax, esmean, rcmult = 0.0001, 5.0, 0.05, 10
    EDrestraint.makeMap(phasedmtz,"FP" ,"FC", "PHIC", 1 , esmin, esmax, esmax)



def makeMapFC(mtzfile):
    from restraints import EDrestraint
    esmin, esmax, esmean, rcmult = 0.0001, 5.0, 0.05, 10
    EDrestraint.makeFCMap(mtzfile,"FC", "PHIC", esmin, esmax, esmax)

def removeChainId(pdbfile):

    from pdbr import line2bfac, isPdbAtomLine, line2bfac, line2atomid, changeBfactor, line2crdstr
    from pdbinfo import segi
    newlines = []


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



def renameResid(pdbfile,outpdb,reference=None):
    removeChainId(pdbfile)

    if reference !=None:
        prot = protein(reference, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    protC = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    resC, residsC, resnumsC, resnsC, chidsC, inscodesC, ptsC = prot2res.readProtRes(protC)

    ip = open(pdbfile, 'r')
    fileList = ip.readlines()
    from pdbr import line2bfac, isPdbAtomLine, line2bfac, line2atomid, changeBfactor, changeResn
    counter = 0

    offset  = int(resnumsC[0]) - 1
    newlines = []
    
    for l in fileList:
        if isPdbAtomLine(l):
            newnum = int(l[23:27])- offset
            if newnum < 10:
                nl = l[0:23]+ "  "+str(newnum)+ " "  + l[27:]
            elif newnum < 100 : 
                
                nl = l[0:23]+ " "+str(newnum)+ " "  + l[27:]
            else :
                nl = l[0:23]+ ""+str(newnum)+ " "  + l[27:]
            
            newlines.append(nl)
        else:
            newlines.append(nl)

    op = open(outpdb, 'w')
    for l in newlines : op.write(l)
    op.close()
    from pdbr import checkRepairMissingSCatoms,missingSCatoms
    from prepareChain import buildAcb, incompleteSCcorrection,PrepareChain
    seq = checkRepairMissingSCatoms(newlines,"Y")
    scMissInds = incompleteSCcorrection(resC, resnsC, ptsC)

    for ri in scMissInds : buildAcb(resC, residsC, resnumsC, resnsC, chidsC, inscodesC, ptsC, ri)


    sys.exit()


    
def main() :

    import optparse ; parser = optparse.OptionParser()
    #parser.add_option("--listfile", action='store', type='string', dest='listfile', help='list of pdb structures to calculate precision', default=None)
    ## for subset analysis
    #parser.adBd_option("--verbose", action='store', type='int', dest='verbose', help='0 means least verbosity etc.', default=0)
    #parser.add_option("--rlist", action='store', type='string', dest='rlist', help='list of pdb structures to calculate precision', default=None)
    #parser.add_option("--clist", action='store', type='string', dest='clist', help='list of pdb structures to calculate precision', default=None)
    #parser.add_option("--plist", action='store', type='string', dest='plist', help='list of pdb structures to calculate precision', default=None)
    #parser.add_option("--blist", action='store', type='string', dest='blist', help='list of pdb structures to calculate precision', default=None)

    #parser.add_option("--mtz", action='store', type='string', dest='mtz', help='cif file to make map', default=None)

    #parser.add_option("--str", action='store', type='string', dest='sheet', help='str file from promotif', default=None)
    #parser.add_option("--hlx", action='store', type='string', dest='helix', help='str file from promotif', default=None)
    #parser.add_option("--pdb", action='store', type='string', dest='pdb', help='pdb file to make map', default=None)
    #parser.add_option("--sg", action='store', type='string', dest='sg', help='sg', default=None)
    #parser.add_option("--resolution", action='store', type='float', dest='resolution', help='resolution', default=None)
    #parser.add_option("--refpdb", action='store', type='string', dest='refpdb', help='reference ', default=None)
    parser.add_option("--loopseq", action='store', type='string', dest='outpdb', help='outfile', default=None)

    
    parser.add_option("--pdb", action='store', type='string', dest='pdb', help='pdbfile', default=None)
    parser.add_option("--start", action='store', type='int', dest='start', help='outfile', default=None)
    parser.add_option("--chainid", action='store', type='string', dest='chainid', help='outfile', default=None)
    parser.add_option("--stop", action='store', type='int', dest='stop', help='counter', default=0)
    #parser.add_option("--noe-nmr",action = 'store', type ='string',dest = 'noeFilename',help = 'file containing NMR restraints',default=None)
    (options, args) = parser.parse_args()
    import misc
    success = 0
    prot = protein(options.pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)


    for k,  v in resids.items():
        #print k,v,type(options.start),type(options.stop)
        if int(resnums[k])  >= options.start and int(resnums[k]) <= options.stop and chids[k] == options.chainid :
            #print k,v,options.start,options.stop
            success  = 1
            #print resnums[k]
            
    if success == 1 :
        print options.pdb,"S"
    else:
        print options.pdb,"F"
    import sys ; sys.exit()
            
        
    
    #getTransform(options.pdb1,options.pdb2)
    #transformNcopy(options.pdb,options.outpdb,options.count)
    #sys.exit()

    #makeMapFC(options.mtz)    


    ## transfor and copy 
    #transformNcopy(options.pdb,options.outpdb,options.count)
    #cryst = "CRYST1   48.00    48.00    48.000  90.00  90.0   90.00 P 1           1  "
    #addHeader(options.outpdb,cryst)
    #sys.exit()
    
    ### for removing hyd
    #removeHyd(options.pdb)
    #out = "1.%s"%options.pdb
    #removeMODEL(out)
    #cryst = "CRYST1   48.00    48.00    48.000  90.00  90.0   90.00 P 1           1  "
    #addHeader(out,cryst)
    #sys.exit()

    
    #removeChainId(options.pdb)
    ## getT(options.pdb,options.refpdb)



    #transformNcopy(options.pdb,options.outpdb,options.count)
    #cryst = "CRYST1   43.00    43.00    43.000  90.00  90.0   90.00 P 1           1  "
    #addHeader(options.outpdb,cryst)





    #removeMODEL(options.outpdb)
    #removeCRYST(options.outpdb)
    #cryst =  "CRYST1  132.00   132.00   132.000  90.00  90.0   90.00 P 1          27  "   
    #addHeader(options.pdb,cryst)
    
    #removeMODEL(options.pdb)
    #getT(options.pdb,options.refpdb)

    array = getBest(options.rlist,options.clist,options.blist,options.plist)
    
    


    
    pdblist = []
    import re
    p = re.compile('\S+')
    for l in open(options.listfile,'r').readlines() :
        m = p.match(l)
        fn =  m.group()
        pdblist.append(fn)
        resetBfac(fn)
        removeCRYST(fn)
        #cryst = "CRYST1   73.700   48.600   43.000  90.00  98.10  90.00 C 1 2 1       8  "
        cryst =  "CRYST1   38.00    38.00    38.000  90.00  90.0   90.00 P 1           1  "   
        addHeader(fn,cryst)
        removeMODEL(fn)
    #joinPDBs("x.pdb",pdblist)




def getTime(sf,pdb,sg):
    a,b,c,alpha,beta,gamma = getCRYST(pdb)
    cif2mtz(sf, "base.mtz", a, b, c, alpha, beta, gamma, sg)
    uniqueify("base.mtz", "rfree.mtz")
    sfall(pdb, "rfree.mtz", "phased.mtz")
    xrayRestGen = []
    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    esmin, esmax, esmean, rcmult, xscoreCutoff = .000, 5., .0, 5, 0.9
    pis = VecInt(resids.keys())
    folabel="FP"
    fclabel="FC"
    philabel="PHIC"
    maptype="2F1-F2"

    if maptype == '2F1-F2' : maptype = 0
    elif maptype == 'F1' : maptype = 1
    else : print "unknown maptype ", maptype ; sys.exit(1)
    er = EDrestraint.makeEDrestraintFromMTZ(pis, "EDrestraint", "phased.mtz", folabel, fclabel, philabel, maptype, esmin, esmax, esmean) 
    score = er.scoreAll(VecVecFloat(pts))

    sys.exit()
    
def getBeets(ss,pdb):
    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    ip = open(ss, 'r')
    fileList = ip.readlines()
    import re,sys

    p = re.compile('(\s+\d+\s+)(\d+)(\s+\d+)(\s+\d+)(\s+\d+\s+)(\d+)(\s+\S+\s+)(\d+\s+)(\d+\s+\d+\s+\d+\s+\d+)(\-*\+*)')

    pps = re.compile('(\d+)\s+(\d+)\s+(\d+)\s+(\d+)')
    sheet = {}
    sp =  {}
    sc = {}
    corrcheck = []
    rev_resnum = {}

    for k,v in resnums.items():
        rev_resnum[int(v)] = int(k)

    sheetnumber = 1
    strand = {}
    for k in fileList:
        params = p.match(k)
        ssParams = list(params.groups())

        #strandNum = ssParams[0]
        #start = ssParams[1]
	#stop = ssParams[2]
        #sheetnumber = int(ssParams[3])
        #numberOfRes = int(ssParams[5])


        #if sheetnumber not in  sheet.keys():
        #    sheet[sheetnumber] = []
        #    sc[sheetnumber] = []
        #    sp[sheetnumber] = []
            
	#numP = ssParams[7]
	#pArray = ssParams[8]
	#direction =ssParams[9]


        #######################################

        strandNum = int(k[0:4])
        start = int(k[6:10])
        stop = int(k[13:17])
        sheetnumber = int(k[20:22])
        numberOfRes = int(k[28:30])

        if sheetnumber not in  sheet.keys():
            sheet[sheetnumber] = []
            sc[sheetnumber] = []
            sp[sheetnumber] = []
            
	numP = int(k[57:59])
        pArray = k[62:75] 
        direction = k[75:79]


        strand[int(strandNum)] = range(int(start),int(stop))
        rstart = resids[rev_resnum[int(start)]]
        rstop = resids[rev_resnum[int(stop)]]

        dd = []
        for dir in direction:
            if dir == '-':dd.append("anti")
            elif dir == '+' : dd.append("parallel")
            
        temp = {}
        temp[int(strandNum)] =  [rstart,rstop]
        sheet[sheetnumber].append(temp)  
        
        temp = {}
        temp[int(strandNum)]  =  dd
        sp[sheetnumber].append(temp)

        temp = {}
        temp[int(strandNum)]  =  pArray 
        sc[sheetnumber].append(temp)

    for k,v in sheet.items():
        sc_to_write = []

        print         
        print
        print "sheet ",
        corrcheck = []
        for dp in v:
            for pair in dp.values():
                print "[%s][%s]"%(pair[0],pair[1]),
        print

        #print "strand-pair ",sp[k]
        #for z in range(1,len(sp[k])):
        #    print sp[k][z].values()[0],


        print "strand-pair ",

        for z in sc[k]:
            for a,b in z.items():
                p1 = int(strand[a][0])
                
                partners = pps.match(b)
                partnerlist = list(partners.groups())
                
                for ppp in range(len(partnerlist)):
                    partner = partnerlist[ppp]
                    ipartner = int(partner)
                    corrDoneA = "%s-%s"%(a,ipartner)
                    corrDoneB = "%s-%s"%(ipartner,a)
                    
                    if ipartner != 0 and (corrDoneA not in corrcheck) and (corrDoneB not in corrcheck):
                        corrcheck.append(corrDoneA)
                        corrcheck.append(corrDoneB)
                        p2 = int(strand[ipartner][0])
                        sc_to_write.append("[%s][%s]"%(resids[rev_resnum[p1]],resids[rev_resnum[p2]]))
                        for part1 in  sp[k]:
                            for key,val in part1.items():
                                if key == a:
                                    print part1[key][ppp],
                                    
                                    
        print
        print "strand-correspondence ",
        for gg in sc_to_write:
            print gg,
        print
        
def getHees(hlx,pdb):
    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    rev_resnum = {}
    for k,v in resnums.items():
        rev_resnum[int(v)] = int(k)
    ip = open(hlx, 'r')
    fileList = ip.readlines()
    import re,sys
    hlxpattern = re.compile('(\s+)(\d+)\s+(\d+)\s+(\d+)\s+(\S)')

    for k in fileList:
        params = hlxpattern.match(k)

        if params != None:
            hhParams = list(params.groups())
            strandNum = hhParams[1]
            start = hhParams[2]
            stop = hhParams[3]
            htype = hhParams[4]
            if htype == 'H':
                helixType = 'helix 4_13'
            elif htype == 'G':
                helixType = 'helix 3_10'
                
            print "%s [%s][%s]"%(helixType,resids[rev_resnum[int(start)]],resids[rev_resnum[int(stop)]])

        
def removeHyd(pdb):
    out = "1.%s"%pdb

    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    from peptidebuild import ModelRenderer
    ModelRenderer(res, resns, chids, resnums, inscodes, [], out).render(pts)



def getRanges(listfile,noefile):
    ip = open(listfile, 'r')
    fileList = ip.readlines()
    
    current = {}
    import re,sys
    p = re.compile('\S+')
    m = p.match(fileList[0])
    fn =  m.group()


    lcurrent = hyd2heavy(fn,noefile)
    ucurrent = hyd2heavy(fn,noefile)
    
    for k in fileList :
        iter  = {}
        m = p.match(k)
        fn =  m.group()
        iter = hyd2heavy(fn,noefile) 

        for a,v in iter.items():
            if iter[a] < lcurrent[a]:
                lcurrent[a] = iter[a]
            if iter[a] > ucurrent[a]:
                ucurrent[a] = iter[a]

    print 
    for k,v in lcurrent.items():
        print "%s,%s,%s"%(k,lcurrent[k],ucurrent[k])


def hyd2heavy(pdb,noeFilename):
    prot = protein(pdb, read_hydrogens=1, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    temp = {}
    for k,v in resns.items():
        if  v in temp.values():
            print "",
        else :
            temp[v] = (res[k])
    carbs = {}
    hyds = {}
    for k,v in temp.items():
        tt = []
        hh = []
        for a,b in v.items():
            if 'C' in a or 'N' in a or 'S' in  a or 'O' in a:
                tt.append(a)
            elif 'H' in a :
                hh.append(a)
        carbs[k] = tt
        hyds[k] = hh


    coder = {}
    for k,v in carbs.items():
        corrH = hyds[k]
        coder[k] = {}
        for c in v :
            match = 0
            for hx in corrH:
                if c[2] in hx and c[2]!= ' ' and c[3] == ' ' and c!= ' OH ':
                    print "",
                    coder[k][hx] = c
                    match  = 1

                elif hx == " H  " and c!= ' OH ':
                    coder[k][' HN '] = ' N  '

                elif c[2:4] in hx and c[2]!= ' ' and c!= ' OH ':
                    print "",
                    coder[k][hx] = c
                    match = 1

                elif hx == " HH " and c== ' OH ':
                    coder[k][hx] = c

    #for a,b in hyds.items():
    #    for k in b:
    #        print a,coder[a][k],k


    noe_ = noe(noeFilename)
    partner1, group1 , partner2,  group2 ,lower ,upper  = noe2restraint.parseNoe(noe_)

    rev_resnum = {}
    for k,v in resnums.items():
        rev_resnum[int(v)] = int(k)

    parsed = {}
    for k,v in partner1.items():
        a1 = group1[k]
        a2 = group2[k]
        if len(a1) == 1:   a1 = " "+a1+ "  "
        elif len(a1) == 2: a1 = " "+a1+ " "
        elif len(a1) == 3: a1 = " "+a1

        if len(a2) == 1:   a2 = " "+a2+ "  "
        elif len(a2) == 2: a2 = " "+a2+ " "
        elif len(a2) == 3: a2 = " "+a2

        n1= partner1[k]
        n2= partner2[k]

        nres1 = rev_resnum[n1]
        nres2 = rev_resnum[partner2[k]]
        
        aares1 = resns[rev_resnum[n1]]
        aares2 = resns[rev_resnum[partner2[k]]]

        if (aares1 == 'PHE' or aares1 == 'LEU' or aares1 == 'LYS' or aares1 == 'GLN' or aares1 == 'SER' or aares1 == 'ASP'  or aares1 == 'ASN' or aares1 == 'HIS' or aares1 == 'TYR'  or aares1 == 'GLU'  or aares1 == 'TRP' or aares1 == 'CYS'  or aares1 == 'ARG'  or aares1 == 'MET' or aares1 == 'PRO')  and a1 == ' HB1' :  a1 = ' HB2' 
        elif (aares1 == 'PHE' or aares1 == 'LEU' or aares1 == 'LYS' or aares1 == 'GLN' or aares1 == 'SER'or aares1 == 'ASP' or aares1 == 'ASN' or aares1 == 'HIS' or aares1 == 'TYR'  or aares1 == 'GLU'  or aares1 == 'TRP' or aares1 == 'CYS'  or aares1 == 'ARG'  or aares1 == 'MET' or aares1 == 'PRO') and a1 == ' HB2' :  a1 = ' HB3'

        elif (aares1 == 'MET'  or aares1 == 'GLU' or aares1 == 'LYS'  or aares1 == 'ARG' or aares1 == 'PRO'  or aares1 == 'GLN')and a1 == ' HG1' :  a1 = ' HG2'
        elif (aares1 == 'MET'  or aares1 == 'GLU' or aares1 == 'LYS'  or aares1 == 'ARG' or aares1 == 'PRO'  or aares1 == 'GLN')and a1 == ' HG2' :  a1 = ' HG3'

        elif (aares1 == 'GLY')and a1 == ' HA1' :  a1 = ' HA2'
        elif (aares1 == 'GLY')and a1 == ' HA2' :  a1 = ' HA3'

        elif (aares1 == 'LYS' )and a1 == ' HE1' :  a1 = ' HE2'
        elif (aares1 == 'LYS')and a1 == ' HE2' :  a1 = ' HE3'

        
        elif (aares1 == 'LYS' or aares1 == 'PRO' or aares1 == 'ARG')and a1 == ' HD1' :  a1 = ' HD2'
        elif (aares1 == 'LYS' or aares1 == 'PRO' or aares1 == 'ARG')and a1 == ' HD2' :  a1 = ' HD3'


        elif (aares1 == 'ILE')and a1 == 'HG11' :  a1 = 'HG12'
        elif (aares1 == 'ILE')and a1 == 'HG12' :  a1 = 'HG13' 
        else :
            print "",
            
        if (aares2 == 'PHE' or aares2 == 'LEU' or aares2 == 'LYS' or aares2 == 'GLN'  or aares2 == 'SER' or aares2 == 'ASP' or aares2 == 'ASN' or aares2 == 'HIS'  or aares2 == 'TYR'  or aares2 == 'GLU'  or aares2 == 'TRP' or aares2 == 'CYS'  or aares2 == 'ARG'  or aares2 == 'MET' or aares2 == 'PRO') and a2 == ' HB1' :  a2 = ' HB2' 
        elif (aares2 == 'PHE' or aares2 == 'LEU' or aares2 == 'LYS' or aares2 == 'GLN'  or aares2 == 'SER' or aares2 == 'ASP' or aares2 == 'ASN' or aares2 == 'HIS'  or aares2 == 'TYR'  or aares2 == 'GLU'  or aares2 == 'TRP' or aares2 == 'CYS' or aares2 == 'ARG'  or aares2 == 'MET' or aares2 == 'PRO')  and a2 == ' HB2' :  a2 = ' HB3'






        elif (aares2 == 'MET'  or aares2 == 'GLU' or aares2 == 'LYS'  or aares2 == 'ARG' or aares2 == 'PRO'  or aares2 == 'GLN')and a2 == ' HG1' :  a2 = ' HG2'
        elif (aares2 == 'MET'  or aares2 == 'GLU'or aares2 == 'LYS'   or aares2 == 'ARG' or aares2 == 'PRO'  or aares2 == 'GLN')and a2 == ' HG2' :  a2 = ' HG3'

        elif (aares2 == 'LYS'  or aares2 == 'PRO' or aares2 == 'ARG' )and a2 == ' HD1' :  a2 = ' HD2'
        elif (aares2 == 'LYS' or aares2 == 'PRO' or aares2 == 'ARG')and a2 == ' HD2' :  a2 = ' HD3'

        elif (aares2 == 'GLY')and a2 == ' HA1' :  a2 = ' HA2'
        elif (aares2 == 'GLY')and a2 == ' HA2' :  a2 = ' HA3'
        elif (aares2 == 'ILE')and a2 == 'HG11' :  a2 = 'HG12'
        elif (aares2 == 'ILE')and a2 == 'HG12' :  a2 = 'HG13'

        elif (aares2 == 'LYS')and a2 == ' HE1' :  a2 = ' HE2'
        elif (aares2 == 'LYS')and a2 == ' HE2' :  a2 = ' HE3'
        else :
            print "",
        
            
        #print aares1,coder[aares1],a1
        #print coder[aares1][a1],aares1

        #print aares2,coder[aares2],a2
        #print coder[aares2][a2],aares2

        p2  =  pts[res[nres2][coder[aares2][a2]]]
        p1 = pts[res[nres1][coder[aares1][a1]]]

        c2  =  coder[aares2][a2]
        c1 =   coder[aares1][a1]
        
        from geom import vec_dist
        d = vec_dist(p1,p2) 
        kiy = "%s,%s,%s,%s"%(partner1[k],c1,partner2[k],c2)
        parsed[kiy] = d
    return  parsed



def getTransform(pdb,refpdb):
    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    
    protf = protein(refpdb, read_hydrogens=0, read_waters=0, read_hets=0)
    resf, residsf, resnumsf, resnsf, chidsf, inscodesf, ptsf = prot2res.readProtRes(protf)

    for j in range(len(pts)):
        x = pts[j][0] - ptsf[j][0]
        y = pts[j][1] - ptsf[j][1]
        z = pts[j][2] - ptsf[j][2]
        print x,y,z
    return x,y,z

def transformNcopy(pdb,out,count=None):
    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    #a,b,c,alpha,beta,gamma = getCRYST(pdb)
    X = 44
    series = ((X,X,X),      #1
              (-X,-X,-X),  #  2

              (-X,+X,+X),   # 3
              (+X,-X,-X),   # 4

              (+X,-X,+X),   # 5  
              (-X,+X,-X),   # 6



              (X,0,0),    # 7
              (-X,0,0), # 8
              (0,X,0),  # 9 
              (0,-X,0), # 10 
              (0,0,X),  # 11
              (0,0,-X), #  12  

              (X,X,0),  # 13
              (0,X,X), # 14 
              (X,0,X), # 15

              (-X,-X,0),  # 16
              (0,-X,-X),  # 17
              (-X,0,-X),   # 18

              (-X,+X,0),  # 19
              (0,-X,+X),  # 20
              (-X,0,+X),   # 21

              (+X,-X,0),  # 22
              (0,+X,-X),  # 23
              (+X,0,-X) ,  # 24            

              (+X,+X,-X),   # 25  
              (-X,-X,+X),   # 26
              (0,0,0),   # 26  
              )

    X = 8
    newpts = []
    for pp in (range(len(pts))):
        newpts.append([0.,0.,0.])
        
    for pp in (range(len(pts))):
        newpts[pp][0] = pts[pp][0] -48.524 # + 5#  + 5  #- 47  #+ 20  #   +5 #+ X #series[count][0]
        newpts[pp][1] = pts[pp][1] -29.858 #+ 2 #+ 5   #+ X #series[count][1]
        newpts[pp][2] = pts[pp][2]  -45.63 # - 5  #+ 2 # - 15 #+ 10  #+ X #series[count][2]
        print pts[pp],newpts[pp]
        
        
    from peptidebuild import ModelRenderer
    ModelRenderer(res, resns, chids, resnums, inscodes, [], out).render(newpts)

def getRfromCNS(pdbfile,mtz,sg,resolution):
    cnsArgs = {}

    for cycle in range(20) : cnsArgs[cycle] = {} ; cnsArgs[cycle]["num_cycles"] = 2 ; cnsArgs[cycle]["temperature"] = 5000

    ## for cycle in range(20) : cnsArgs[cycle] = {} ; cnsArgs[cycle]["num_cycles"] = 1 ; cnsArgs[cycle]["temperature"] = 1
    
    a,b,c,alpha,beta,gamma = getCRYST(pdbfile)
    cnsout =      "dummy.pdb"
    phasedmtz =   "dummy.mtz"
    cycle = 0
    uniqueify(mtz, "rfree.mtz")
    
    cnsRefinement("rfree.mtz", pdbfile, phasedmtz, cnsout,a, b, c, alpha, beta, gamma, sg, resolution,cnsArgs, cycle)




def addHeader(pdbfile,header):
    ip = open(pdbfile, 'r')
    fileList = ip.readlines()
    newline = []
    for l in fileList:
        newline.append(l)

    ip.close()
    
    op = open(pdbfile, 'w')
    op.write(header)
    op.write("\n")
    for l in newline:
        op.write(l)
    op.close()





def getBest(rlist,clist,blist,plist):
    ip = open(rlist, 'r')
    fileList = ip.readlines()
    import re,sys
    p = re.compile('\S+')
    rama = {}
    irama = {}
    ip.close()
    for l in fileList:
        m = p.match(l)
        fn =  m.group()
        x= p.findall(l)
        rama[x[5]] = x[0]
        irama[x[0]] = x[5]




    ip = open(plist, 'r')
    fileList = ip.readlines()
    import re,sys
    p = re.compile('\S+')
    pack = {}
    ipack = {}
    ip.close()
    for l in fileList:
        m = p.match(l)
        fn =  m.group()
        x= p.findall(l)
        pack[x[6]] = x[0]#
        ipack[x[0]] = x[6]


    ip = open(clist, 'r')
    fileList = ip.readlines()
    import re,sys
    p = re.compile('\S+')
    chi = {}
    ichi = {}
    ip.close()
    for l in fileList:
        m = p.match(l)
        fn =  m.group()
        x= p.findall(l)
        chi[x[5]] = x[0]
        ichi[x[0]] = x[5]

        

    ip = open(blist, 'r')
    fileList = ip.readlines()
    import re,sys
    p = re.compile('\S+')
    bbc = {}
    ibbc = {}
    ip.close()
    for l in fileList:
        m = p.match(l)
        fn =  m.group()
        x= p.findall(l)
        bbc[x[4]] = x[0]
        ibbc[x[0]] = x[4]

    rkey = rama.keys()
    rkey.sort()

    ckey = chi.keys()
    ckey.sort()


    pkey = pack.keys()
    pkey.sort()


    bkey = bbc.keys()
    bkey.sort()

    pmode = pkey[50]
    cmode = ckey[50]
    bmode = bkey[50]
    rmode = rkey[50]

    rbest , pbest , cbest , bbest = [],[],[],[]
    for xx in rkey:
        if xx > rmode:
            rbest.append(rama[xx])

    print "\n***********************************************"
    for xx in ckey:
        if xx > cmode:
            cbest.append(chi[xx])
            
    print "\n***********************************************"
    
    for xx in pkey:
        if xx > pmode:
            pbest.append(pack[xx])
            
    print "\n***********************************************"
    for xx in bkey:
        if xx > bmode:
            bbest.append(bbc[xx])
    print "\n***********************************************"


    for bb in bbest:
        if bb in pbest and bb in cbest and bb in rbest:
            print bb
            print irama[bb],ipack[bb],ichi[bb],ibbc[bb]
            
def atom2PRO(pdbfile,outpdb):
    from pdbr import line2bfac, isPdbAtomLine, line2bfac, line2atomid, changeBfactor, line2crdstr,changeALTLOC, line2resid,line2resnum,line2allid,isOXTLine,OXT2OT

    from pdbinfo import segi
    newlines = []


    ip = open(pdbfile, 'r')
    fileList = ip.readlines()
    for l in fileList:

        if isOXTLine(l):
            l = OXT2OT(l)
            
        if isPdbAtomLine(l):
            l  =  l[0:68] + '    ' + "PRO \n"
            newlines.append(l)
        else :
            print l
            #newlines.append(l)
            
    op = open(outpdb, 'w')
    for l in newlines : op.write(l)
    op.close()


if __name__ == "__main__" : main()
