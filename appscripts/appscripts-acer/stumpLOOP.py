
#from pref import removeMODEL,removeCRYST
import sys
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


import prepareChain


def getCRYST(pdbfile):
    a,b,c,alpha,beta,gamma = None,None,None,None,None,None
    import re
    x = []
    p = re.compile('\S+')
    for l in open(pdbfile,'r').readlines() :
        chara = len(l)
        if l[0:5] in "CRYST1":
            pl = l[8:chara]
            x = p.findall(pl)
            print x
        else : continue
    return float(x[0]),float(x[1]),float(x[2]),float(x[3]),float(x[4]),float(x[5])



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
    print "scMissInds",scMissInds
    for ri in scMissInds : buildAcb(resC, residsC, resnumsC, resnsC, chidsC, inscodesC, ptsC, ri)






def sortedDictValues1(adict):
    items = adict.items()
    items.sort()
    return [value for key, value in items]

def sortedDictValues2(adict):
    keys = adict.keys()
    keys.sort()
    return [dict[key] for key in keys]

def parseLoopres(pdbfile,loopres,chains):
    loopresids = []
    from pref import parseResnums
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    lines = open(loopres, 'r').readlines()
    import re
    allis = {} ; chainids = {}
    bigset = {}

    for c in chains :
        bigset[c] = {}
    stop = None
    for li in range(len(lines)) :
        l = re.sub("\n", "", lines[li])
        if l == '' or l[0] == '#' : continue
        if l[0:5] == "loop " :
            
            rids,chain = parseResnums(l)
            #print rids,chain 
            dict2add = bigset[chain]
            dict2add[int(rids[0])] = int(rids[1])

            #allis[int(rids[0])] = int(rids[1])
            #chainids[int(rids[0])] = chain

            
    #allis = sortedDictValues1(allis)
    #stop = 0
#    print bigset
    for cc , allis in bigset.items():
        stop = None 
        kk = allis.keys()
        kk.sort()
        #print "load ", pdbfile
        if stop == None :
            stop = allis[kk[0]]
            


        for k in range(1,len(kk)):
            akey = kk[k]
            print "loop '%s' [%d][%d] "%(cc,stop+1,akey-1)
            #print "color red, resi %d:%d "%(stop+1,akey-1)
            stop  = allis[akey]
            
            
        
def parseResnums(l) :
    start, resnums = None, []
    for i in range(len(l)) :
        if l[i] == '[' : start = i
        elif l[i] == ']' :
            resnums.append( l[start+1:i] )
    for i in range(len(l)) :
        if l[i] == '\'' and l[i+2] == '\'' :
            chain = l[i+1]
        
#print resnums
    #sys.exit()
    return resnums,chain

            

    
def main() :
    import re,sys
    import optparse ; parser = optparse.OptionParser()
    ## for subset analysis

    parser.add_option("--pml", action='store', type='string', dest='pml', help='pymol file', default=None)
    parser.add_option("--pdb", action='store', type='string', dest='pdb', help='pymol file', default=None)
    parser.add_option("--chain", action='store', type='string', dest='chain', help='pymol file', default=None)
    

    (options, args) = parser.parse_args()
    parseLoopres(options.pdb,options.pml,options.chain)


    

    


def getTime(sf,pdb):
    xrayRestGen = []
    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    esmin, esmax, esmean, rcmult, xscoreCutoff = .000, 5., .0, 5, 0.9
    pis = VecInt(resids.keys())
    er = EDrestraint.makeEDrestraintFromMap(pis, "Xranker", phasedmtz, esmin, esmax, smax)
    score = er.scoreAll(VecVecFloat(pts))
    print score

    
def getBeets(ss,pdb,corr):
    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    ip = open(ss, 'r')
    fileList = ip.readlines()
    import re,sys
    
    p = re.compile('(\s+\d+\s+)(\d+)(\s+\d+)(\s+\d+)(\s+\d+\s+)(\d+)(\s+\S+\s+)(\d+\s+)(\d+\s+\d+\s+\d+\s+\d+)(\-*\+*)')
    
    pps = re.compile('(\s*)(\d+)\s+(\d+)\s+(\d+)\s+(\d+)')

    sheet = {}
    sp =  {}
    sc = {}
    corrcheck = []
    rev_resnum = {}

    for k,v in resnums.items():
        rev_resnum[int(v)] = int(k)

    sheetnumber = 1
    strand = {}
    numPartnerArray = {}
    partnerArray = {}
    for k in fileList:

        strandNum = int(k[0:4])
        start = int(k[6:10])
        stop = int(k[13:17])
        sheetnumber = int(k[20:22])
        numberOfRes = int(k[28:30])

        if sheetnumber not in  sheet.keys():
            sheet[sheetnumber] = []
            sc[sheetnumber] = {}
            sp[sheetnumber] = {}
            
	numP = int(k[57:59])
        pArray = k[61:75] 
        direction = k[75:79]

        

        strand[int(strandNum)] = range(int(start),int(stop))

        rstart = resids[rev_resnum[int(start)]]
        rstop = resids[rev_resnum[int(stop)]]

        numPartnerArray[int(strandNum)] = numP
        

        
        partners = pps.match(pArray)
        partnerlist = list(partners.groups())

        for x in partnerlist :
            if x == ' ' or x == '':
                partnerlist.remove(x)
            else:
                int(x)

        partnerArray[int(strandNum)] = partnerlist
        dd = {}
        assert(len(direction)==4)

        for d in range(len(direction)):
            dir = direction[d]
            if dir == ' ':
                continue
            
            corrPartner = partnerlist[d]
            
            if dir == '-':dd[corrPartner] = "anti"
            elif dir == '+' : dd[corrPartner] = "parallel"
            
        temp = {}
        temp[int(strandNum)] =  [start,stop]
        sheet[sheetnumber].append(temp)  

        
        
        temp = {}
        temp[int(strandNum)]  =  dd
        sp[sheetnumber][strandNum] = dd
        

        temp = {}
        temp[int(strandNum)]  =  partnerlist
        sc[sheetnumber][int(strandNum)] = partnerlist



    ## order sheets
    newsheet = {}
    for sheetnum,strands in sheet.items():
        newstrands,tempstrands = [],[]



        for strandy in strands:
            tempstrands.append(strandy) 
            newstrands.append({})



        firststrandnotdone = -1
        laststrandnotdone = -1
        for strandy in strands:
            for strandnum, strandinfo in strandy.items():
                np = numPartnerArray[strandnum]
                if np == 1 and firststrandnotdone == -1:
                    newstrands[0] = strandy
                    firststrandnotdone = strandy.keys()[0]
                    #print "ADD 1",strandy
                    tempstrands.remove(strandy)
                    continue
                if np == 1 and laststrandnotdone == -1:
                    newstrands[len(strands)-1] = strandy
                    laststrandnotdone = strandy.keys()[0]
                    #print "ADD L ",strandy
                    tempstrands.remove(strandy)

        counter = 0
        while len(tempstrands) != 0:
            for strandy in tempstrands:
            
                strandk = strandy.keys()[0]

                ##print "sk",strandk,partnerArray[firststrandnotdone],firststrandnotdone
                if str(strandk) in partnerArray[firststrandnotdone]:
                    newstrands[1+counter] = strandy
                    firststrandnotdone = strandk
                    counter = counter + 1
                    tempstrands.remove(strandy)
                    break
                    
        for xx in newstrands:
            if len(xx.keys()) == 0 :
                print newstrands

        newsheet[sheetnum] = newstrands

   
    for sheetnumber,strandinfo in newsheet.items():
        print
        print "sheet ",
        for strandss in strandinfo:
            for pair in strandss.values():
                rstart = resids[rev_resnum[int(pair[0])]]
                rstop = resids[rev_resnum[int(pair[1])]]
                print "[%s][%s]"%(pair[0],pair[1]),
        print

        print "strand-pair ",
        newpairorder,pairorder = [],[]
        for it in range(0,len(strandinfo)-1):
            currStrand = strandinfo[it]
            nextStrand = strandinfo[it+1]
            cs = currStrand.keys()[0]
            ##6   106    107    1   0   2  YW                        1   3   0   0   0+   99  0 GHGTTHGT   106  107
               

            ns = nextStrand.keys()[0]
            ## print sp[sheetnumber][cs],ns,sheetnumber,cs
            pairorder.append(sp[sheetnumber][cs][str(ns)])
            

        newpairorder.append(pairorder[0])
        print newpairorder[0],
        for it in range(1,len(pairorder)):
            currentStrand = pairorder[it]
            prevStrand = newpairorder[it-1]

            if prevStrand == 'parallel' and currentStrand == 'anti':
                print "anti",
                newpairorder.append('anti')
            if prevStrand == 'anti' and currentStrand == 'parallel':
                print "anti",
                newpairorder.append('anti')
            if prevStrand == 'parallel' and currentStrand == 'parallel':
                print "parallel",
                newpairorder.append('parallel')
            if prevStrand == 'anti' and currentStrand == 'anti':
                print "parallel",
                newpairorder.append('parallel')

#        for it in newpairorder:
#            print it,

        print 
        print "strand-correspondence ",

        fstrandss = strandinfo[0]
        v1 = int(fstrandss.values()[0][0])
        ##print "[%s]"%(v1),

        
        for cc in range(len(strandinfo)-1):
            strandss = strandinfo[cc]
            nextstrandss = strandinfo[cc+1]
            currStrand = strandinfo[cc]
            nextStrand = strandinfo[cc+1]
            cs = currStrand.keys()[0]
            ns = nextStrand.keys()[0]

            strandcorr = newpairorder[cc]
            
            s1start = int(strandss.values()[0][0])
            s1stop = int(strandss.values()[0][1])
            
            s2start = int(nextstrandss.values()[0][0])
            s2stop = int(nextstrandss.values()[0][1])
            
            #print "WW",s1start,s1stop, s2start, s2stop
            #print "R1", range(s1start,s1stop+1)
            #print "R2",range(s2start,s2stop+1),corr[51]
            matchfound = 0
            for ss in range(s1start,s1stop+1):
                
                if ss in corr.keys():
                    p1 = ss
                    listp2 = corr[ss]
                    for p2 in listp2:
                        if p2 in range(s2start,s2stop+1) and matchfound == 0 :
                            print "[%s][%s]"%(p1,p2),
                            matchfound = 1
                            break
                        
        print

#            if strandcorr == 'anti':
#                v1 = int(nextstrandss.values()[0][1])
#            else :
#                v1 = int(nextstrandss.values()[0][0])
#                
#            if cc == len(strandinfo)-2:
#                print "[%s]"%(v1),
#
#            else :
#                print "[%s][%s]"%(v1,v1),
                
#        for cc in range(len(strandinfo)-1):
#            strandss = strandinfo[cc]
#            nextstrandss = strandinfo[cc+1]
#            currStrand = strandinfo[cc]
#            nextStrand = strandinfo[cc+1]
#            cs = currStrand.keys()[0]
#            ns = nextStrand.keys()[0]
#            if sp[sheetnumber][cs][str(ns)] == "anti":
#                v1 = resids[rev_resnum[int(strandss.values()[0][0])]]
#                v2 = resids[rev_resnum[int(nextstrandss.values()[0][1])]]#
#
#                v1 = int(strandss.values()[0][0])
#                v2 = int(nextstrandss.values()[0][1])#
#
#
#            else :
#                v1 = resids[rev_resnum[int(strandss.values()[0][0])]]
#                v2 = resids[rev_resnum[int(nextstrandss.values()[0][0])]]#
#
#                v1 = int(strandss.values()[0][0])
#                v2 = int(nextstrandss.values()[0][0])#
#


        

        
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


    


    pdbstart = int(resnums[0])
    pdbstop = int(resnums[len(resnums)-1])


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

            if int(start) < pdbstart+2 :
                start = pdbstart+2

            if int(stop) > pdbstop-2 :
                stop = pdbstop -2
            #print "%s [%s][%s]"%(helixType,resids[rev_resnum[int(start)]],resids[rev_resnum[int(stop)]])
            print "%s [%s][%s]"%(helixType,start,stop)



def getHees2(hlx,pdb):
    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    rev_resnum = {}
    for k,v in resnums.items():
        rev_resnum[int(v)] = int(k)
    ip = open(hlx, 'r')
    fileList = ip.readlines()
    oplist = []
    pdbstart = int(resnums[0])
    pdbstop = int(resnums[len(resnums)-1])

    infoStart = 0
    infoStop = 0
    for k in fileList:
        if "data" in k:
            infoStart = 1
        if "Interactions" in k :
            infoStop = 1

        if len(k) > 20 and infoStart == 1 and infoStop == 0:
            
            
            strandNum = k[0:3]
            chain = k[4]
            start = k[8:11]
            stop = k[15:18]
            htype = k[20]
            if htype == 'H':
                helixType = 'helix 4_13'
            elif htype == 'G':
                helixType = 'helix 3_10'
                
            if int(start) < pdbstart+2 :
                start = pdbstart+2
                
            if int(stop) > pdbstop-2 :
                stop = pdbstop -2
        
            print "%s [%s][%s]"%(helixType,start,stop)
            l = [helixType,start,stop]
            oplist.append(l)
    return oplist

def removeHyd(pdb):
    out = "1.%s"%pdb
    print out
    prot = protein(pdb, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    from peptidebuild import ModelRenderer
    ModelRenderer(res, resns, chids, resnums, inscodes, [], out).render(pts)



def getCorres(hbdfile):
    ip = open(hbdfile, 'r')
    fileList = ip.readlines()
    corr = {}
    for k in fileList :
        p1 = int(k[5:10])
        p2 = int(k[23:27])
        
        if p1 in corr.keys():
            val = corr[p1]
            val.append(p2)
            corr[p1] = val
        else : 
            corr[p1] = [p2]

            
        if p2 in corr.keys():
            val = corr[p2]
            val.append(p1)
            corr[p2] = val
        else : 
            corr[p2] = [p1]
        

    return corr
        
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





def getBest(rlist,clist,blist,plist,corr):
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
