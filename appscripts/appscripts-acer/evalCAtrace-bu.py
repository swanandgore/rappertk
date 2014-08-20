''' this script evaluates the CA-trace performance. 
    given a template and traced conformers, it calculates ca rmsd, mainchain rmsd, all-atom rmsd, chi-deviation per conformer
'''

import math
from pdbr import protein
import prot2res
from geometry import calcDihed, VecFloat, dihedDiff, calcDist
from data import chiRes

def findChis(res, pts, ind, resname) :
    if not resname in chiRes.keys() :
        return []
    chis = []
    for chiANs in chiRes[resname] :
        
        chis.append( calcDihed( VecFloat(pts[ res[ind][chiANs[0]] ]), VecFloat(pts[ res[ind][chiANs[1]] ]), VecFloat(pts[ res[ind][chiANs[2]] ]), VecFloat(pts[ res[ind][chiANs[3]] ]) ) )
    return chis

def findPhiPsiOmega(res, pts, ind) :
    phi = calcDihed( VecFloat(pts[ res[ind-1][' C  '] ]), VecFloat(pts[ res[ind][' N  '] ]), VecFloat(pts[ res[ind][' CA '] ]), VecFloat(pts[ res[ind][' C  '] ]) )
    psi = calcDihed( VecFloat(pts[ res[ind][' N  '] ]), VecFloat(pts[ res[ind][' CA '] ]), VecFloat(pts[ res[ind][' C  '] ]), VecFloat(pts[ res[ind+1][' N  '] ]) )
    ome = calcDihed( VecFloat(pts[ res[ind][' CA '] ]), VecFloat(pts[ res[ind][' C  '] ]), VecFloat(pts[ res[ind+1][' N  '] ]), VecFloat(pts[ res[ind+1][' CA '] ]) )
    return phi,psi,ome

def rmsHelper(res, pts, ind, res1, pts1, ind1) :
    cadiff, cacnt, mcdiff, mccnt, aadiff, aacnt = 0.,0,0.,0,0.,0
    for name, ai in res[ind].items() :
        if name in res1[ind1].keys() :
            diff = calcDist(VecFloat(pts[res[ind][name]]),VecFloat(pts1[res1[ind1][name]]))
        else : continue
        if name == ' CA ' :
            cadiff, cacnt = cadiff + diff, cacnt + 1
        if name in [' N  ',' CA ',' C  ', ' O  '] :
            mcdiff, mccnt = mcdiff + diff, mccnt + 1
        aadiff, aacnt = aadiff + diff, aacnt + 1
    return cadiff, cacnt, mcdiff, mccnt, aadiff, aacnt

def comparePhiPsiOmegaChi(refprot, testprot, mconly=None, givenResids=None) :
    assert mconly == None or mconly == 1
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(refprot)
    resTest, residsTest, resnumsTest, resnsTest, chidsTest, inscodesTest, ptsTest = prot2res.readProtRes(testprot)
    resid2ind = {}
    for k,v in residsTest.items() : resid2ind[v] = k
    ## report whats in one prot and not other and vice versa
    rids, ridsTest = resids.values(), residsTest.values()
    for rid in rids :
        if not rid in ridsTest : print "WARNING : [%s] is in reference but not in test" % rid
    for rid in ridsTest :
        if not rid in rids : print "WARNING : [%s] is in test but not in reference" % rid
    ## if both proteins have residues with same resid, compare. else report the disparity and continue.
    ## if same chain exists for ind+1,ind-1 in both, compare phi-psi-omega. else just compare chi
    chi1diff, chi1diffcnt, chi1accu = 0, 0, 0
    chi12diff, chi12diffcnt, chi12accu = 0, 0, 0
    ri2aarms = {}
    rms = None
    for ind,rid in resids.items() :
        if givenResids != None and not rid in givenResids : continue
        if not rid in resid2ind.keys() :
            continue
        indTest = resid2ind[rid]
        cmpPSO, pso, psoTest = None, (9999.99,9999.99,9999.99), (9999.99,9999.99,9999.99)
        try :
            if chids[ind+1] == chids[ind] and chids[ind-1] == chids[ind] : cmpPSO = 1
        except KeyError : pass
        if cmpPSO :
            pso = findPhiPsiOmega(res, pts, ind)
            psoTest = findPhiPsiOmega(resTest, ptsTest, indTest)
        trms = rmsHelper(res, pts, ind, resTest, ptsTest, indTest)
        if rms == None : rms = list(trms)
        else :
            for ri in range(len(trms)) : rms[ri] = rms[ri] + trms[ri]
        ri2aarms[ind] = trms[4]/trms[5]
        chis, chisTest = [], []
        if not mconly :
            print resids[ind], "---------------------------------------------------------------------"
            print chiRes[resns[ind]]
            miss = 0
            for chiANs in chiRes[resns[ind]] :
                
                if chiANs[0] not in res[ind].keys() or  chiANs[1] not in res[ind].keys() or chiANs[2] not in res[ind].keys() or chiANs[3] not in res[ind].keys()  :
                    miss = 1
            if miss == 1 :
                chis = []
                chisTest = []

                
            else :
                chis = findChis(res, pts, ind, resns[ind])
                chisTest = findChis(resTest, ptsTest, indTest, resnsTest[indTest])

            if len(chis) > 0 :
                chi1diff += math.fabs( dihedDiff(chis[0],chisTest[0]) )
                chi1diffcnt += 1
                if math.fabs( dihedDiff(chis[0],chisTest[0]) ) < 40 :
                    chi1accu += 1
            if len(chis) > 1 :
                chi12diff += math.fabs( dihedDiff(chis[0],chisTest[0]) ) + math.fabs( dihedDiff(chis[1],chisTest[1]) )
                chi12diffcnt += 1
                if math.fabs( dihedDiff(chis[0],chisTest[0]) ) < 40 and math.fabs( dihedDiff(chis[1],chisTest[1]) ) < 40 : chi12accu += 1
        print ( "[%10s] %8.2f %8.2f %8.2f" + " %8.2f"*len(chis) )  %  tuple( [resids[ind]] + list(pso) + list(chis) )
        print ( "             %8.2f %8.2f %8.2f" + " %8.2f"*len(chisTest) )  %  tuple( list(psoTest) + list(chisTest) )
    if not mconly :
        print "CHI___1 %4d %5.2f (percentage of sidechains with <40 error, average error)" % (chi1accu*100/chi1diffcnt, chi1diff/chi1diffcnt) #, chi1diff, chi1diffcnt
        print "CHI_1+2 %4d %5.2f" % (chi12accu*100/chi12diffcnt, chi12diff/(2*chi12diffcnt)) #, chi12diff, chi12diffcnt
    print "CA rmsd", rms[0]/rms[1]
    print "MC rmsd", rms[2]/rms[3]
    if not mconly : print "AA rmsd", rms[4]/rms[5]
    ris, rmsds = list( ri2aarms.keys() ), list( ri2aarms.values() )
    for i in range(len(ris)) :
        for k in range(i+1,len(ris)) :
            if rmsds[i] < rmsds[k] :
                temp = rmsds[i] ; rmsds[i] = rmsds[k] ; rmsds[k] = temp
                temp = ris[i] ; ris[i] = ris[k] ; ris[k] = temp
    top10bigAArmsdRes = [] ; numresids = 10
    if len(ris) < 10 : numresids = len(ris)
    for i in range(numresids) : print i, numresids ; top10bigAArmsdRes.append( resids[ris[i]] )
    return chi1accu*100/chi1diffcnt, chi1diff/chi1diffcnt, chi12accu*100/chi12diffcnt, chi12diff/(2*chi12diffcnt), rms[0]/rms[1], rms[2]/rms[3], rms[4]/rms[5], top10bigAArmsdRes

if __name__ == "__main__" :
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--refpdb", action='store', type='string', dest='refpdbfile', help='reference PDB file containing one structure')
    parser.add_option("--testpdb", action='store', type='string', dest='testpdbfile', help='structures in this pdb file are to be compared to reference structure')
    parser.add_option("--mconly", action='store', type='int', dest='mconly', help='to compare mainchains only, set to 1. else all-atom comparison will be attempted.', default=None)

    (options, args) = parser.parse_args()
    if options.mconly != 1 :
        options.mconly = None
        

    
    refprot = protein(options.refpdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    alllines = open(options.testpdbfile, 'r').readlines()
    curlines = []
    for li in range(len(alllines)) :
        l = alllines[li]
        if not l : break
        if l[0:5] in ['MODEL','REMAR','USER '] : continue
        curlines.append(l)
        if l[0:6] == 'ENDMDL' or li == len(alllines)-1 :
            testprot = protein(curlines, read_hydrogens=0, read_waters=0, read_hets=0)
            curlines = []
            comparePhiPsiOmegaChi(refprot, testprot, options.mconly)
            #import sys; sys.exit(0)
