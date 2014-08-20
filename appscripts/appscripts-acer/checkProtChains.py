from pdbr import protein
import prot2res
from data import resAtoms , three2one
from geometry import calcDist, VecFloat
import sys

def check(res, resids, resnums, resns, chids, inscodes, pts, cadistCutoff=4.,null=[],printsum=1,mconly=None) :
    chains = set(chids.values())
    chad = []
        
    mcmiss, scmiss, chbr = set(), set(), set()

    
    
    for chid in chains :
        keys = res.keys() ; keys.sort()
        start, end = None, keys[len(keys)-1]
        
        for k in keys :
            if (chids[k] == chid and start ==None) :
                start = k
                
            elif chids[k] != chid and start!=None and end == keys[len(keys)-1] :

                end = k-1
        if printsum == 1:
            print
            print "------------------ [%s] to [%s] ----------------" % (resids[start], resids[end]) ;

        #assert(start != None  and end !=None)
        if start == None or end == None :
            print "Error reading chains from coordinate file "
            import sys ; sys.exit()
        for k in range(start,end+1) :
            assert chids[k] == chid
            mcok, scok = [], []

            
            for an in resAtoms[resns[k]] :
                if not res[k].has_key(an) :
                    if resns[k] not in three2one.keys():
                        continue
                    if an in [' N  ',' CA ',' C  ',' O  '] :
                        mcok.append(an)
                    else :
                        scok.append(an)
            if len(mcok) > 0 :
                mcmiss.add(k) ;
                if printsum == 1:
                    print "Missing mainchain atoms in", resids[k], mcok
            if len(scok) > 0 :
                scmiss.add(k) ;
                if printsum == 1 and mconly != 1:
                    print "Missing sidechain atoms in", resids[k], scok

            if k == end :
                continue

            if resns[k] not in three2one.keys():
                continue
            
            if not res[k].has_key(' CA ') or not res[k+1].has_key(' CA ') :
                chbr.add(k) ;
                print "Adding chainbrak1",resids[k],resids[k+1],res[k]
                if printsum == 1:
                    print "\nCannot determine chainbreak due to missing CA atoms in ", resids[k], resids[k+1]

            elif calcDist( VecFloat(pts[res[k][' CA ']]), VecFloat(pts[res[k+1][' CA ']]) ) > cadistCutoff and resids[k+1] not in null and resids[k] not in null :
                
                chbr.add(k) ;
                print "Adding chainbrak2",resids[k]
                if printsum == 1:
                    print "\nChainbreak between ", resids[k], resids[k+1]
                
                chad.append(resids[k])
            #if  (int(resnums[k]) + 1) != int(resnums[k+1]) and chids[k] == chids[k+1] :
            #    if k not in chbr:
            #        chbr.add(k) 
            #        print "Adding chainbrak3",resids[k]
            #        if printsum == 1:
            #            print "\nPossible Chainbreak between", resids[k], resids[k+1]
            #    if k+1 not in chbr :
            #        chbr.add(k+1)
            #        print "Adding chainbrak4",resids[k+1]
    lmcmiss , lscmiss , lchbr = [] , [], []
#    print null

    for x in mcmiss :
        lmcmiss.append(x)
    for y in scmiss :
        lscmiss.append(y)
    for z in chbr : 
        lchbr.append(z)
#    import sys ; sys.exit()
    return lmcmiss, lscmiss, lchbr

if __name__ == "__main__" :
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='template PDB file')
    parser.add_option("--cacaCutoff", action='store', type='float', dest='pdbfile', help='min dist between adjacent CA to detect a chain-break', default=4.)
    (options, args) = parser.parse_args()
    prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    check(res, resids, resnums, resns, chids, inscodes, pts, options.cacaCutoff)
