from pdbr import protein
import prot2res
from data import resAtoms , three2one
from geometry import calcDist, VecFloat
import sys


def check(res, resids, resnums, resns, chids, inscodes, pts, cadistCutoff=4.,null=[],printsum=1,mconly=None) :
    chains = set(chids.values())
    chad = []
        
    mcmiss, scmiss, chbr = set(), set(), set()
    
    inbreak = "0"
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
            

            if  res[k].has_key(' CA ') and inbreak  == "1" :
                inbreak = "0"
                chbr.add(k)



            if  res[k].has_key(' CA ') and not res[k+1].has_key(' CA ') :
                if inbreak == "0" : 
                    chbr.add(k) ;
                    inbreak = "1"


            if   resids[k+1] in null and resids[k] not in null :
                chbr.add(k) ; 
                
            if   resids[k] in null and resids[k+1] not in null :
                chbr.add(k+1) ; 


            if   res[k].has_key(' CA ') and  res[k+1].has_key(' CA ') and calcDist( VecFloat(pts[res[k][' CA ']]), VecFloat(pts[res[k+1][' CA ']]) ) > cadistCutoff and resids[k+1] not in null and resids[k] not in null :
                
                chbr.add(k) ; chbr.add(k+1) ;


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
    for x in mcmiss :
        if x not in lmcmiss : 
            lmcmiss.append(x)
    for y in scmiss :
        if y not in lscmiss : 
            lscmiss.append(y)
    for z in chbr : 
        if z not in lchbr :
            lchbr.append(z)

    return lmcmiss, lscmiss, lchbr




def getbands(res, resids, resnums, resns, chids, inscodes, pts, cadistCutoff=4.,null=[],printsum=1,mconly=None) :
    chains = set(chids.values())
    chad = []
    bands = {} 
    mcmiss, scmiss, chbr = set(), set(), set()
    for xd in range(1000):
        bands[xd] = [None,None]
    
    inbreak = "0"
    nextband  = -1 ;
    for chid in chains :
        keys = res.keys() ; keys.sort()
        start, end = None, keys[len(keys)-1]
        
        for k in keys :
            if (chids[k] == chid and start ==None) :
                start = k
            elif chids[k] != chid and start!=None and end == keys[len(keys)-1] :
                end = k-1

        if start == None or end == None :
            print "Error reading chains from coordinate file "
            import sys ; sys.exit()


        foundAnchor = 0 ; anchor  = start
        for rr in range(start,end):
            if ' CA ' in res[rr].keys() and  ' CA ' in res[rr+1].keys() :
                foundAnchor = 1
                anchor = rr
                break
        if foundAnchor == 1 :
            print "Anchor residues for bootstrapping OK",resids[anchor]
        else:
            print "Error : Chain cannot be bootstrapped"

            
        nextband = nextband + 1
        bands[nextband] = [anchor,None]

        
        for k in range(anchor,end+1) :
            assert chids[k] == chid
            mcok, scok = [], []

            if k == end or  resns[k] not in three2one.keys():
                continue

            ########
            if  res[k].has_key(' CA ') and inbreak  == "1" :
                inbreak = "0"
                nextband = nextband + 1 
                bands[nextband] =  [None,None]
                bands[nextband][0] = k


            if  res[k].has_key(' CA ') and not res[k+1].has_key(' CA ') :
                if inbreak == "0" : 
                    inbreak = "1"
                    bands[nextband][1] = k
                    
                

            if   resids[k+1] in null and resids[k] not in null :
                bands[nextband][1] = k

            if   resids[k] in null and resids[k+1] not in null :
                nextband = nextband + 1 
                bands[nextband] =  [None,None]
                bands[nextband][0] = k + 1
######
                
            if   res[k].has_key(' CA ') and  res[k+1].has_key(' CA ') and calcDist( VecFloat(pts[res[k][' CA ']]), VecFloat(pts[res[k+1][' CA ']]) ) > cadistCutoff and resids[k+1] not in null and resids[k] not in null :
                
                bands[nextband][1] = k 
                nextband = nextband + 1 
                bands[nextband] =  [None,None]
                bands[nextband][0] = k +1

        if bands [nextband][1] == None : 
            bands[nextband][1] = end

    return bands



def check3(res, resids, resnums, resns, chids, inscodes, pts, cadistCutoff=4.,null=[],printsum=1,mconly=None,badresids=[],modelIn=None,hresids = []) :
    print badresids
    prot = protein(modelIn, read_hydrogens=0, read_waters=0, read_hets=0)
    ores, oresids, oresnums, oresns, ochids, oinscodes, opts , oatomids , obfac = prot2res.readProtRes2(prot)
    dex = {}
    for k , v in oresids.items() :
        dex[v] = k
    for k , v in resids.items() :
        
        if v in badresids : 
            if   res[k].has_key(' CA ') and  res[k+1].has_key(' CA ')  :
                dist =  calcDist( VecFloat(pts[res[k][' CA ']]), VecFloat(pts[res[k+1][' CA ']]) ) 
                if dist < 3.8 : 
                    print "*", resnums[k] , dist,
                else :
                    print resnums[k] , dist,
                if v in oresids.values():
                    kd = dex[v]
                    print resnums[kd] , calcDist( VecFloat(opts[ores[kd][' CA ']]), VecFloat(opts[ores[kd+1][' CA ']]) )  , 
                    if  v in hresids :
                        print "H"
                    else :
                        print "C"
                else :
                    if  v in hresids :
                        print "H"
                    else :
                        print "C"
                    

def check2(res, resids, resnums, resns, chids, inscodes, pts, cadistCutoff=4.,null=[],printsum=1,mconly=None) :
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

        inbreak = "0"
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

            
            if  ' CA ' in res[k].keys() and inbreak=="1" :
                inbreak = "0"
                chbr.add(k)




            if not res[k+1].has_key(' CA ') :
                if inbreak == "0" : 
                    chbr.add(k) ;
                    inbreak = "1"

                    if printsum == 1:
                        print "\nCannot determine chainbreak due to missing CA atoms in ", resids[k], resids[k+1]
            else :
                print "KK", resids[k+1],res[k]

            
            if  res[k].has_key(' CA ') and res[k+1].has_key(' CA ') :
                if calcDist( VecFloat(pts[res[k][' CA ']]), VecFloat(pts[res[k+1][' CA ']]) ) > cadistCutoff and resids[k+1] not in null and resids[k] not in null :
                
                    chbr.add(k) ;
                    chbr.add(k+1) ;

                    if printsum == 1:
                        print "\nChainbreak between ", resids[k], resids[k+1]
                
                    chad.append(resids[k])

    lmcmiss , lscmiss , lchbr = [] , [], []


    for x in mcmiss :
        lmcmiss.append(x)
    for y in scmiss :
        lscmiss.append(y)
    for z in chbr : 
        lchbr.append(z)
    import sys; sys.exit()
    return lmcmiss, lscmiss, lchbr

if __name__ == "__main__" :
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--pdb", action='store', type='string', dest='pdbfile', help='template PDB file')
    parser.add_option("--refpdb", action='store', type='string', dest='refpdbfile', help='template PDB file')
    parser.add_option("--ssfile", action='store', type='string', dest='ssfile', help='template PDB file')
    parser.add_option("--cacaCutoff", action='store', type='float', dest='ca', help='min dist between adjacent CA to detect a chain-break', default=4.)
    (options, args) = parser.parse_args()
    prot = protein(options.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    #check(res, resids, resnums, resns, chids, inscodes, pts, options.cacaCutoff)

    if options.ssfile != None :
        ## only residues in loop regions are rebuilt
        from betaBuilderAK import parseSSfileAK
        helices , sheets = [], [] 
        helices , sheets = parseSSfileAK(options.ssfile,resids,resnums)
        hresids = []
        hids  = []
        starts, ends = [], []
        for strands, dirs, corrInds in sheets : 
            for beg,end in strands :
                starts.append(beg) ;
                ends.append(end)

        
        for hh in  helices:
            hstart = hh[1][0]
            hstop = hh[1][1]
            for hresidue in range(hstart,hstop+1):
                hresids.append(resids[hresidue])

            for si in range(len(starts)):
                for sresidue in range(starts[si],ends[si]):
                    hresids.append(resids[sresidue])

        print "adding strand", hresids
                    
        check3(res, resids, resnums, resns, chids, inscodes, pts, options.ca ,[],0,None,resids.values(),options.refpdbfile,hresids)

    else :
        check3(res, resids, resnums, resns, chids, inscodes, pts, options.ca,[],0,None,resids.values(),options.refpdbfile,[])
        
