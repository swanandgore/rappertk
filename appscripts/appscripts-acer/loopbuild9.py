import math, os, sys, optparse, random

from misc import verbose
from peptidebuild import ModelRenderer, BufferModelRenderer
from PopulationStrategy import PopulationStrategy
from prepareChainV5 import removeSC, mergeBRlists, addNdummyGly, addCdummyGly, printResAtoms, PrepareChain, prepareRNAchain, incompleteSCcorrection , incompleteMCSCcorrection  , incompleteMCcorrection , incompleteSCcorrection2
from pdbr import protein
import prot2res
from geometry import CAtraceGH
import data
from data import vdwr
from builders import VecVecFloat, VecInt

# locate bad regions, class them as loops, Cter, Nter, NtoC
# randomize their order and return

def locateRegionsRandomize3(resids, chids, badresids,modelN2C=1,buildType="notpoor") :
    #    tempresids = {} ; tempchids = {}
    #    for k,v in 
    bad = [None] * len( resids.keys() )
    for ri in resids.keys() :
        if resids[ri] in badresids : bad[ri] = 1
    #print bad
    keys = list( resids.keys() )
    keys.sort()
    loops, looptypes = [], []
    loopOn = None
    debug = 1
    for k in keys :
        #print k, bad[k],
        if bad[k] and not loopOn :
            #print 'here 1'
            loopOn = 1
            loops.append( [k,k] )
            if k+1 == len(keys) or chids[k+1] != chids[k] : ## check for ending singleton
                looptypes.append("Cter") ; loopOn = None ; continue
            if k-1 < 0 or chids[k-1] != chids[k] :
                if debug == 1 : print "T1"
                looptypes.append("Nter")
            else : looptypes.append("loop")
        elif loopOn and not bad[k] :

            loops[ len(loops)-1 ][1] = k-1 ; loopOn = None
        elif loopOn and ( k+1==len(keys) or chids[k+1] != chids[k] ) :

            if looptypes[ len(looptypes)-1 ] == "Nter" : looptypes[ len(looptypes)-1 ] = "NtoC"
            else : looptypes[ len(looptypes)-1 ] = "Cter"
            loops[ len(loops)-1 ][1] = k ; loopOn = None
#    print loops, "\n", looptypes
    for li in range(len(loops)) : ## extend loops on both sides by 1 when loopsize is 1
        if looptypes[li] != "loop" : continue
        loops[li][0] -= 1 ; loops[li][1] += 1
        ## danger commented amk
    print loops, "\n", looptypes
    for li in range(len(loops)) : ## loop -> ter for loops close to terminals
        if looptypes[li] != "loop" : continue
        for k in [1,2] :
            if (loops[li][0]-k < 0 or chids[loops[li][0]-k] != chids[loops[li][0]]) and looptypes[li] == "loop" :
                loops[li][0] -= (k-1) ; looptypes[li] = "Nter"
            if (loops[li][1]+k >= len(resids) or chids[loops[li][1]+k] != chids[loops[li][1]]) and looptypes[li] == "loop" :
                loops[li][1] += (k-1) ; looptypes[li] = "Cter"
    print loops, "\n", looptypes
    #for li in range(len(loops)) : ## Cter -> NtoC if only 1 good residue at chain start
    #    if looptypes[li] == "Nter" :
    #        refri = loops[li][1]
    #        if not chids.has_key(refri+1) or chids[refri] != chids[refri+1] : looptypes[li] = "NtoC"
    #        elif not chids.has_key(refri+2) or chids[refri] != chids[refri+2] : loops[li][1] += 1 ; looptypes[li] = "NtoC"
    #    if looptypes[li] == "Cter" :
    #        refri = loops[li][0]
    #        if not chids.has_key(refri-1) or chids[refri] != chids[refri-1] : looptypes[li] = "NtoC"
    #        elif not chids.has_key(refri-2) or chids[refri] != chids[refri-2] : loops[li][1] -= 1 ; looptypes[li] = "NtoC"
    #print loops, "\n", looptypes


    while 1 :
        canbreak = 1
        for li in range(len(loops)-1) : ## remove consecutive Cters or Nters with same chid
            if chids[loops[li][1]] != chids[loops[li+1][0]] : continue
            if looptypes[li] != looptypes[li+1] : continue
            if not looptypes[li] in ["Nter","Cter"] : continue
            loops[li][1] = loops[li+1][1]
            loops = loops[0:li+1] + loops[li+2:]
            looptypes = looptypes[0:li+1] + looptypes[li+2:]
            canbreak = None ; break
        if canbreak : break
#        print loops, "\n", looptypes
#    print loops, "\n", looptypes
    ## merge regions that are separated by 1 position in the chain
    while 1 :
        canbreak = 1
        for li in range(len(loops)-1) :
            if loops[li][1]+2 >= loops[li+1][0] and chids[ loops[li][1] ] == chids[ loops[li+1][0] ] :
                loops[li][1] = loops[li+1][1] ;
                if looptypes[li] == "Nter" and looptypes[li+1] == "Cter" : looptypes[li] = "NtoC"
                elif looptypes[li] == "Nter" and looptypes[li+1] == "loop" : looptypes[li] = "Nter"
                elif looptypes[li] == "loop" and looptypes[li+1] == "Cter" : looptypes[li] = "Cter"
                elif looptypes[li] == "Cter" or looptypes[li+1] == "Nter" : assert None
                loops = loops[0:li+1] + loops[li+2:]
                looptypes = looptypes[0:li+1] + looptypes[li+2:]
                canbreak = None ; break
        if canbreak : break


    order = []
    while len(order) < len(loops) :
        ri = int(math.floor( random.random() * len(loops) ))
        if not ri in order : order.append(ri)
    #print "order", order
    rloops, rlooptypes = [], []
    nrloops, nrlooptypes = [], []

    for oi in order :
        rloops.append( loops[oi] ) ;
        rlooptypes.append( looptypes[oi] )
        


    if len(loops) == 1  and modelN2C == 1:
        looptypes[0] = "NtoC"
        
    for i in range(len(loops)):
        if looptypes[i]  == 'loop' and buildType  != "poor" :
            nrloops.append([loops[i][0]+1,loops[i][1]-1])
            nrlooptypes.append(looptypes[i])
            print "2BL",resids[loops[i][0]+1] , resids[loops[i][1]-1] , looptypes[i] , buildType

        elif looptypes[i]  == 'cCter':
            nrloops.append([loops[i][0]+1,loops[i][1]])
            nrlooptypes.append(looptypes[i])
            print resids[loops[i][0]] , resids[loops[i][1]] , looptypes[i]

        elif looptypes[i]  == 'cNter':
            nrloops.append([loops[i][0],loops[i][1]])
            nrlooptypes.append(looptypes[i])
            print resids[loops[i][0]] , resids[loops[i][1]] , looptypes[i]            
            
        else :
            nrloops.append(loops[i])
            nrlooptypes.append(looptypes[i])





    if len(nrloops) == 0 :
        print "Cannot build.Please report problem with Error 210 to ak459@cam.ac.uk"
        import sys ; sys.exit()

    

    return nrloops, nrlooptypes

def locateRegionsRandomize(resids, chids, badresids,modelN2C=1,buildType="notpoor") :
    #    tempresids = {} ; tempchids = {}
    #    for k,v in 
    bad = [None] * len( resids.keys() )
    for ri in resids.keys() :
        if resids[ri] in badresids : bad[ri] = 1
    #print bad
    keys = list( resids.keys() )
    keys.sort()
    loops, looptypes = [], []
    loopOn = None
    debug = 1
    for k in keys :
        #print k, bad[k],
        if bad[k] and not loopOn :
            #print 'here 1'
            loopOn = 1
            loops.append( [k,k] )
            if k+1 == len(keys) or chids[k+1] != chids[k] : ## check for ending singleton
                looptypes.append("Cter") ; loopOn = None ; continue
            if k-1 < 0 or chids[k-1] != chids[k] :
                if debug == 1 : print "T1"
                looptypes.append("Nter")

            else : looptypes.append("loop")
        elif loopOn and not bad[k] :

            loops[ len(loops)-1 ][1] = k-1 ; loopOn = None
        elif loopOn and ( k+1==len(keys) or chids[k+1] != chids[k] ) :

            if looptypes[ len(looptypes)-1 ] == "Nter" : looptypes[ len(looptypes)-1 ] = "NtoC"
            else : looptypes[ len(looptypes)-1 ] = "Cter"
            loops[ len(loops)-1 ][1] = k ; loopOn = None
#    print loops, "\n", looptypes
    for li in range(len(loops)) : ## extend loops on both sides by 1 when loopsize is 1
        if looptypes[li] != "loop" : continue
        loops[li][0] -= 1 ; loops[li][1] += 1
        ## danger commented amk

    # for li in range(len(loops)) : ## loop -> ter for loops close to terminals
    #    if looptypes[li] != "loop" : continue
    #   for k in [1,2] :
    #      if (loops[li][0]-k < 0 or chids[loops[li][0]-k] != chids[loops[li][0]]) and looptypes[li] == "loop" :
    #                if debug == 1 : print "T2"
    #         loops[li][0] -= (k-1) ; looptypes[li] = "Nter"
    #    if (loops[li][1]+k >= len(resids) or chids[loops[li][1]+k] != chids[loops[li][1]]) and looptypes[li] == "loop" :
    #       loops[li][1] += (k-1) ; looptypes[li] = "Cter"
    #    print loops, "\n", looptypes

    #for li in range(len(loops)) : ## Cter -> NtoC if only 1 good residue at chain start
    #    if looptypes[li] == "Nter" :
    #      if debug == 1 : print "T3"
    #       refri = loops[li][1]
    #     if not chids.has_key(refri+1) or chids[refri] != chids[refri+1] : looptypes[li] = "NtoC"
    #    elif not chids.has_key(refri+2) or chids[refri] != chids[refri+2] : loops[li][1] += 1 ; looptypes[li] = "NtoC"
    #        if looptypes[li] == "Cter" :
    #           refri = loops[li][0]
    #          if not chids.has_key(refri-1) or chids[refri] != chids[refri-1] : looptypes[li] = "NtoC"
    #         elif not chids.has_key(refri-2) or chids[refri] != chids[refri-2] : loops[li][1] -= 1 ; looptypes[li] = "NtoC"
    #    print loops, "\n", looptypes
    while 1 :
        canbreak = 1
        for li in range(len(loops)-1) : ## remove consecutive Cters or Nters with same chid
            if chids[loops[li][1]] != chids[loops[li+1][0]] : continue
            if looptypes[li] != looptypes[li+1] : continue
            if not looptypes[li] in ["Nter","Cter"] : continue
            loops[li][1] = loops[li+1][1]
            loops = loops[0:li+1] + loops[li+2:]
            looptypes = looptypes[0:li+1] + looptypes[li+2:]
            canbreak = None ; break
        if canbreak : break
#        print loops, "\n", looptypes
#    print loops, "\n", looptypes
    ## merge regions that are separated by 1 position in the chain
    while 1 :
        canbreak = 1
        for li in range(len(loops)-1) :
            if loops[li][1]+2 >= loops[li+1][0] and chids[ loops[li][1] ] == chids[ loops[li+1][0] ] :
                loops[li][1] = loops[li+1][1] ;
                if looptypes[li] == "Nter" and looptypes[li+1] == "Cter" : looptypes[li] = "NtoC"
                elif looptypes[li] == "Nter" and looptypes[li+1] == "loop" : looptypes[li] = "Nter"
                elif looptypes[li] == "loop" and looptypes[li+1] == "Cter" : looptypes[li] = "Cter"
                elif looptypes[li] == "Cter" or looptypes[li+1] == "Nter" : assert None
                loops = loops[0:li+1] + loops[li+2:]
                looptypes = looptypes[0:li+1] + looptypes[li+2:]
                canbreak = None ; break
        if canbreak : break


    order = []
    while len(order) < len(loops) :
        ri = int(math.floor( random.random() * len(loops) ))
        if not ri in order : order.append(ri)
    #print "order", order
    rloops, rlooptypes = [], []
    nrloops, nrlooptypes = [], []

    for oi in order :
        rloops.append( loops[oi] ) ;
        rlooptypes.append( looptypes[oi] )
        


    if len(loops) == 1  and modelN2C == 1:
        looptypes[0] = "NtoC"
        
    for i in range(len(loops)):
        if looptypes[i]  == 'loop' and buildType  != "poor" :
            nrloops.append([loops[i][0]+1,loops[i][1]-1])
            nrlooptypes.append(looptypes[i])
            print "2BL",resids[loops[i][0]+1] , resids[loops[i][1]-1] , looptypes[i] , buildType

        elif looptypes[i]  == 'cCter':
            nrloops.append([loops[i][0]+1,loops[i][1]])
            nrlooptypes.append(looptypes[i])
            print resids[loops[i][0]] , resids[loops[i][1]] , looptypes[i]

        elif looptypes[i]  == 'cNter':
            nrloops.append([loops[i][0],loops[i][1]])
            nrlooptypes.append(looptypes[i])
            print resids[loops[i][0]] , resids[loops[i][1]] , looptypes[i]            
            
        else :
            nrloops.append(loops[i])
            nrlooptypes.append(looptypes[i])





    if len(nrloops) == 0 :
        print "Cannot build.Please report problem with Error 210 to ak459@cam.ac.uk"
        import sys ; sys.exit()

    

    return nrloops, nrlooptypes



def locateRegionsRandomize2(resids, chids, badresids,hetids = {}) :
    tempresids = {} ; tempchids = {}
    for k,v in resids.items():
        if hetids[k] == 1 :
            continue
        else:
            tempresids[k] = v
            tempchids[k] = chids[k]

    bad = [None] * len( tempresids.keys() )
    for ri in tempresids.keys() :
        if tempresids[ri] in badresids : bad[ri] = 1
    #print bad
    keys = list( tempresids.keys() )
    keys.sort()
    loops, looptypes = [], []
    loopOn = None
    for k in keys :
        #print k, bad[k],
        if bad[k] and not loopOn :
            #print 'here 1'
            loopOn = 1
            loops.append( [k,k] )
            if k+1 == len(keys) or tempchids[k+1] != tempchids[k] : ## check for ending singleton
                looptypes.append("Cter") ; loopOn = None ; continue
            if k-1 < 0 or tempchids[k-1] != tempchids[k] : looptypes.append("Nter")
            else : looptypes.append("loop")
        elif loopOn and not bad[k] :
            #print 'here 3'
            loops[ len(loops)-1 ][1] = k-1 ; loopOn = None
        elif loopOn and ( k+1==len(keys) or tempchids[k+1] != tempchids[k] ) :
            #print 'here 2'
            if looptypes[ len(looptypes)-1 ] == "Nter" : looptypes[ len(looptypes)-1 ] = "NtoC"
            else : looptypes[ len(looptypes)-1 ] = "Cter"
            loops[ len(loops)-1 ][1] = k ; loopOn = None
#    print loops, "\n", looptypes
    for li in range(len(loops)) : ## extend loops on both sides by 1 when loopsize is 1
        if looptypes[li] != "loop" : continue
        loops[li][0] -= 1 ; loops[li][1] += 1
#    print loops, "\n", looptypes
    for li in range(len(loops)) : ## loop -> ter for loops close to terminals
        if looptypes[li] != "loop" : continue
        for k in [1,2] :
            if (loops[li][0]-k < 0 or tempchids[loops[li][0]-k] != tempchids[loops[li][0]]) and looptypes[li] == "loop" :
                loops[li][0] -= (k-1) ; looptypes[li] = "Nter"
            if (loops[li][1]+k >= len(resids) or tempchids[loops[li][1]+k] != tempchids[loops[li][1]]) and looptypes[li] == "loop" :
                loops[li][1] += (k-1) ; looptypes[li] = "Cter"
#    print loops, "\n", looptypes
    for li in range(len(loops)) : ## Cter -> NtoC if only 1 good residue at chain start
        if looptypes[li] == "Nter" :
            refri = loops[li][1]
            if not tempchids.has_key(refri+1) or tempchids[refri] != tempchids[refri+1] : looptypes[li] = "NtoC"
            elif not tempchids.has_key(refri+2) or tempchids[refri] != tempchids[refri+2] : loops[li][1] += 1 ; looptypes[li] = "NtoC"
        if looptypes[li] == "Cter" :
            refri = loops[li][0]
            if not tempchids.has_key(refri-1) or tempchids[refri] != tempchids[refri-1] : looptypes[li] = "NtoC"
            elif not tempchids.has_key(refri-2) or tempchids[refri] != tempchids[refri-2] : loops[li][1] -= 1 ; looptypes[li] = "NtoC"
#    print loops, "\n", looptypes
    while 1 :
        canbreak = 1
        for li in range(len(loops)-1) : ## remove consecutive Cters or Nters with same chid
            if tempchids[loops[li][1]] != tempchids[loops[li+1][0]] : continue
            if looptypes[li] != looptypes[li+1] : continue
            if not looptypes[li] in ["Nter","Cter"] : continue
            loops[li][1] = loops[li+1][1]
            loops = loops[0:li+1] + loops[li+2:]
            looptypes = looptypes[0:li+1] + looptypes[li+2:]
            canbreak = None ; break
        if canbreak : break
##        print loops, "\n", looptypes
#    print loops, "\n", looptypes
    ## merge regions that are separated by 1 position in the chain
    while 1 :
        canbreak = 1
        for li in range(len(loops)-1) :
            if loops[li][1]+2 >= loops[li+1][0] and tempchids[ loops[li][1] ] == tempchids[ loops[li+1][0] ] :
                loops[li][1] = loops[li+1][1] ;
                if looptypes[li] == "Nter" and looptypes[li+1] == "Cter" : looptypes[li] = "NtoC"
                elif looptypes[li] == "Nter" and looptypes[li+1] == "loop" : looptypes[li] = "Nter"
                elif looptypes[li] == "loop" and looptypes[li+1] == "Cter" : looptypes[li] = "Cter"
                elif looptypes[li] == "Cter" or looptypes[li+1] == "Nter" : assert None
                loops = loops[0:li+1] + loops[li+2:]
                looptypes = looptypes[0:li+1] + looptypes[li+2:]
                canbreak = None ; break
        if canbreak : break
#        print loops, "\n", looptypes
    #sys.exit(0)
    order = []
    while len(order) < len(loops) :
        ri = int(math.floor( random.random() * len(loops) ))
        if not ri in order : order.append(ri)
    #print "order", order
    rloops, rlooptypes = [], []
    nrloops, nrlooptypes = [], []

    for oi in order :
        rloops.append( loops[oi] ) ;
        rlooptypes.append( looptypes[oi] )
        

    for i in range(len(loops)):
        if looptypes[i]  == 'loop':

            nrloops.append([loops[i][0]+1,loops[i][1]-1])
            
            nrlooptypes.append(looptypes[i])
            print resids[loops[i][0]]
        else :
            nrloops.append(loops[i])
            nrlooptypes.append(looptypes[i])

    import sys ; sys.exit()
    return nrloops, nrlooptypes


def get_and_verify_bootstrap_restraint(i , pts ,res , totp = 2):
    if totp == 1 :
        if ' CA ' in res[i].keys():
            ithCa = res[i][' CA ']
            pti = pts[ithCa]
            return True
        else :
            return False

    else :
        if ' CA ' in res[i].keys() and  ' CA ' in res[i+1].keys():
            
            return True
        else :
            return False



def makeKnownPos(numpts, blist) :
    knownPositions = range(numpts)
    for b in blist :
        for i in range(b.getOP().size()) : knownPositions.remove( b.getOP()[i] )
    return knownPositions



def break_chain_in_bands( start, stop, pts, min_band_size,max_band_size,res,resns,resids):
    
    last_band = []
    band_size = 0
    bands,bandtypes = [], []
    band_count = 0

    length = stop - start  + 1
    i = 0
    while (i <  length ) : 
        n_total_remaining = length - i
        n_remaining = n_total_remaining


        
        if (n_remaining == 0 and  n_total_remaining != 0) : 
            band_size = 0
            really_build = False

            
            while (band_size < n_total_remaining - 1 and
                   (get_and_verify_bootstrap_restraint(start + i + band_size - 2 , pts , res) == True  or
                    get_and_verify_bootstrap_restraint(start + i + band_size + 1 -2 , pts, res) == True) ) : 
                band_size = band_size + 1 
                
                
        else :

            band_size = random.randint( min_band_size , max_band_size)

            if ( n_remaining - band_size < min_band_size ):
                band_size = n_remaining


            really_build = True


            while (band_size > min_band_size and (i + band_size) < length - 1 and
                   (get_and_verify_bootstrap_restraint(start + i + band_size -2 , pts,res , 2) == False) ) : 
                band_size = band_size -1
            

            while (band_size < n_remaining and band_size < n_remaining - 1 and
                   (get_and_verify_bootstrap_restraint(start + i + band_size -2 , pts ,res ,2 ) == False )):

                band_size =  band_size  + 1


            if ( n_remaining - band_size < min_band_size  ):
                band_size = n_remaining

#        print get_and_verify_bootstrap_restraint(start + i + band_size -2 , pts ,res ,1 )
#        if (band_size < n_total_remaining - 1 and get_and_verify_bootstrap_restraint(start + i + band_size -2 , pts , res, 1)) :
#            
#            
#            print "No viable bootstrap found"
                   

#        print "  subband starting at %s of length %d into bands\n"%( resids[start + i], band_size)

        if band_count == 0 :
            band_type = "Nter"
            band_type = "NtoC"

        elif  i + band_size == length :
            band_type = "Cter"
        elif (band_count%2 == 0) or (band_count==1)  :
            band_type = "NtoC"
        else :
            band_type = "loop"
            band_type = "NtoC"
        if band_type == "Cter" :
            band_type = "NtoC"
            bands.append([start+i,band_size+start+i-1])

        else :
            bands.append([start+i,band_size+start+i-1])
        bandtypes.append(band_type)
        band_count = band_count + 1
        i = i + band_size

#    for r in range(len(bands)):
#        print "Banded",bands[r],bandtypes[r]
#        print "BANDIT",resids[bands[r][0]],resids[bands[r][1]],bandtypes[r]

    return bands,bandtypes





def bandBuild(loops, looptypes, bandInfrastructure, res, resids, resnums, resns, chids, inscodes, pts, extraRestraints, extraRestrOpt,        backtrack, popsize, radii, gridHelper, ranker,missingpts,natt) :

    if len(loops)  == 0 :
        print "Nothing to build..exiting now"
        import sys ; sys.exit()



    print "Band building starts ........\n\n" 
    unbuiltpts= []
    for bandIndex in range(len(loops)):
        print "BAND %d *** [%s] [%s] of type %s" % (bandIndex, resids[ loops[bandIndex][0] ] , resids[ loops[bandIndex][1] ], looptypes[bandIndex])

    bandinds = range( len(loops) )
    built = [] ; unbuilt = []
    blist = []
    for bl,nt,bf,rl,reorderBuilders in bandInfrastructure :
        for b in bl : blist.append(b)
    knownPositions = makeKnownPos(len(pts), blist) ## knownPositions are those that are not built by builders
    debug = 1
    if debug == 1 :
        pts2res = {}
        for a1,b1 in res.items():
            for a2,b2 in b1.items():
                pts2res[b2] = resnums[a1]+a2
    #    print "known",knownPositions
    #    for kp in knownPositions :
    #        print "HP",pts2res[kp],

    while( len(bandinds) > 0 ) :
        bandIndex = bandinds[ int( math.floor(random.random() * len(bandinds)) ) ]
        print "BAND***(%d remain)********** [%s] [%s] %s" % (len(bandinds), resids[ loops[bandIndex][0] ] , resids[ loops[bandIndex][1] ], looptypes[bandIndex])
        bl,nt,bf,rl, reorderBuilders = bandInfrastructure[bandIndex]
        allbop = []
        restrs = [] ## add extra restraints
        for xr in extraRestraints : restrs.append(xr)
        for r in rl : restrs.append(r)
        if bandIndex < 0 : break
        buffRenderer = BufferModelRenderer(pts)
        
        strategy = PopulationStrategy(backtrack, popsize, pts, radii, knownPositions, gridHelper, bl, nt, bf, reorderBuilders, restrs, extraRestrOpt, buffRenderer, natt, 1, None)
        #strategy.snapRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [], "snap.pdb")
        strategy.ranker = ranker

        if strategy.execute() == 0 : ## have to restart, restore knownPositions OR just let the band be same as given
            #bandinds = range(len(loops)) ; knownPositions = makeKnownPos(len(pts), blist)
            #print "BANDBUILDING RESTARTS" ; continue
            print "BAND CANT BE BUILT! leaving it as it was"
            unbuilt.append([resids[ loops[bandIndex][0] ] , resids[ loops[bandIndex][1] ], looptypes[bandIndex] ])
            print "missin",missingpts
            #import sys; sys.exit()
            for b in bl : ## modify knownPositions and pts
                for i in range(b.getOP().size()) :
                    if  b.getOP()[i]  not in missingpts:
                        knownPositions.append( b.getOP()[i] )
                    unbuiltpts.append( b.getOP()[i] )
        else:
            built.append( [resids[ loops[bandIndex] [0] ] , resids[ loops[bandIndex][1] ], looptypes[bandIndex] ] )
            for b in bl : ## modify knownPositions and pts
                for i in range(b.getOP().size()) : knownPositions.append( b.getOP()[i] )
        #ModelRenderer(res, resns, chids, resnums, inscodes, [], "band%d.pdb"%bandIndex).render(pts)
        #print "wrote band", "band%d.pdb"%bandIndex
        bandinds.remove(bandIndex)


    return built, unbuilt,unbuiltpts

class Multiloop :
    def __init__(s, pdbfile, badresids, mconly, caRad, scRad, scReduction, guidedSampling, popsize, backtrack, nmodels, outpdb, restrGen, prepC,sconly=None, cacaCutoff = 5.,null=[],modelN2C=0,natt=5,missingpts = {},loopclosure = 1 , allOpt=0, edOpt=0):

        s.pdbfile, s.badresids, s.mconly, s.caRad, s.scRad, s.scReduction, s.guidedSampling, s.popsize, s.backtrack, s.nmodels, s.outpdb, s.restrGen \
            = pdbfile, badresids, mconly, caRad, scRad, scReduction, guidedSampling, popsize, backtrack, nmodels, outpdb, restrGen
        if prepC == None :
            s.prepC = PrepareChain("PRL")
        else :
            s.prepC = prepC

        s.natt = natt

        s.sconly = sconly
        if not mconly in [1,None] : print mconly ; assert None
        s.ranker = None
        s.cellsym = None
        s.cacaCutoff  = cacaCutoff
        s.null = null
        s.modelN2C = modelN2C
        s.missingpts = missingpts
        
        s.poor = "notpoor"
        s.loopstartnum  = None ; s.loopstopnum  = None
        s.loopclosure = loopclosure
        s.allopt = allOpt
        s.edopt  = edOpt

    # given the resids to be rebuilt, multiloop will identify the bands or bands can be provided directly
    # multiloop will then create necessary builders/restraints etc and build

    
    def run(s) :

        from prepareChainV5 import buildAcb
        import checkProtChainsV4

        print "Reading coordinate file", s.pdbfile
        prot = protein(s.pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

        unbuiltpts = [] ; guidedSamplingRadius = None
        if s.guidedSampling != None : guidedSamplingRadius = s.caRad
        blist, numtrials, bfoll, rlist, dummies = [], [], {}, [], {} ; bandInfrastructure = []
                

        if type(s.badresids) == type([1,2]) and len(s.badresids)==2 and type(s.badresids[0]) == type([1,2]) : # bands are given 
            loops, looptypes = s.badresids[0], s.badresids[1]
        else :
            if s.badresids == None :
                s.badresids =  [];             badids = []
                from data import three2one
                for k , v in resids.items():
                    if resns[k] in three2one.keys() :
                        s.badresids.append(resids[k])
                        badids.append(k)
                    else :
                        print "WARNING!! %s not recognised as standard amino acid"%resns[k]

                print "AA RESIDUES TO BE REBUILT ARE",
                if len(s.badresids) ==0:
                    print "NONE"
                    import sys; sys.exit()

                for bd in s.badresids :
                    print
                    print "          ========>",bd

            else :
                for bd in s.badresids :
                    print
                    print "          ========>",bd

                badids = []
                for k ,v in resids.items():
                    if v in s.badresids :
                        badids.append(k)


            if s.mconly :
                res, pts = removeSC(res, pts,res.keys(),badids)


            print "\nGrouping residues into bands...."
            loops, looptypes = locateRegionsRandomize3(resids, chids, s.badresids,s.modelN2C,s.poor) # locate loops and order them randomly
            print loops, looptypes

        for k,v in res.items():
            print resids[k],k , v
#        import sys ; sys.exit()            
        print "\n\nChecking for missing atoms and chain breaks....."
        print "\n\n=====================  MISSING ATOMS SUMMARY ====================================================="
        print
        mcmiss, scmiss, chainBreaks = checkProtChainsV4.check(res, resids, resnums, resns, chids, inscodes, pts, s.cacaCutoff,s.null,1,s.mconly)

        if s.caRad != None :
            bands = checkProtChainsV4.getbands(res, resids, resnums, resns, chids, inscodes, pts, s.cacaCutoff,s.null,1,s.mconly)

        print "\n\n================== END MISSING ATOMS SUMMARY ====================================================="
        print
        print

        
        #for b in bands :
         #   print resids[bands[b][0]] , resids[bands[b][1]]
#        import sys; sys.exit()

        chainBreaks = list(chainBreaks) ;         chainBreaks.sort() ;
        wholenewloops = [];        wholenewlooptypes = []
        modloops , modlooptypes ,  notbuilding = [], [], []
        
        for k in range(len(loops)):
            
            loopstart = loops[k][0] ;             loopend = loops[k][1]

            print "bbl",resids[loopstart],resids[loopend],looptypes[k]
            if s.caRad != None :
                for dol in range(loops[k][0],loops[k][1]+1):
                    if ' CA ' not in res[dol].keys():
                        print "If using Ca-restraints, all Calpha atoms need to be present"


                for bb in bands.values() :
                    bstart = bb[0] ; bstop = bb[1]
                    
                    print "Going through band",resids[bb[0]],resids[bb[1]] , bb[0] , bb[1]
                    

                    if s.caRad != None :
                        if loopstart == bstart and (bstop - bstart) < 2 and loopend >= bstop :
                            for xx in range(loops[k][0] ,bstop + 1) : 
                                notbuilding.append(resids[xx])
                                print "1 Notbuilding ... ",resids[xx]
                        elif loopstart == bstart and (bstop - bstart) >= 2 and loopend >= bstop :

                            wholenewlooptypes.append("NtoC")
                            wholenewloops.append([bstart, bstop])
                            print "A1"
                            


                        elif loopstart == bstart and  loopend < bstop and  (bstop - bstart > 0):
                            wholenewlooptypes.append("Nter")
                            wholenewloops.append([bstart, loopend])                            
                            print "A2"
                        elif loopstart == bstart and  loopend < bstop and  (bstop - bstart == 0):
                            for xx in range(loopstart , loopend+1) : 
                                notbuilding.append(resids[xx])
                                print "2 Notbuilding ...",resids[xx]
                            



                        elif loopstart > bstart  and loopend >= bstop and loopstart <= bstop and (bstop - bstart > 0):
                            wholenewlooptypes.append("Cter")
                            wholenewloops.append([loopstart,bstop])                            
                            print "A3"

                        elif loopstart > bstart  and loopend >= bstop and loopstart <= bstop and (bstop - bstart == 0):
                            for xx in range(loopstart , bstop+1) : 
                                notbuilding.append(resids[xx])
                                print "3 Notbuilding ...",resids[xx]




                        
                        elif loopstart > bstart and loopend < bstop:
                            wholenewlooptypes.append("loop")
                            wholenewloops.append([loopstart,loopend])
                            print "A4"



                        elif loopstart < bstart and loopend == bstart and (bstop - bstart > 0) :
                            wholenewlooptypes.append("Nter")
                            wholenewloops.append([bstart,bstart])
                            print "A5"
                            
                        elif loopstart < bstart and loopend == bstart and (bstop - bstart == 0) :
                            for xx in range(bstart , bstart+1) : 
                                notbuilding.append(resids[xx])
                                print " 4 Notbuilding ...",resids[xx]



                        elif loopstart < bstart and loopend < bstop and loopend in range(bstart,bstop) and (bstop - bstart > 0):
                            wholenewlooptypes.append("Nter")
                            wholenewloops.append([bstart,loopend])
                            print "A6"
                            
                        elif loopstart < bstart and loopend < bstop and loopend in range(bstart,bstop) and (bstop - bstart == 0):
                            for xx in range(bstart , loopend+1) : 
                                notbuilding.append(resids[xx])
                                print " 5 Notbuilding ...",resids[xx]


                            
                        elif loopstart < bstart and loopend >= bstop and bstop-bstart >= 2:
                            wholenewlooptypes.append("NtoC")
                            wholenewloops.append([bstart,bstop])
                            print "A7"
                        elif loopstart < bstart and loopend >= bstop and bstop-bstart < 2:
                            for xx in range(bstart ,bstop + 1) : 
                                notbuilding.append(resids[xx]);
                                print " 6 Notbuilding ...",resids[xx]
                                
                        else :
                            print "Skipping ... "
                            
                        
                    print "--------------------------------------------------------------------------------------"
                    print
                    print

            else :
                if looptypes[k] == "loop" and len(loops) == 1 and s.poor != "poor":


                    nanchor = int(s.loopstartnum) - 1
                    canchor = int(s.loopstopnum) + 1

                    nanchor1 = int(s.loopstartnum) - 2
                    canchor1 = int(s.loopstopnum) + 2

                    print nanchor,canchor,nanchor1,canchor1 

                    nfound , cfound  = None , None
                    nfound1 , cfound1  = None , None

                    loopchain = []
                    lchid = chids[loops[k][0]]
                    for k1 , v1 in chids.items():
                        if v1  ==  lchid : 
                            loopchain.append(k1)
                    for k1 in loopchain :
                        if int(resnums[k1]) == nanchor :
                            if ' CA ' in res[k1].keys() :
                                nfound = 1
                        if int(resnums[k1]) == canchor :
                            if ' CA ' in res[k1].keys() :
                                cfound = 1                                
                        if int(resnums[k1]) == canchor1 :
                            cfound1  = 1

                        if int(resnums[k1]) == nanchor1 :
                            nfound1  = 1                            

                    if nfound != None and cfound != None  and nfound1 != None and cfound1 != None  :
                        looptypes[k] = "loop"

                    elif (nfound == None or nfound1 == None) and cfound != None and cfound1 !=None :
                        print "WARNING !! N-anchor residue not found"
                        print "Loop closure restraints at N-term will not be used"
                        looptypes[k] = "Nter"    

                    elif nfound != None and (cfound == None or cfound1 == None) and nfound1 != None :
                        print "WARNING !! C-anchor residue not found"
                        print "Loop closure restraints at C-term will not be used"
                        looptypes[k] = "Cter"

                    else :
                        
                        if nfound == None and  nfound1 == None :
                            print "Anchor residues: Nterm %d [%s] (CA-atom) and %d [%s] not found in coordinates "%(nanchor,lchid,nanchor1,lchid)
                            print "Cannot build loop"
                            print "Exiting now"
                            import sys ; sys.exit()

                        if cfound == None and  cfound1 == None :
                            print "Anchor residues: Cterm %d [%s] (CA-atom) and %d [%s] not found in coordinates "%(canchor,lchid,canchor1,lchid)
                            print "Cannot build loop"
                            print "Exiting now"
                            import sys ; sys.exit()

                        else :
                            print "Anchor residues not found in coordinates "
                            print "Cannot build loop"
                            print "Exiting now"
                            import sys ; sys.exit()                            

        if s.caRad != None : 
            loops = wholenewloops
            looptypes = wholenewlooptypes

        print "dfffsdasf",
        print loops,looptypes

        min_band_size = 5 ; max_band_size = 500
        for kx in range(len(loops)):
            if loops[kx][1] - loops[kx][0] > 700  and (looptypes[kx] == "NtoC" or looptypes[kx]== 'Cter' or looptypes[kx] == 'Nter') and s.caRad != None :
                modloops,modlooptypes = break_chain_in_bands( loops[kx][0], loops[kx][1], pts, min_band_size ,max_band_size ,res,resns,resids)
                loops.remove(loops[kx])
                looptypes.remove(looptypes[kx])
                for ml in modloops:
                    loops.append(ml)
                for mlt in modlooptypes :
                    looptypes.append(mlt)
            else :
                if loops[kx][1] - loops[kx][0] > 700  and s.caRad == None and looptypes[kx] != "loop":
                    print "Rappertk cannot build a segment greater than 700 residues"
                    print "Either (1) set --use-ca-restraints True "
                    print "                  or                    "
                    print "       (2) Restrict residue range to < 700 residues"
                    
                    import sys ; sys.exit()
                    
        print "loops",loops
        print "looptypes",looptypes
        

        mcMissInds = []; scMissInds = [] ; caMissInds = []; s.badresids = []  ; badids = []
        mcmiss, scmiss, chainBreaks1 = checkProtChainsV4.check(res, resids, resnums, resns, chids, inscodes, pts, s.cacaCutoff,[],0)
        
        for startindex,endindex in loops :
            for ri in range(startindex, endindex+1) :
                s.badresids.append( resids[ri] )
                badids.append(ri)

################# Addmissing atoms ###################################################                

        for bd in badids :
            if resids[bd] in notbuilding:
                print "Error in banding", notbuilding,resids[bd]
                import sys ; sys.exit()

        if s.mconly != 1:
            for sc in scmiss :
                if sc in badids : 
                    scMissInds.append(sc)
                
        for mc in mcmiss :
            if mc in badids : 
                if ' CA ' not in res[mc].keys() :
                    caMissInds.append(mc)

                mcMissInds.append(mc)
                misspts = incompleteMCcorrection(res, resns, pts , mc,resids)
                for km,vm in misspts.items() :
                    if km in s.missingpts.keys() :
                        for vvm in vm : 
                            s.missingpts[km].append(vvm)
                    else :
                        s.missingpts[km] = []
                        for vvm in vm : 
                            s.missingpts[km].append(vvm)



                        
        if s.mconly != 1:
            for ri in scMissInds :
                misspts = incompleteSCcorrection2(res, resns, pts ,ri,resids)
                for km,vm in misspts.items() :
                    if km in s.missingpts.keys() :
                        for vvm in vm : 
                            s.missingpts[km].append(vvm)
                    else :
                        s.missingpts[km] = []
                        for vvm in vm : 
                            s.missingpts[km].append(vvm)
                            

                if ri not in mcMissInds : 
                    print "Building temporary sidechain for",resids[sc]
                    buildAcb(res, resids, resnums, resns, chids, inscodes, pts, ri)


#            vvpts = VecVecFloat(pts)
#            for ri in scMissInds :
#                s.prepC.makeChiBuilder(ri,ri,ri, res,resns,resids).buildSample(vvpts, 0)
#            for pi in range(len(pts)) : pts[pi] = [ vvpts[pi][0], vvpts[pi][1], vvpts[pi][2], ]

            for r in range(len(loops)):
                if looptypes[r] == "Nter":
                    for sc in scmiss :
                        if sc == int(loops[r][1]) + 1  :
                            misspts = incompleteSCcorrection2(res, resns, pts ,sc,resids)
                            scMissInds.append(sc)
                            for km,vm in misspts.items() :
                                if km in s.missingpts.keys() :
                                    for vvm in vm : 
                                        s.missingpts[km].append(vvm)
                                else :
                                    s.missingpts[km] = []
                                    for vvm in vm : 
                                        s.missingpts[km].append(vvm)
                                        
                                
                            if sc not in mcmiss : 
                                print "Building temporary sidechain for",resids[sc] #,sc,mcMissInds
                                buildAcb(res, resids, resnums, resns, chids, inscodes, pts, sc)

#                            vvpts = VecVecFloat(pts)
#                            s.prepC.makeChiBuilder(sc,sc,sc, res,resns,resids).buildSample(vvpts, 0)
#                            for pi in range(len(pts)) :pts[pi] = [ vvpts[pi][0], vvpts[pi][1], vvpts[pi][2], ]
#                            print "Addining missing sidechain atoms:",sc
                            
                if looptypes[r] == "Cter":
                    for sc in scmiss :
                        if sc == int(loops[r][0]) - 1  :
                            misspts = incompleteSCcorrection2(res, resns, pts ,sc,resids)
                            scMissInds.append(sc)
                            for km,vm in misspts.items() :
                                if km in s.missingpts.keys() :
                                    for vvm in vm : 
                                        s.missingpts[km].append(vvm)
                                else :
                                    s.missingpts[km] = []
                                    for vvm in vm : 
                                        s.missingpts[km].append(vvm)

                                        
                            if sc not in mcmiss: 
                                print "Building temporary sidechain for",resids[sc]
                                buildAcb(res, resids, resnums, resns, chids, inscodes, pts, sc)
                                #                            vvpts = VecVecFloat(pts)                                
                                #                          s.prepC.makeChiBuilder(sc,sc,sc, res,resns,resids).buildSample(vvpts, 0)
                                #                           for pi in range(len(pts)) : pts[pi] = [ vvpts[pi][0], vvpts[pi][1], vvpts[pi][2], ]

                if looptypes[r] == "loop":
                    for sc in scmiss :
                        if sc == int(loops[r][0]) - 1  :
                            misspts = incompleteSCcorrection2(res, resns, pts ,sc,resids)
                            scMissInds.append(sc)
                            for km,vm in misspts.items() :
                                if km in s.missingpts.keys() :
                                    for vvm in vm : 
                                        s.missingpts[km].append(vvm)
                                else :
                                    s.missingpts[km] = []
                                    for vvm in vm : 
                                        s.missingpts[km].append(vvm)

                                        
                            if sc not in mcmiss : 
                                print "Building temporary sidechain for",resids[sc]
                                buildAcb(res, resids, resnums, resns, chids, inscodes, pts, sc)

#                            vvpts = VecVecFloat(pts)
#                            s.prepC.makeChiBuilder(sc,sc,sc, res,resns,resids).buildSample(vvpts, 0)
#                            for pi in range(len(pts)) :    pts[pi] = [ vvpts[pi][0], vvpts[pi][1], vvpts[pi][2], ]

                        if sc == int(loops[r][1]) + 1  :
                            misspts = incompleteSCcorrection2(res, resns, pts ,sc,resids) 
                            scMissInds.append(sc)
                            for km,vm in misspts.items() :
                                if km in s.missingpts.keys() :
                                    for vvm in vm : 
                                        s.missingpts[km].append(vvm)
                                else :
                                    s.missingpts[km] = []
                                    for vvm in vm : 
                                        s.missingpts[km].append(vvm)

                            if sc not in mcmiss : 
                                print "Building temporary sidechain for",resids[sc]
                                buildAcb(res, resids, resnums, resns, chids, inscodes, pts, sc)

                            ###vvpts = VecVecFloat(pts)
                            #s.prepC.makeChiBuilder(sc,sc,sc, res,resns,resids).buildSample(vvpts, 0)
                            #for pi in range(len(pts)) : pts[pi] = [ vvpts[pi][0], vvpts[pi][1], vvpts[pi][2], ]                            
                            
                            
                            ## print to check: What abt restraints






################# End - - - - Add missing atoms ###################################################                

        rev_resids = {} ; missingpointids = []
        for king, queen in resids.items():
            rev_resids[queen] = king


        keyt = s.missingpts.keys()
        print keyt
        for xd in keyt :
            missid = rev_resids[xd]
            for aid in s.missingpts[xd]:
                print resids[missid]
                ptid  = res[missid][aid]
                missingpointids.append(ptid)
                
            


        for li in range(len(loops)) :
            print "BANDS ARE ",resids[loops[li][0]],resids[loops[li][1]],looptypes[li]
            startindex, endindex = loops[li]
            chid = chids[startindex]
            assert chid == chids[endindex]
            bl,bf,nt,rl,dum, reorderBuilders = None,None,None,None,None, 1

            
            if s.caRad == None : 
                for x in range(startindex,endindex+1):
                    caMissInds.append(x)

                    

            if s.scRad == None:
                for x in range(startindex,endindex+1):
                    scMissInds.append(x)
                s.scRad = 1.


            if s.sconly == 1 :
                if looptypes[li] == "NtoC" :
                    for rr in res[startindex].values():
                        unbuiltpts.append(rr)
                    for rr in res[endindex].values():
                        unbuiltpts.append(rr)                        

                    bl,bf,nt,rl,dum = s.prepC.prepareSConly(startindex+1, endindex-1,  res, resids, resnums, resns, chids, inscodes, pts, s.scRad, scMissInds,caMissInds,mcMissInds)

                elif looptypes[li] == "Nter" :
                    for rr in res[startindex].values():
                        unbuiltpts.append(rr)
                        
                    bl,bf,nt,rl,dum = s.prepC.prepareSConly(startindex+1, endindex,  res, resids, resnums, resns, chids, inscodes, pts, s.scRad, scMissInds,caMissInds,mcMissInds)
                elif looptypes[li] == "Cter" :
                    for rr in res[endindex].values():
                        unbuiltpts.append(rr)
                    bl,bf,nt,rl,dum = s.prepC.prepareSConly(startindex, endindex-1,  res, resids, resnums, resns, chids, inscodes, pts, s.scRad, scMissInds,caMissInds,mcMissInds)
                elif looptypes[li] == "loop" :
                    bl,bf,nt,rl,dum = s.prepC.prepareSConly(startindex, endindex,  res, resids, resnums, resns, chids, inscodes, pts, s.scRad, scMissInds,caMissInds,mcMissInds)

                    

            elif looptypes[li] == "NtoC" :
                if resns[ loops[li][0] ] in [ "  A", "  T", "  C", "  G", "  U", ] :
                    bl,nt,bf,rl,dum = prepareRNAchain(chid, res, resids, resnums, resns, chids, inscodes, pts, s.caRad, s.scRad)
                else :
            
                    rev_resnum = {} ; chainkey = []
                    chainid = chids[startindex]
                    for k,v in chids.items():
                        if chids[k] == chainid:
                            chainkey.append(k)
                    for k in chainkey :
                        rev_resnum[int(resnums[k])] = int(k)

                    startResidue = int(resnums[startindex]) ; endResidue = int(resnums[endindex])

                    if (startindex not in resnums.keys()) or resids[startindex] in s.null :
                        print "N-ter residue [%s] for building chain needs to be present in input coordinate file"%resids[startindex]
                        import sys  ; sys.exit()

                    if ' CA ' not in res[startindex].keys():
                        print "C-alpha atom of anchor N residue [%s] for chain to be modelled needs to be present in the pdbfile"%resids[startindex]
                        import sys  ; sys.exit()


                    #if ((startResidue+1) not in rev_resnum.keys()) or resids[pstartindex] in s.null :
                    #    print "Anchor N residue [%s] for chain to be modelled needs to be present in the pdbfile"%startResidue
                    #    import sys  ; sys.exit()

                    #if ' CA ' not in res[pstartindex].keys():
                    #    print "C-alpha atom of anchor N residue [%s] for chain to be modelled needs to be present in the pdbfile"%resids[pstartindex]
                    #    import sys  ; sys.exit()
                        
                    

                    bl,nt,bf,rl,dum = s.prepC.preparePeptideChain(startindex,endindex,chid, res, resids, resnums, resns, chids, inscodes, pts, s.mconly, 1, guidedSamplingRadius, s.caRad, s.scRad, scMissInds,caMissInds)

            elif looptypes[li] == "Nter" :

                rev_resnum = {} ; chainkey = []
                chainid = chids[startindex]
                for k,v in chids.items():
                    if chids[k] == chainid:
                        chainkey.append(k)
                for k in chainkey :
                    rev_resnum[int(resnums[k])] = int(k)

                #print startindex,endindex
                #import sys ; sys.exit()
                if((int(resnums[endindex])+1) not in rev_resnum.keys()):
                    print " Anchor residue %s  needs to be present in the coordinate file"%resids[endindex+1]
                    import sys ; sys.exit()


                resNanchorindex = endindex + 1

                if " CA " not in res[resNanchorindex].keys() :
                    print " C-alpha atom of anchor residue %s  needs to be present in the coordinate file"%resids[endindex+1]
                    
                    import sys ;  sys.exit()
                    
                dum, firstindex = addNdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts)


                bl,bf,nt,rl = s.prepC.prepareChainTerminal("Nterm", startindex, endindex, firstindex, res, resids, resnums, resns, chids, inscodes, pts, s.mconly, s.caRad, s.scRad, guidedSamplingRadius, scMissInds,caMissInds)


            elif looptypes[li] == "Cter" :
                rev_resnum = {} ; chainkey = []
                chainid = chids[startindex]
                for k,v in chids.items():
                    if chids[k] == chainid:
                        chainkey.append(k)
                for k in chainkey :
                    rev_resnum[int(resnums[k])] = int(k)

                if int(resnums[startindex-1]) not in rev_resnum.keys():
                    print "Anchor residue %s for C-terminal to be modelled needs to be present in the pdbfile"%resids[startindex-1]
                    import sys ;  sys.exit()

                resCanchorindex = startindex-1

                if " CA " not in res[resCanchorindex].keys():
                    print "CA atom of Anchor  C residue  for C-terminal to be modelled needs to be present in the pdbfile"
                    import sys ;  sys.exit()

                    
                dum, lastindex = addCdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts)
                bl,bf,nt,rl = s.prepC.prepareChainTerminal("Cterm", startindex, endindex, lastindex, res, resids, resnums, resns, chids, inscodes, pts, s.mconly, s.caRad, s.scRad, guidedSamplingRadius, scMissInds,caMissInds)



            elif looptypes[li] == "loop" :
                rev_resnum = {} ; chainkey = []
                chainid = chids[startindex]
                for k,v in chids.items():
                    if chids[k] == chainid:
                        chainkey.append(k)
                for k in chainkey :
                    rev_resnum[int(resnums[k])] = int(k)
                    
                if((startindex-1) not in resnums.keys()):
                    print "Anchor N residue for loop to be modelled needs to be present in the pdbfile"
                    import sys ;  sys.exit()

                if((startindex-2) not in resnums.keys()):
                    print "Anchor N-terminal residues,start-1 and start-2,  for loop to be modelled needs to be present in the pdbfile"
                    import sys ;  sys.exit()                    
        
                if((endindex + 1 ) not in resnums.keys()):
                    print "Anchor C-terminal residue for loop to be modelled needs to be present in the pdbfile"
                    import sys ; sys.exit()

                if((endindex + 2 ) not in resnums.keys()):
                    print "Anchor C-terminal residue for loop to be modelled needs to be present in the pdbfile"
                    import sys ; sys.exit()                    
                    
                resNanchorindex = startindex-1
                resCanchorindex = endindex+1

                if " CA " not in res[resNanchorindex].keys() or  " CA " not in res[resCanchorindex].keys():
                    print "CA atom of Anchor N + C residue  for loop to be modelled needs to be present in the pdbfile"
                    import sys ; sys.exit()

                
#                if endindex-startindex < 5 or random.random() > 1.5 :
                 if endindex-startindex <= 12  : #or random.random() > 1.5 :

                    bl,bf,nt,rl = s.prepC.preparePeptideLoop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, s.mconly, s.caRad, s.scRad, scMissInds, guidedSamplingRadius,None,caMissInds,s.loopclosure)


                else :
                    
                    bl,bf,nt,rl = s.prepC.preparePeptideMidloop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, s.mconly, s.caRad, s.scRad, scMissInds, guidedSamplingRadius,caMissInds,s.loopclosure)
                    
                    reorderBuilders = None
            

            else :
                print "Dont understand looptype", looptypes[li] ; assert None


            print loops[li],looptypes[li]
            print bl
            print dummies
            print dum
            print "================================"
            mergeBRlists(blist, bfoll, numtrials, rlist, bl, bf, nt, rl, dummies, dum)
            #bandInfrastructure.append( (bl,nt,bf,[], reorderBuilders) )
            bandInfrastructure.append( (bl,nt,bf,rl, reorderBuilders) )

        #import sys ; sys.exit()

        ai2pos, aiSC = [-999] * len(pts), [-999] * len(pts) ## flag atoms
        for index,val in res.items() :
            for name, ai in val.items() :
                ai2pos[ai] = index
                if name == ' N  ' : aiSC[ai] = 0
                elif name == ' CA ' : aiSC[ai] = 1
                elif name == ' C  ' : aiSC[ai] = 2
                elif name == ' O  ' : aiSC[ai] = 3
                elif name == ' SG ' and resns[index] == 'CYS' : aiSC[ai] = 5 # for disulphide vdw reduction
                elif name == ' CD ' and resns[index] == 'PRO' : aiSC[ai] = 6 # PRO CD shdnt be clash-checked with prev res's C
                else : aiSC[ai] = 4
        for k in dummies.keys() :
            for an in dummies[k] :
                aiSC[res[k][an]] = -1
#                print "makind fummy", resids[k],an
                
        radii = [0] * len(pts) ## make radii for all atoms including dummie
        for index, val in res.items() :
            for aname, pi in val.items() :
                if resns[index] in vdwr.keys() :
                    if aname not in  vdwr[resns[index]].keys() :
                        print "UNKNOWN atom types [%s] %s"%(aname , resids[index])
                        import sys ; sys.exit()
                    radii[pi] = vdwr[resns[index]][aname]
                else :
                    try : radii[pi] = vdwr['XXX'][aname[1]]
                    except KeyError : radii[pi] = vdwr['XXX']['C']
        if len(radii) != len(pts) :
            print "Van der Waals radii for certain atoms could not be found"
            import sys ; sys.exit()

        gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), s.scReduction, data.PROBE_DISULFIDE_OVERLAP_MARGIN)
        if s.cellsym != None : gridHelper.setCellsym(s.cellsym[0], s.cellsym[1], s.cellsym[2], s.cellsym[3], s.cellsym[4], s.cellsym[5], s.cellsym[6])

        extraRestraints, extraRestrOpt = [], []
        if not type(s.restrGen) == type([1,2]) : s.restrGen = [s.restrGen]
        for rg in s.restrGen :
            if rg == None : continue
            xrlist, optional = rg.generate(blist, aiSC, res, resids, pts)
            rsize = len(extraRestraints)
            for xr in xrlist : extraRestraints.append(xr)
            for ori in optional : extraRestrOpt.append(rsize + ori)

#        if allopt == 1 : 
#        for r in rlist :
            
        #print extraRestraints
        #print extraRestrOpt
        #import sys ;
        #sys.exit()
        modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, dummies.keys(), s.outpdb)
        nbuilt = 0

        all_built = {};        all_unbuilt  = {} ; #unbuiltpts = []
        for nm in range(s.nmodels) :
            all_built[nm] = [] ;              all_unbuilt[nm] = []
            built, unbuilt,unbuiltpts1 = bandBuild(loops, looptypes, bandInfrastructure, res, resids, resnums, resns, chids, inscodes, pts, extraRestraints, extraRestrOpt,    s.backtrack, s.popsize, radii, gridHelper, s.ranker,missingpointids,s.natt)
            for r in unbuiltpts1 :
               unbuiltpts.append(r)
            print "Missingptsids",missingpointids
            print "s.miss",s.missingpts

            modelRenderer.render(pts,None,missingpointids,unbuiltpts) ; nbuilt += 1
            if len(unbuiltpts1) != 0 :
                print "BANDI failed",s.outpdb
            else : 
                print "BANDI succ",s.outpdb                
            all_built[nm].append(built);        all_unbuilt[nm].append(unbuilt)

          
        return nbuilt,all_built,all_unbuilt








def main(pdbfile, startResid, endResid, mconly, caRad, scRad, scReduction, guidedSampling, popsize, backtrack, nmodels, outpdb, restrGen=None) :
    prot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=0)

    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)

    startindex, endindex = None, None
    for index in range(len(resids)) :
        print "----%s-----%s------" % (resids[index], startResid)
        if resids[index] == startResid : startindex = index
        if resids[index] == endResid : endindex = index

    assert startindex <= endindex and startindex != None and endindex != None
    assert chids[startindex] == chids[endindex] == chids[startindex-2] == chids[endindex+2]

    if mconly : res, pts = removeSC(res, pts, res.keys())

    knownPositions = range(len(pts))
    for ind in range(startindex-1,endindex+2) :
        for an,pi in res[ind].items() :
            knownPositions.remove(pi)
    knownPositions.append(res[startindex-1][' N  '])
    knownPositions.append(res[startindex-1][' CA '])
    knownPositions.append(res[endindex+1][' CA '])
    knownPositions.append(res[endindex+1][' C  '])
    knownPositions.append(res[endindex+1][' O  '])

    guidedRadius = None
    if guidedSampling : guidedRadius = caRad
    blist, bfoll, numtrials, irlist = preparePeptideLoop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, [], guidedRadius)

    print "Residues and atom-numbers-----------------------"
    keys = res.keys() ; keys.sort()
    for k in keys : print resids[k], res[k]
    print "------------------------------------------------\n\n\n"

    for r in irlist : print r.name()
    for b in blist : b.describe() ; print ''
    print "------------------------------------------------\n\n\n"

    radii = [0] * len(pts) ## make radii for all atoms including dummies
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else :
                try : radii[pi] = vdwr['XXX'][aname[1]]
                except KeyError : radii[pi] = vdwr['XXX']['C']

    for index in res.keys() :
        for aname,pi in res[index].items() :
            if pi in knownPositions : continue
            pts[pi] = (0.,0.,0.)
    pts = VecVecFloat(pts)

    ai2pos, aiSC = [-999] * len(pts), [-999] * len(pts) ## flag atoms
    for index,val in res.items() :
        for name, ai in val.items() :
            ai2pos[ai] = index
            if name == ' N  ' : aiSC[ai] = 0
            elif name == ' CA ' : aiSC[ai] = 1
            elif name == ' C  ' : aiSC[ai] = 2
            elif name == ' O  ' : aiSC[ai] = 3
            elif name == ' SG ' and resns[index] == 'CYS' : aiSC[ai] = 5 # for disulphide vdw reduction
            elif name == ' CD ' and resns[index] == 'PRO' : aiSC[ai] = 6 # PRO CD shdnt be clash-checked with prev res's C
            else : aiSC[ai] = 4

    optRestraints = []
    if restrGen :
        xrlist, optional = restrGen.generate(blist, aiSC)
        irsize = len(irlist)
        for xr in xrlist : irlist.append(xr)
        for ori in optional : optRestraints.append(irsize + ori)

    scReduction, ssReduction = scReduction, data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    from peptidebuild import ModelRenderer
    modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, [], outpdb)


    assert len(radii) == len(pts)
    
    from PopulationStrategy import PopulationStrategy
    strategy = PopulationStrategy(backtrack, popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, 1, irlist, optRestraints, modelRenderer, nmodels*100, nmodels)
    strategy.execute()

def callmain() :
    from commonOptions import makeParser, parseOptions, addXrayOptions
    parser = makeParser()
    parser = addXrayOptions(parser)
    parser.add_option("--startres", action='store', type='string', dest='startResid', help='resid of start loop residue')
    parser.add_option("--endres", action='store', type='string', dest='endResid', help='resid of start loop residue')
    parser.remove_option("--buildN2C")
    options = parseOptions(parser)

    xrayRestGen = None
    if options.mtzfn :
        from prepareChainV5 import XrayRestraintsGenerator
        xrayRestGen = XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, options.sigXmin, options.sigXmax, options.sigXmean)

    main(options.pdbfile, options.startResid, options.endResid, options.mconly,
            options.caRad, options.scRad, options.scReduction, options.guidedSampling,
            options.popsize, options.backtrack, options.nmodels, options.outpdb,
            xrayRestGen)

def testMultiloop() :
    pdbfile = "1AAC.pdb"
    badresids = ["MET   28 "]
    badresids.append("LYS   29 ")
    badresids.append("THR   32 ")
    badresids.append("PRO   33 ")
    badresids.append("GLU   34 ")
    badresids.append("ASP    1 ")
    badresids.append("LYS    2 ")
    badresids.append("VAL  104 ")
    badresids.append("GLU  105 ")
    mconly = None ; caRad = 2 ; scRad = 2000 ; scReduction = 0.75 ; guidedSampling = 1 ; popsize = 100 ; backtrack = None ; nmodels = 1 ; outpdb = "out.pdb"
    import misc ; misc.setVerbosity(6)
    multiloop(pdbfile, None, mconly, caRad, scRad, scReduction, guidedSampling, popsize, backtrack, nmodels, outpdb, restrGen=None)


def testme() :
    badresids = []
    badresids.append("VAL   58 ")
    badresids.append("PRO   96 ")
    badresids.append("TYR   78 ")
    badresids.append("ALA   26 ")
    badresids.append("ALA   85 ")
    badresids.append("ASN   54 ")
    badresids.append("ASP   18 ")
    badresids.append("SER   79 ")
    badresids.append("PHE   97 ")
    badresids.append("LEU   35 ")
    badresids.append("GLU   64 ")
    badresids.append("MET   98 ")
    badresids.append("GLN   76 ")
    badresids.append("GLU   49 ")
    badresids.append("PRO   52 ")
    badresids.append("ASN   47 ")
    badresids.append("ILE    5 ")
    badresids.append("ALA    3 ")
    badresids.append("LYS   29 ")
    badresids.append("VAL   22 ")
    badresids.append("LYS    2 ")
    badresids.append("VAL   43 ")
    badresids.append("ILE   21 ")
    badresids.append("PRO   70 ")
    badresids.append("LYS   73 ")
    badresids.append("ARG   48 ")
    badresids.append("HIS   95 ")
    badresids.append("ILE   46 ")
    badresids.append("GLU   75 ")
    badresids.append("TRP   45 ")
    badresids.append("THR   87 ")
    badresids.append("VAL   39 ")
    badresids.append("PHE   82 ")
    badresids.append("GLU   15 ")
    badresids.append("VAL   37 ")
    badresids.append("VAL   16 ")
    badresids.append("ASP   41 ")
    badresids.append("VAL  104 ")
    badresids.append("LEU   67 ")
    badresids.append("ILE   25 ")
    badresids.append("GLY   40 ")
    badresids.append("ARG   99 ")
    badresids.append("PRO    6 ")
    multiloop("cns0.pdb", badresids, None, 3, 2000, 0.75, 1, 100, None, 1, "out.pdb", None)

def testme1() :
    badresids = None
    import prepareChainV5
    import misc ; misc.setVerbosity(5)
    xrayRestGen = [ prepareChainV5.XrayRestraintsGenerator("phased0.mtz", "FP", "FC", "PHIC", "2F1-F2", 0.1, 2, .5) ] ##loose
    multiloop("cnsC.pdb", badresids, None, 2, 2000, 0.5, 1, 100, "2X2", 1, "out.pdb", xrayRestGen)


if __name__ == "__main__" :
    testme1() ; sys.exit(0)
    callmain() ; sys.exit()
    testMultiloop()

    resids, chids, badresids = {}, {}, []
    resids[0] = "00" ; chids[0] = "A"
    resids[1] = "11" ; chids[1] = "A" ;
    resids[2] = "22" ; chids[2] = "A"
    resids[3] = "33" ; chids[3] = "A" ;
    resids[4] = "44" ; chids[4] = "B"
    resids[5] = "55" ; chids[5] = "B" ;
    resids[6] = "66" ; chids[6] = "B"
    resids[7] = "77" ; chids[7] = "B" ;
    resids[8] = "88" ; chids[8] = "B"

    badresids.append("00")
    #badresids.append("11")
    badresids.append("22")
    #badresids.append("33")
    badresids.append("44")
    #badresids.append("55")
    badresids.append("66")
    #badresids.append("77")
    badresids.append("88")

    loops, looptypes = locateRegionsRandomize(resids, chids, badresids)
    print loops, looptypes
