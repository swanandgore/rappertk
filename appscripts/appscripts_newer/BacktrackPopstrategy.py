from PopulationStrategy import PopulationStrategy
from builders import VecVecFloat, VecInt, VecFloat
from geometry import Grid
from random import random
import misc
import sys, time
from math import floor

class BacktrackPopstrategy (PopulationStrategy) :
    def execute_1_backup(s) :
        s.grids, s.ptslist = [], [] ## prepare grids and pts for all population members
        for pi in range(s.popsize) :
            s.ptslist.append( VecVecFloat(s.pts) )
            grid = Grid(3, s.ptslist[pi], s.radii, s.overlapReductions)
            grid.justAdd( VecInt(s.kp) )
            s.grids.append(grid)
        bi = 0
        maxBacktracks = 10 ; nfallbacks = {}
        while bi < len(s.builders) :
            ## remeber the valid parents before this try.
            ## we'r assuming that if build has built 0 children, there has been no copying
            if not s.build(bi) :
                if not bi in nfallbacks.keys() : nfallbacks[bi] = maxBacktracks
                if nfallbacks[bi] == 0 : return None
                newbi = s.builderFallIndex[bi]
                for bti in range(maxBacktracks, nfallbacks[bi], -1) :
                    if newbi < 0 : continue
                    newbi = s.builderFallIndex[newbi]
                if newbi < 0 : return None
                for bk in range(bi, newbi-1, -1) :
                    cbop = s.builders[bk].getOP()
                    for i in range(cbop.size()) :
                        for gi in range(s.popsize) : s.grids[gi].remove(cbop[i])
                nfallbacks[bi] = nfallbacks[bi] - 1
                print "BUILDER %d falling back onto %d" % (bi, newbi)
                bi = newbi
            else : bi = bi + 1
        ## choose a valid member of population and render
        pr = int(floor( random() * s.popsize ))
        s.mr.render(s.ptslist[pr])
        return 1

    def findKthFallback(s, bi, k) :
        tbi = bi
        for i in range(k) :
            tbi = s.builderFallIndex[tbi]
            if tbi < 0 : assert 0
        return tbi

    def findMaxFallbacks(s, bi) :
        maxPossFallbacks = 0
        tbi = bi
        while 1 :
            tbi = s.builderFallIndex[tbi]
            if tbi < 0 : return maxPossFallbacks
            maxPossFallbacks += 1
        
    def execute_1(s) :
        s.grids, s.ptslist = [], [] ## prepare grids and pts for all population members
        for pi in range(s.popsize) :
            s.ptslist.append( VecVecFloat(s.pts) )
            grid = Grid(3, s.ptslist[pi], s.radii, s.gridHelper)
            grid.justAdd( VecInt(s.kp) )
            s.grids.append(grid)
        print "START TIME", time.asctime()
        bi = 0
        maxBacktrackSteps = 2 ; stepsize = 5; nfallbacks = {}
        while bi < len(s.builders) :
            ## remeber the valid parents before this try.
            ## we'r assuming that if build has built 0 children, there has been no copying
            if s.build(bi) != 1 :
                print "BUILDER finding fallback"
                if not bi in nfallbacks.keys() : ## init fallback count
                    nfallbacks[bi] = []
                    ulim = maxBacktrackSteps * stepsize
                    maxPossFallbacks = s.findMaxFallbacks(bi)
                    if ulim > maxPossFallbacks : ulim = maxPossFallbacks
                    for k in range(ulim, 0, -1*stepsize) : nfallbacks[bi].append(k)
                    #nfallbacks[bi].append(0)
                if len(nfallbacks[bi]) == 0 :
                    #for pr in range(s.popsize) : s.mr.render(s.ptslist[pr]) #XXX
                    #sys.exit(1) #XXX
                    return None ## all fallbacks exhausted
                newbi = s.findKthFallback(bi, nfallbacks[bi].pop())
                for bk in range(bi, newbi-1, -1) :
                    cbop = s.builders[bk].getOP()
                    for i in range(cbop.size()) :
                        for gi in range(s.popsize) : s.grids[gi].remove(cbop[i])
                print "BUILDER %d falling back onto %d" % (bi, newbi)
                bi = newbi
            if bi in s.bfollowers.keys() :
                bi = bi + len(s.bfollowers[bi])
            bi = bi + 1
        ## choose a valid member of population and render
        #for pr in range(s.popsize) : s.mr.render(s.ptslist[pr]) #XXX
        #sys.exit(0) #XXX
        pr = int(floor( random() * s.popsize ))
        s.mr.render(s.ptslist[pr])
        return 1
