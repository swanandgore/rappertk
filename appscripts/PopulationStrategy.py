from BasicStrategy import BasicStrategy
from builders import VecVecFloat, VecInt, VecFloat
from geometry import Grid, MapIntIntFloat
from random import random, seed
seed(1973) # for Paul!
from math import floor, sqrt
import time, math, string, copy

from misc import verbose

class PopulationStrategy (BasicStrategy) :
    '''this is rapper-like population-based search strategy.
        a pool of built models is maintained at every step of build process.
        if a model is not extendable, its slot is taken up by an extendable (read fitter) model.
        unlike BasicStrategy, multiple grids and pointsets need to be maintained.
        python documenting within triple-quotes is delightfully simple and useful!
    '''
    def __init__(s, backtrack, popSize, *args) :
        BasicStrategy.__init__(s, *args)
        s.popsize = popSize ; assert s.popsize > 0
        s.numBacktrackSteps, s.backtrackStepsize = None, None
        if backtrack != None : ## this is supposed to be aXb, numstepsXstepsize, eg 4X5
            xpos = backtrack.find('X')
            assert xpos > 0
            s.numBacktrackSteps = string.atoi(backtrack[0:xpos])
            s.backtrackStepsize = string.atoi(backtrack[xpos+1:])
        s.useEqOpp = None

    def execute(s) :
        for i in range(s.nattempts) : # XXX
            startTime = time.time()
            s.snap = None
            if s.execute_1() == 1 :
                s.nmodels = s.nmodels-1
                print "ATTEMPT %d successfully generated population" % i, time.time()-startTime
            else : print "ATTEMPT %d of generating a population failed" % i, time.time()-startTime
            if s.nmodels == 0 : break

    def initRestraintFailures(s, bi) :
        bis = [bi]
        if bi in s.bfollowers.keys() : bis = bis + s.bfollowers[bi]
        for i in bis :
            for ri in s.bldr2restraints[i] :
                if not ri in s.restraintFailures.keys() : s.restraintFailures[ri] = 0
                s.restraintFailures[ri] = 0

    def printRestraintFailures(s, bi) :
        #return
        bis, totfail = [bi], 0
        if bi in s.bfollowers.keys() : bis = bis + s.bfollowers[bi]
        for i in bis :
            for ri in s.bldr2restraints[i] :
                print "%8d failures for " % s.restraintFailures[ri], s.restraints[ri].name()
                totfail += s.restraintFailures[ri]
        print "%8d total failures" % totfail

    def discardOptionalRestraints(s, bi) :
        bis = [bi]
        if bi in s.bfollowers.keys() : bis = bis + s.bfollowers[bi]
        discarded = None
        for i in bis :
            for ri in s.bldr2restraints[i] :
                if ri in s.optionalRestr and s.restraintFailures[ri] > 0 :
                    s.discardedRestraints.add(ri)
                    print "discarding optional restraint", s.restraints[ri].name(), "with %d failures" % s.restraintFailures[ri]
                    discarded = 1
        return discarded

    def writeSnapshot(s, bi) :
        if s.snap == None : return
        s.snap += 1
        bptis = set()
        for b in s.builders[0:bi+1] :
            bop = b.getOP()
            for i in range(bop.size()) : bptis.add(bop[i])
        newmr = copy.deepcopy(s.mr)
        newmr.f.close()
        newmr.f = open("snap%d.pdb" % s.snap, 'w')
        for pr in s.validParents :
            #s.mr.render(s.ptslist[pr], bptis)
            newmr.render(s.ptslist[pr], bptis)

    def execute_1(s) :
        s.restraintFailures = {}
        s.discardedRestraints = set()
        s.validParents = range(s.popsize)
        s.grids, s.ptslist = [], [] ## prepare grids and pts for all population members
        for pi in range(s.popsize) :
            s.ptslist.append( VecVecFloat(s.pts) )
            grid = Grid(3, s.ptslist[pi], s.radii, s.gridHelper)
            grid.justAdd( VecInt(s.kp) )
            s.grids.append(grid)
        print "START TIME", time.asctime()
        bi = 0 ; numFallbacks = {}
        while bi < len(s.builders) :
            succeeded = 1
            if s.build(bi) != 1 :
                succeeded = None
                if s.discardOptionalRestraints(bi) :
                    if s.build(bi) == 1 : succeeded = 1
                    #bptis = set()
                    #for b in s.builders[0:bi+1] :
                        #bop = b.getOP()
                        #for i in range(bop.size()) : bptis.add(bop[i])
                    #for pr in s.validParents : s.mr.render(s.ptslist[pr], bptis) #XXX
                    #import sys ; sys.exit(1) #XXX
                    #return None
            if succeeded :
                if bi in s.bfollowers.keys() :
                    bi = bi + len(s.bfollowers[bi])
                bi = bi + 1
            elif s.backtrackStepsize != None :
                if not numFallbacks.has_key(bi) : numFallbacks[bi] = 0
                if numFallbacks[bi] >= s.numBacktrackSteps : return None
                numFallbacks[bi] += 1 ; oldbi = bi
                for i in range(numFallbacks[bi]*s.backtrackStepsize) :
                    if bi == 0 : break
                    bi = s.builderFallIndex[bi]
                assert bi >= 0
                print "BUILDER", oldbi, "falling onto", bi
            else : return None
        ## choose any member of valid parents and render
        pr = s.validParents[ int(floor( random() * len(s.validParents) )) ]
        s.mr.render(s.ptslist[pr])
        return 1

    ## return 1 on successul build, -1 on SphPos-restraint failure, -2 on clash failure, -3 if builder refuses, -4 on RATrestraint failure, -5 on other restraint failures
    def tryBuilderOnMember(s, bldr, bi, pi) :
        if verbose(7) : print "TRYING builder on member %d" % (pi)
        #print "TRYING builder", bldr.name(), "on member %d" % (pi), bldr.sessionSize
        if not bldr.build(s.ptslist[pi]) : ## call builder
            if verbose(7) : print "builder refused to build"
            return -3
        for ri in s.bldr2restraints[bi] :
            restr = s.restraints[ri]
            if ri in s.discardedRestraints : continue
            if not restr.satisfied(s.ptslist[pi]) :
                s.restraintFailures[ri] += 1
                if verbose(7) : print "restraint failed"
                rtype = "%s" % type(restr)
                if 'SphPos' in rtype : return -1
                elif 'RAT' in rtype : return -4
                else : return -5
        opinds = bldr.getOP()
        if s.grids[pi].add(opinds) == 0 : # clash-check XXX
            if verbose(7) : print "clashcheck failed"
            return -2
        for opi in range(opinds.size()) : s.grids[pi].remove(opinds[opi]) # remove points
        if verbose(7) : print "TRIED builder %d on member %d :-)" % (bi,pi)
        return 1

    def recursiveBuild(s, bclones, i, bis, pi) :
        #print "REC", i, bis
        ret = s.tryBuilderOnMember(bclones[bis[i]][pi], bis[i], pi)
        if ret != 1 : return ret
        if i+1 == len(bis) : return ret
        bclones[bis[i+1]][pi].endSession()
        bclones[bis[i+1]][pi].startSession(1000000)
        numUniqueTries, count, maxcount = 0, 0, 10*s.numtrials[bis[i+1]]
        while numUniqueTries < s.numtrials[bis[i+1]] and count < maxcount :
            ret = s.recursiveBuild(bclones, i+1, bis, pi)
            if ret == 1 : return ret
            elif ret != -3 : numUniqueTries += 1
            count += 1
        return -5

    ## try extending the population by building bi'th builder
    def build(s,bi) :
        buildRefusal, gridClashes = 0, 0
        s.initRestraintFailures(bi)
        opinds = set()
        opis = s.builders[bi].getOP()
        for i in range(opis.size()) : opinds.add(opis[i])
        bclones, bis = {bi:[]}, [bi]
        if bi in s.bfollowers.keys() :
            for bk in s.bfollowers[bi] :
        	opis = s.builders[bk].getOP()
        	for i in range(opis.size()) : opinds.add(opis[i])
                bclones[bk] = [] ; bis.append(bk)
        for pi in range(s.popsize) :
            for bk in bclones.keys() :
                popbldr = s.builders[bk].makeCopy()
                popbldr.endSession()
                popbldr.startSession(10000000) #s.numtrials[bi])
                bclones[bk].append(popbldr)
        opinds = VecInt(list(opinds))
        numTriesExceptBuildfails = s.numtrials[bi] * len(s.validParents) - 1
        eqOppTries = numTriesExceptBuildfails * 95/100
        maxChildren = s.popsize
        children, totch = {}, 0
        #print numTriesExceptBuildfails, eqOppTries, maxChildren
        maxpasses = int( s.numtrials[bi] * ((0.+s.popsize) / len(s.validParents)) )
        for passnum in range(maxpasses) :
            #print "passnum", passnum, "of", maxpasses, s.numtrials[bi], s.popsize, len(s.validParents)
            if totch > maxChildren : break
            for pi in s.validParents :
                if totch > maxChildren : break
                if s.useEqOpp and pi in children.keys() and numTriesExceptBuildfails > eqOppTries : continue
                buildStatus = -3
                if numTriesExceptBuildfails >= 0 :
                    while buildStatus == -3  : buildStatus = s.recursiveBuild(bclones, 0, bis, pi)
                else :
                    buildStatus = s.recursiveBuild(bclones, 0, bis, pi)
                if buildStatus == -3 : buildRefusal += 1
                else : numTriesExceptBuildfails -= 1 # num of non-build-refusal attempts
                if buildStatus == -2 : gridClashes += 1
                if buildStatus != 1 : continue
                child = VecVecFloat()
                for opi in range(opinds.size()) : child.push_back( s.ptslist[pi][opinds[opi]] )
                if not pi in children.keys() : children[pi] = [] # make space for this parent's children
                children[pi].append(child) ; totch = totch + 1
        for pi in range(s.popsize) :
            for bk in bis : bclones[bk][pi].endSession()
        if verbose(6) : s.printRestraintFailures(bi)
        print "BUILDER %4d : %4d parents : %4d children : (r%6d c%6d) : " % (bi, len(children.keys()), totch, buildRefusal, gridClashes),
        print s.builders[bi].name(), "-------------------------"
        fertile = list(children.keys())
        if len(fertile) == 0 : return None
        barren = []
        for pi in range(s.popsize) :
            if not pi in fertile : barren.append(pi)
        s.validParents = []
        while 1 : ## fill barren slots
            if len(barren) == 0 or totch == len(children) : break
            for fi in fertile :
                if len(barren) == 0 or totch == len(children) : break
                if len(children[fi]) == 1 : continue
                bri = barren.pop() ; s.validParents.append(bri)
                assert fi != bri
                child = children[fi].pop()
                s.grids[bri].popCopyGridPoints( s.grids[fi] )
                totch = totch - 1
                for opi in range(opinds.size()) :
                    s.ptslist[bri][opinds[opi]] = child[opi] #VecFloat(child[opi])
                assert s.grids[bri].add(opinds) != 0
        for fi in fertile : ## fill fertile slots
            assert len(children[fi]) >= 1
            child = children[fi][0]
            for opi in range(opinds.size()) :
                s.ptslist[fi][opinds[opi]] = child[opi]
            assert s.grids[fi].add(opinds) != 0
            s.validParents.append(fi)
        return 1
