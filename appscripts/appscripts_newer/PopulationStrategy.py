from BasicStrategy import BasicStrategy
from builders import VecVecFloat, VecInt, VecFloat
from geometry import Grid, MapIntIntFloat
from random import random
from math import floor, sqrt
import time, math, string, copy

from misc import verbose


def samecrds(crds1, crds2, crdinds=None) :
    if crdinds == None :
        try : numcrds = crds1.size()
        except : numcrds = len(crds1)
        crdinds = range(numcrds)
    for i in crdinds :
        if abs( crds1[i][0] - crds2[i][0] ) > 1e-3 : return None
        if abs( crds1[i][1] - crds2[i][1] ) > 1e-3 : return None
        if abs( crds1[i][2] - crds2[i][2] ) > 1e-3 : return None
    return 1

def nextOrder(ord, maxs) :
    for i in range(len(ord)) :
        if ord[i] + 1 == maxs[i] :
            ord[i] = 0
            if i+1 == len(ord) : return None
        else : ord[i] += 1 ; break
    return ord
def makeOrders(maxsizes) :
    for ms in maxsizes :
        if ms <= 0 : return []
    allzero = [0,] * len(maxsizes)
    orders = [ allzero ] ## set of orders to return
    while(1) :
        ord = nextOrder( list(orders[ len(orders)-1 ]) , maxsizes )
        if ord == None : break
        orders.append(ord)
    #for ord in orders : print ord
    return orders
def midOrders(nextlevel, maxsizes) :
    ms0 = list(maxsizes)
    ms1 = list(maxsizes)
    for i in range(len(ms0)) :
        if ms0[i] > nextlevel-1 : ms0[i] = nextlevel-1
        if ms1[i] > nextlevel : ms1[i] = nextlevel
    ord0 = makeOrders(ms0)
    ord1 = makeOrders(ms1)
    #print ms0, len(ord0), ms1, len(ord1)
    for ord in ord0 : ord1.remove(ord)
    return ord1

class PopulationStrategy (BasicStrategy) :
    '''this is rapper-like population-based search strategy.
        a pool of built models is maintained at every step of build process.
        if a model is not extendable, its slot is taken up by an extendable (read fitter) model.
        unlike BasicStrategy, multiple grids and pointsets need to be maintained.
        python documenting within triple-quotes is delightfully simple and useful!
    '''
    def __init__(s, backtrack, popSize, *args) :
        BasicStrategy.__init__(s, *args)
        s.ranker = None
        s.popsize = popSize ; assert s.popsize > 0
        s.numBacktrackSteps, s.backtrackStepsize = None, None
        if backtrack != None : ## this is supposed to be aXb, numstepsXstepsize, eg 4X5
            xpos = backtrack.find('X')
            assert xpos > 0
            s.numBacktrackSteps = string.atoi(backtrack[0:xpos])
            s.backtrackStepsize = string.atoi(backtrack[xpos+1:])
        s.useEqOpp = None
        s.snap = 0 ; s.snapRenderer = None

    def execute(s) :
        numSucc = 0
        for i in range(s.nattempts) : # XXX
            startTime = time.time()
            s.snap = 0
            if s.execute_1() == 1 :
                s.nmodels = s.nmodels-1 ; numSucc += 1
                print "ATTEMPT %d successfully generated population" % i, time.time()-startTime
            else : print "ATTEMPT %d of generating a population failed" % i, time.time()-startTime
            if s.nmodels == 0 : break
        return numSucc

    def initRestraintFailures(s, bi) :
        bis = [bi]
        if bi in s.bfollowers.keys() : bis = bis + s.bfollowers[bi]
        for i in bis :
            for ri in s.bldr2restraints[i] :
                if not ri in s.restraintFailures.keys() : s.restraintFailures[ri] = [0,0] ## failed, checked
                s.restraintFailures[ri] = [0,0]

    def printRestraintFailures(s, bi) :
        #return
        bis, totfail = [bi], 0
        if bi in s.bfollowers.keys() : bis = bis + s.bfollowers[bi]
        for i in bis :
            for ri in s.bldr2restraints[i] :
                #if not ri in s.optionalRestr : continue
                print "%8d/%8d failures for " % (s.restraintFailures[ri][0],s.restraintFailures[ri][1]), s.restraints[ri].name()
                totfail += s.restraintFailures[ri][0]
        print "%8d total failures" % totfail

    def discardOptionalRestraints(s, bi) :
        bis = [bi]
        if bi in s.bfollowers.keys() : bis = bis + s.bfollowers[bi]
        discarded = None
        riMaxviol = -1
        for i in bis :
            for ri in s.bldr2restraints[i] :
                if ri in s.optionalRestr and s.restraintFailures[ri][0] > 0 and ri not in s.discardedRestraints :
                    if riMaxviol == -1 : riMaxviol = ri
                    elif s.restraintFailures[ri][0]/(0.+s.restraintFailures[ri][1]) > s.restraintFailures[riMaxviol][0]/(0.+s.restraintFailures[riMaxviol][1]) : riMaxviol = ri
        if riMaxviol > -1 : ## discard the one with highest violations
            s.discardedRestraints.add(riMaxviol)
            print "discarding optional restraint", s.restraints[riMaxviol].name(), "with %d/%d failures" % (s.restraintFailures[riMaxviol][0],s.restraintFailures[riMaxviol][1])
            discarded = 1
        return discarded

    def writeSnapshot(s, bi) :
        if s.snapRenderer == None : return
        s.snap += 1
        bptis = set()
        curbldrs = s.builders[0:bi+1]
        if bi in s.bfollowers.keys() :
            for bk in s.bfollowers[bi] : curbldrs.append( s.builders[bk] )
        for b in curbldrs :
            bop = b.getOP()
            for i in range(bop.size()) : bptis.add(bop[i])
        newmr = copy.deepcopy(s.snapRenderer)
        newmr.f = open("snap%d.pdb" % s.snap, 'w')
        print "SNAP", s.snap
        for pr in s.validParents :
            newmr.render(s.ptslist[pr], bptis)

    def execute_1(s) :
        s.restraintFailures = {}
        s.discardedRestraints = set()
        s.validParents = [0] #range(s.popsize)
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
                while succeeded==None and s.discardOptionalRestraints(bi) :
                    if s.build(bi) == 1 : succeeded = 1
                    #bptis = set()
                    #for b in s.builders[0:bi+1] :
                        #bop = b.getOP()
                        #for i in range(bop.size()) : bptis.add(bop[i])
                    #for pr in s.validParents : s.mr.render(s.ptslist[pr], bptis) #XXX
                    #import sys ; sys.exit(1) #XXX
                    #return None
            if succeeded :
                s.writeSnapshot(bi)
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
                for bk in range(oldbi, bi-1, -1) : ## remove op points belonging to builders between oldbi and bi from grids
                    cbop = s.builders[bk].getOP()
                    for i in range(cbop.size()) :
                        for gi in range(s.popsize) : s.grids[gi].remove(cbop[i])
            else : return None
        ## choose any member of valid parents and render, or if ranker is given, choose the lowest scoring
        if s.ranker == None :
            pr = s.validParents[ int(floor( random() * len(s.validParents) )) ]
        else :
            ptinds, crds = [], []
            for b in s.builders :
                bop = b.getOP()
                for i in range(bop.size()) : ptinds.append(bop[i])
            ## init the score with given structure, ie if given is better, ensemble will be ignored
            pr = -1 ; minscore = 1e10
            if s.ranker.rankGivenEnsemble :
                for pti in ptinds : crds.append( s.pts[pti] )
                minscore = s.ranker.score(crds)
            for pi in s.validParents :
                crds = []
                for pti in ptinds : crds.append( s.ptslist[pi][pti] )
                sc = s.ranker.score(crds)
                if minscore > sc : minscore = sc ; pr = pi
                #print pr, minscore, sc
        if pr >= 0 :
            s.mr.render(s.ptslist[pr])
            if s.ranker != None and s.ranker.rankGivenEnsemble : print "Ensemble contains a better member than given structure"
            elif s.ranker != None : print "Choosing best member from ensemble"
        else : print "Given structure was better than all ensemble members."
        return 1

    ## return 1 on successul build, -1 on SphPos-restraint failure, -2 on clash failure, -3 if builder refuses, -4 on RATrestraint failure, -5 on other restraint failures
    def tryBuilderOnMember(s, bldr, bi, pi, rotind=None) :
        if verbose(7) : print "TRYING builder on member %d" % (pi)
        #print "TRYING builder", bldr.name(), "on member %d" % (pi), bldr.sessionSize
        if rotind : bstat = bldr.buildSample(s.ptslist[pi],rotind) ## call builder
        else : bstat = bldr.build(s.ptslist[pi])
        if not bstat :
            if verbose(7) : print "builder refused to build"
            return -3
        for ri in s.bldr2restraints[bi] :
            restr = s.restraints[ri]
            if ri in s.discardedRestraints : continue
            s.restraintFailures[ri][1] += 1
            if not restr.satisfied(s.ptslist[pi]) :
                s.restraintFailures[ri][0] += 1
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

    def sc3(s, bclones, bis, pi) :
        scores, rinds = {}, {}
        for bi in bis[1:] : ## rank the rotamers of given builders in ascending order of ranker score
            bclones[bi][pi].endSession() ; bclones[bi][pi].startSession(1000000)
            energies, ris = [], [] # lower energy is better energy
            for ri in range(100000) :
                if bclones[bi][pi].buildSample(s.ptslist[pi], ri) == 0 : break
                bclones[bi][pi].endSession() ; bclones[bi][pi].startSession(1000000)
                if s.tryBuilderOnMember(bclones[bi][pi], bi, pi, ri) != 1 : continue
                child = [] ; opinds = bclones[bi][pi].getOP()
                for ci in range(opinds.size()) : child.append(s.ptslist[pi][ opinds[ci] ])
                energies.append( s.ranker.score(child) ) ; ris.append(ri)
            for i in range(len(energies)) :
                for k in range(i+1, len(energies)) :
                    if energies[i] > energies[k] :
                        temp = energies[i] ; energies[i] = energies[k] ; energies[k] = temp
                        temp = ris[i] ; ris[i] = ris[k] ; ris[k] = temp
            scores[bi], rinds[bi] = energies, ris
        ## go on trying various combinations of rotamers in ascending energy order and break when successful
        order, maxsizes = [], []
        for bi in bis[1:] : order.append(0) ; maxsizes.append( len(scores[bi]) )
        ret = -3
        for oi in range(1,1000) :
            orders, ens = midOrders(oi, maxsizes), [] ## find orders and order on combined energy
            for ord in orders :
                en = 0.
                for i in range(len(bis)-1) :
                    bi = bis[i+1]
                    en += scores[bi][ord[i]]
                ens.append(en)
            for i in range(len(orders)) :
                for k in range(i+1,len(orders)) :
                    if ens[i] > ens[k] :
                        ten = ens[i] ; ens[i] = ens[k] ; ens[k] = ten
                        tord = orders[i] ; orders[i] = orders[k] ; orders[k] = tord
            for order in orders :
                ret = -3
                for i in range(len(bis)-1) :
                    bi = bis[i+1]
                    bclones[bi][pi].endSession() ; bclones[bi][pi].startSession(1000000)
                    ret = s.tryBuilderOnMember(bclones[bi][pi], bi, pi, rinds[bi][order[i]])
                    if ret != 1 : break
                    assert 0 != s.grids[pi].add(bclones[bi][pi].getOP())
                for bi in bis[1:] : s.grids[pi].removes(bclones[bi][pi].getOP())
                if ret == 1 : break ; print oi, len(orders), order ; break
            if ret == 1 : break
        return ret

    # returns a map of pi->[pj1,pj2,...] s.t. pj's have same input coordinates for bldr bi as pi
    def findDistinctValidParents(s, bi) :
        notdistinct = [] ## check if input points are identical for any 2 parents and if so, prevent one of them from making children
        biip, pis = s.builders[bi].getIP(), []
        if biip.size() == 0 : return list(s.validParents), {}
        for i in range(biip.size()) : pis.append(biip[i])
        same = {}
        for pi in s.validParents : same[pi] = []
        for pi in s.validParents :
            for pj in s.validParents :
                if pi == pj or pj in notdistinct or pi in notdistinct : continue
                aresame = 1
                aresame = samecrds(s.ptslist[pi], s.ptslist[pj], pis)
                if aresame==1 : same[pi].append(pj) ; notdistinct.append(pj)
        for pi,v in same.items() :
            if len(v) == 0 : del same[pi]
        allsame = {}
        for pi,v in same.items() :
            samegroup = list( [pi,] + list(v) )
            for i in range(len(samegroup)) :
                key, vals = samegroup[i], list(samegroup)
                vals.remove( samegroup[i] )
                allsame[key] = vals
                #print key, vals
        vps = list(s.validParents)
        for pi in notdistinct : vps.remove(pi)
        return vps, allsame

    def recursiveBuild(s, bclones, i, bis, pi) :
        #print "REC", i, bis
        ret = s.tryBuilderOnMember(bclones[bis[i]][pi], bis[i], pi)
        if ret != 1 : return ret
        if i+1 == len(bis) : return ret
        dorecurse = None
        if s.ranker and s.ranker.rankChildren and s.ranker.rankRecurse : dorecurse = 1
        assert 0 != s.grids[pi].add( bclones[bis[i]][pi].getOP() )
        bclones[bis[i+1]][pi].endSession()
        bclones[bis[i+1]][pi].startSession(1000000)
        if "PeptideBridgeBuilder" in s.builders[bis[0]].name() and i == 0 and dorecurse : # special treatment for 3-sc case of PeptideBridgeBuilder
            ret = s.sc3(bclones, bis, pi)
            s.grids[pi].removes( bclones[bis[i]][pi].getOP() ) ; return ret
        numUniqueTries, count, maxcount = 0, 0, 10*s.numtrials[bis[i+1]]
        if dorecurse :
            opinds = [] ; bestchild = None
            for bi in bis[i+1:] :
                opis = s.builders[bi].getOP()
                for ci in range(opis.size()) : opinds.append(opis[ci])
        while numUniqueTries < s.numtrials[bis[i+1]] and count < maxcount :
            ret = s.recursiveBuild(bclones, i+1, bis, pi)
            if ret != -3 : numUniqueTries += 1
            if ret == 1 and dorecurse :
                child = []
                for crdi in opinds : child.append(s.ptslist[pi][crdi])
                if bestchild == None or s.ranker.score(child) < s.ranker.score(bestchild) : bestchild = child
            elif ret == 1 : s.grids[pi].removes( bclones[bis[i]][pi].getOP() ) ; return ret
            count += 1
        if dorecurse and bestchild != None :
            for ci in range(len(opinds)) : s.ptslist[pi] [ opinds[ci] ] = bestchild[ci]
            s.grids[pi].removes( bclones[bis[i]][pi].getOP() ) ; return 1
        s.grids[pi].removes( bclones[bis[i]][pi].getOP() ) ; return -5

    ## try extending the population by building bi'th builder
    def build(s,bi) :
        buildRefusal, gridClashes = 0, 0
        s.initRestraintFailures(bi)
        opinds = list()
        opis = s.builders[bi].getOP()
        for i in range(opis.size()) : opinds.append(opis[i])
        bclones, bis = {bi:[]}, [bi]
        if bi in s.bfollowers.keys() :
            for bk in s.bfollowers[bi] :
                opis = s.builders[bk].getOP()
                for i in range(opis.size()) : opinds.append(opis[i])
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
        if s.ranker != None and s.ranker.rankChildren != None :
            maxChildren *= s.ranker.rankChildren ; s.ranker.curbest, s.ranker.cutoff = 1e10, 1e10
            if "PeptideBridgeBuilder" in s.builders[bi].name() : maxChildren *= 100
        children, totch = {}, 0
        #print numTriesExceptBuildfails, eqOppTries, maxChildren
        vps, sameparents = s.findDistinctValidParents(bi)
        #print "validparents", len(vps), len(s.validParents)
        maxpasses = int( s.numtrials[bi] * ((0.+s.popsize) / len(s.validParents)) )
        for passnum in range(maxpasses) :
            #import sys ; print ".", ; sys.stdout.flush()
            #print "passnum", passnum, "of", maxpasses, s.builders[bi].name(), s.numtrials[bi], s.popsize, len(s.validParents), totch, maxChildren
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
                ## parents with identical input to this builder are likey to provide same children. is this a diff child really?
                newchild = 1
                if len(vps) == 1 and len(s.validParents) > 1 :
                    for pj in sameparents[pi] :
                        if not newchild : break
                        if not children.has_key(pj) : continue
                        for ch in children[pj] :
                            if samecrds(ch, child) : newchild = None ; break
                if not newchild : continue ; print "rejecting child", pi, pj ; continue
                if not pi in children.keys() : children[pi] = [] # make space for this parent's children
                children[pi].append(child) ; totch = totch + 1
        nch = 0
        for pi in children.keys() : nch += len(children[pi])
        if s.ranker != None and s.ranker.rankChildren != None and totch > s.popsize : ## remove low-scoring children
            chscores = {} ; allscores = []
            for pi in children.keys() : ## lower score is better score
                if len(children[pi]) == 0 : continue
                chscores[pi] = []
                for child in children[pi] :
                    if s.ranker.rankLeaderBuilderOnly : allscores.append ( s.ranker.score(child, s.builders[bi].getOP().size()) )
                    else : allscores.append ( s.ranker.score(child) )
                    chscores[pi].append( allscores[len(allscores)-1] )
            allscores.sort() ; chCut = allscores[s.popsize]
            #print allscores, chCut
            totch = 0
            for pi in children.keys() :
                newch = []
                for i in range(len(chscores[pi])) :
                    if chscores[pi][i] < chCut : newch.append( children[pi][i] ) ; totch += 1
                children[pi] = newch
            for pi in children.keys() :
                if len(children[pi]) == 0 : del children[pi]
        for pi in range(s.popsize) :
            for bk in bis : bclones[bk][pi].endSession()
        if verbose(6) : s.printRestraintFailures(bi)
        print "BUILDER %4d : %4d parents : %4d/%6d children : (r%6d c%6d) : " % (bi, len(children.keys()), totch, nch, buildRefusal, gridClashes),
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