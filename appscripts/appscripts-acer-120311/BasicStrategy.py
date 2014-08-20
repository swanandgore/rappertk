from geometry import Grid, VecVecInt, VecInt
from builders import VecVecFloat, VecFloat
from misc import verbose

class Strategy :
    '''the strategy class which all other strategies shd extend from'''
    def __init__(s, pts, radii, knownPositions, gridHelper, blist, numtrials, bfollowers, rlist, optR, modelRenderer, natt, nmod) :
        s.builders, s.numtrials, s.bfollowers, s.restraints, s.optionalRestr = blist, numtrials, bfollowers, rlist, optR
        s.mr = modelRenderer
        s.pts, s.radii, s.kp = VecVecFloat(pts), VecFloat(radii), knownPositions
        s.nattempts, s.nmodels = natt, nmod
        s.gridHelper = gridHelper

    def execute(s) : pass

class BasicStrategy (Strategy) :
    '''this is the most basic strategy of exhaustive search through the tree of conformational options. useful for small tasks and testing.'''
    def __init__(s, pts, radii, knownPositions, gridHelper, blist, numtrials, bfollowers, reorderBuilders, rlist, optRs, modelRenderer, natt, nmod, knownORbuilt=1) :
        Strategy.__init__(s, pts, radii, knownPositions, gridHelper, blist, numtrials, bfollowers, rlist, optRs, modelRenderer, natt, nmod)
        s.knownORbuilt = knownORbuilt # verify that each atom is known or built
        if reorderBuilders :
            s.reorderBuilders()
        else :
            import cbutils
            cbutils.builderSanityCheck(blist, s.pts.size(), knownPositions, s.knownORbuilt)
            cbutils.checkBuilderOrder(blist, knownPositions)
            s.builderFallIndex = [-1] * len(s.builders)
            bfolls, blead = [], []
            for k,v in s.bfollowers.items() : bfolls += list(v)
            for bi in range(len(s.builders)) :
                if not bi in bfolls : blead.append(bi)
            blead.sort()
            for bli in range(len(blead)-1, 0, -1) :
                s.builderFallIndex[ blead[bli] ] = blead[bli-1]
            for bli in range(len(blead)-1, -1, -1) :
                bl = blead[bli]
                if not bl in s.bfollowers.keys() : continue
                for fi in s.bfollowers[bl] : s.builderFallIndex[fi] = bl
            s.builderFallIndex[0] = -1
            #print s.builderFallIndex ;# assert None
        s.reorderRestraints()

    def reorderBuilders(s) :
        '''topo-sort and reorder the builders'''
        from cbutils import findBuilderOrder
        bwt = [-1] * len(s.builders)
        for bi in s.bfollowers.keys() :
            for fi in s.bfollowers[bi] : bwt[fi] = 0
        builderOrder, builderFallback, builderFallIndex = findBuilderOrder(s.builders, bwt, s.pts.size(), s.kp, s.knownORbuilt, writeDigraph='dotfile')
#        if verbose(6) : print '---------------final order of builders----------------------------------------------'
        for i in range(len(builderOrder)) :
            bi = builderOrder[i]
#            if verbose(6) : print bi, " ",
            currBuilder = s.builders[bi]
#            if verbose(6) : currBuilder.describe()
            cbip, cbop = currBuilder.getIP(), currBuilder.getOP()
#            if verbose(6) :
#                print '--- fallback %d %d --- ' % (builderFallback[i], builderFallIndex[i]),
#                for i in range(cbip.size()) : print cbip[i],
#                print ' ---> ',
#                for i in range(cbop.size()) : print cbop[i],
#                print ''
#        if verbose(6) : print '------------------------------------------------------------------------------------\n\n\n'
        newbuilders, newnt, fallback, bfoll = [], [], [], {}
        for i in range(len(builderOrder)) :
            bi = builderOrder[i]
            newbuilders.append( s.builders[bi] )
            newnt.append( s.numtrials[bi] )
            if not bi in s.bfollowers.keys() : continue # ds this builder have followers?
            bfoll[i] = []
            # is follower order maintained? are s.bfollowers[bi] same as ones immediately following i in builderOrder ?
            givenFoll = list(s.bfollowers[bi]) ; givenFoll.sort()
            neword = list(builderOrder[i+1:i+len(s.bfollowers[bi])+1]) ; neword.sort()
            if givenFoll != neword :
                print "Follower order is not obeyed by builder reordering : ", s.builders[bi].name(), givenFoll, neword
                assert None
            for bfi in range(len(s.bfollowers[bi])) : bfoll[i].append(i+bfi+1)
            #print bfoll[i], "follow", i
        s.builders, s.numtrials, s.bfollowers = newbuilders, newnt, bfoll
        s.builderFallIndex = builderFallIndex
#        if verbose(6) :
#            print "\n\n\n"
#            for bi in range(len(s.builders)) :
#                if bi in s.bfollowers.keys() : print s.builders[bi].name(), "leads", s.builders[s.bfollowers[bi][0]].name()
#            print "\n\n\n"

    def reorderRestraints(s) :
        '''find which restraints to check after which builder'''
        newrestraints = []
        for i in range(len(s.builders)) : newrestraints.append( set([]) )
        # find source builder for each point
        ptinds = set() ## find number of points
        for b in s.builders :
            for i in range(b.getIP().size()) : ptinds.add( b.getIP()[i] )
            for i in range(b.getOP().size()) : ptinds.add( b.getOP()[i] )
        srcBldr = {}
        for k in ptinds : srcBldr[k] = -999
        for bi in range(len(s.builders)) :
            b = s.builders[bi]
            bop = b.getOP()
            for i in range(bop.size()) :
                assert srcBldr[ bop[i] ] == -999
                srcBldr[ bop[i] ] = bi
        
        for ri in range(len(s.restraints)) :
            r = s.restraints[ri]
            #print r, srcBldr, r.name()
            rinds, highestBI = r.getInds(), -999
            #print "size", rinds.size()
            for i in range(rinds.size()) :
                if s.knownORbuilt == None and not srcBldr.has_key(rinds[i]) : continue
                bi = srcBldr[ rinds[i] ]
                #print rinds[i], rinds.size(), i, bi
                if highestBI < bi : highestBI = bi
            if s.knownORbuilt == None and highestBI == -999 : continue
            assert highestBI != -999 and highestBI >= 0
            newrestraints[highestBI].add(ri)
        s.bldr2restraints = newrestraints
        #print 'RLIST----------------------'
        #for rs in s.restraints :
        #    print len(rs),
        #    for r in rs : print r,
        #    print ''
        #print '---------------------------\n\n\n'

    def execute(s) :
        ## start model building
        for nm in range(s.nattempts) :
            print "starting MODEL ", nm
            ## add starting points to grid unconditionally
            s.grid = Grid(3, s.pts, s.radii, s.overlapReductions)
            s.grid.justAdd( VecInt(s.kp) )
            #grid.justAdd(blist[ builderOrder[0] ].getIP())
            ## start building
            s.numtrialsUsed = [0] * len(s.numtrials)
            bout = s.build()
            print "BUILD STATUS", bout
            if bout == 1 :
                s.mr.render()
                s.nmodels -= 1
            if s.nmodels == 0 : break

    def findFallback(s, bi) :
        fi = s.builderFallIndex[bi]
        candidate = -1
        if fi == -1 : # have no-one to fall upon, find earlier bi with fi = -1
            for bk in range(bi-1, -1, -1) :
                if s.builderFallIndex[bk] == -1 :
                    candidate = bk
                    break
        else : candidate = fi
        
        if candidate == -1 : return candidate
        if s.numtrialsUsed[candidate] >= s.numtrials[candidate] : 
            return s.findFallback(candidate)
        return candidate

    def build(s) :
        bi = 0
        while bi < len(s.builders) :
            done = 0
            if s.numtrialsUsed[bi] == 0 : s.builders[bi].endSession() ; s.builders[bi].startSession()
            for ti in range(s.numtrialsUsed[bi] , s.numtrials[bi]) :
                s.numtrialsUsed[bi] = s.numtrialsUsed[bi] + 1
                print "----------BUILDER index %d %d/%d ----------" % (bi, ti, s.numtrials[bi]),
                s.builders[bi].describe()
                if not s.builders[bi].build(s.pts) : # call builder
                    print "Builder could not build" ; continue
                restraintSatisfied = 1 # check restraints
                for ri in range(len( s.restraints[bi] )) :
                    restraintSatisfied = s.restraints[bi][ri].satisfied(s.pts)
                    if restraintSatisfied == 0 : break
                if restraintSatisfied == 0 :
                    print "restraints unsatisfied" ; continue
                if s.grid.add(s.builders[bi].getOP(), s.clashExclusions[bi]) == 0 : # clash-check
                    print "clashcheck failed" ; continue
                done = 1
                break
            if done == 0 :
                print "fallback of", bi, "=", s.builderFallIndex[bi]
                newbi = s.findFallback(bi)
                for bk in range(bi, newbi, -1) :
                    cbop = s.builders[bk].getOP()
                    for i in range(cbop.size()) : s.grid.remove(cbop[i])
                    s.numtrialsUsed[bk] = 0
                if newbi < 0 : return 0
                bi = newbi
            else : bi = bi + 1
        return 1
