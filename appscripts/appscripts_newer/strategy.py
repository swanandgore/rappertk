from cbutils import findBuilderOrder, findBuilderRestraintOrder, Build
from geometry import Grid, VecVecInt, VecInt

def makeClashExcl_1st2ndCovNbr(covconn, blist) :
    clashExclInds = []
    print "starting clasExcl--------------------------------------------"
    if len(clashExclInds) == 0 : # 1st & 2nd cov nbrs are to be excluded
        for bi in range(len(blist)) :
            bop = blist[bi].getOP()
            excl = []
            for i in range(bop.size()) :
                ai = bop[i]
                nbrs = set()
                if not ai in covconn.keys() : continue
                for ai1 in covconn[ai] :
                    nbrs.add(ai1)
                    for ai2 in covconn[ai1] : nbrs.add(ai2)
                print "CLASH_EXCL", ai, list(nbrs)
                excl.append( list(nbrs) )
            clashExclInds.append(excl)
    for ci in range(len(clashExclInds)) :
        clashExclInds[ci] = VecVecInt( clashExclInds[ci] )
    print "Auto-generation of clashExclInds done------------------------\n\n\n"
    return clashExclInds

def strategy(blist, pts, radii, overlapReductions, knownPositions, clashExclInds, numtrials, irlist, modelRenderer) :
    builderOrder, builderFallback, builderFallIndex = findBuilderOrder(blist, pts.size(), knownPositions, "dotout")
    print '---------------final order of builders----------------------------------------------'
    for i in range(len(builderOrder)) :
        bi = builderOrder[i]
        print bi, " ",
        blist[bi].describe()
        cbip, cbop = blist[bi].getIP(), blist[bi].getOP()
        print '--- fallback %d %d --- ' % (builderFallback[i], builderFallIndex[i]),
        for i in range(cbip.size()) : print cbip[i],
        print ' ---> ',
        for i in range(cbop.size()) : print cbop[i],
        print ''
    print '------------------------------------------------------------------------------------\n\n\n'
    newbuilders, newclashexcl, newnt, fallback = [], [], [], []
    for bi in builderOrder :
        newbuilders.append( blist[bi] )
        newclashexcl.append( clashExclInds[bi] )
        newnt.append( numtrials[bi] )
    blist, clashExclInds, numtrials = newbuilders, newclashexcl, newnt


    rlist = findBuilderRestraintOrder(blist, irlist)
    print 'RLIST----------------------'
    for rs in rlist :
        print len(rs),
        for r in rs : print r,
        print ''
    print '---------------------------\n\n\n'

    ## start model building
    nattepmpts, nmodels = 1000, 100
    kp = VecInt(knownPositions)
    for nm in range(nattepmpts) :
        print "starting MODEL ", nm
        ## add starting points to grid unconditionally
        grid = Grid(3, pts, radii, overlapReductions)
        grid.justAdd(kp)
        #grid.justAdd(blist[ builderOrder[0] ].getIP())
        ## start building
        numtrialsUsed = [0] * len(numtrials)
        bout = Build(0, grid, blist, builderFallIndex, numtrials, numtrialsUsed, rlist, clashExclInds)
        print "BUILD STATUS", bout
        if bout == 1 :
            modelRenderer.render()
            nmodels -= 1
        if nmodels == 0 : break
