
#from data import vdwr, resAtoms, mcConn, scConn, consts

## verify that builder order given is ok, ie a builder builds upon known or built atoms
def checkBuilderOrder(blist, knownPositions) :
    kp = set(knownPositions)
    for b in blist :
        bop = b.getOP()
        bip = b.getIP()
        for i in range(bip.size()) :
            if not bip[i] in kp :
                print "builder %s does not build upon known/built point %d" % (b.name(), bip[i])
                assert None
        for i in range(bop.size()) :
            kp.add(bop[i])

## detects if a point is output more than a builder, or no builder, or is known still built
def builderSanityCheck(blist, numpts, knownPositions) :
    srcBldr = [-999] * numpts
    for bi in range(len(blist)) :
        b = blist[bi]
        bop = b.getOP()
        for i in range(bop.size()) :
            #print bop[i]
            if srcBldr[ bop[i] ] != -999 :
                print bop[i], "is probably being built twice!", b.name(), blist[srcBldr[ bop[i] ]].name() ; assert 0
            srcBldr[ bop[i] ] = bi
            if bop[i] in knownPositions : 
                print bop[i], 'shdnt be a known point! as it is built by', b.name() ; assert 0 
    for i in range(len(srcBldr)) :
        if srcBldr[i] == -999 :
            if i in knownPositions : pass
            else : print i, "is an unknown point without builder" ; assert(0)

## from the list of parents of nodes, write a dot digraph
def writeDotGraph(parents, nodenames, order=None) :
    dotstr = "digraph { "
    for ni in range(len(parents)) :
        for pi in parents[ni] :
            dotstr = dotstr + ' "%s" -> "%s" ; ' % (nodenames[ni], nodenames[pi])
    if order != None :
        for ni in range( len(nodenames)-1 ) :
            dotstr = dotstr + ' "%s" -> "%s" [color=red] ; ' % (nodenames[ni], nodenames[ni+1])
    dotstr = dotstr + " }"
    return dotstr

# find builder dependences on other builders and builds topo-sort of builders
# execute after builderSanityCheck
from toposort import findTopologicalOrder
def findBuilderOrder(blist, numpts, knownPositions, writeDigraph=None) : # writeDigraph can be filename also
    builderSanityCheck(blist, numpts, knownPositions)
    opOf, ipTo = [-999] * numpts, []
    for i in range(numpts) : ipTo.append([])
    builderDep = {}
    for bi in range(len(blist)) :
        builderDep[bi] = set([])
        b = blist[bi]
        bop, bip = b.getOP(), b.getIP()
        for i in range(bip.size()) : ipTo[ bip[i] ].append(bi)
        for i in range(bop.size()) : opOf[ bop[i] ] = bi
    #print "OPOF", opOf
    #print "IPTO", ipTo
    for i in range( numpts ) :
        for k in range(len(ipTo[i])) :
            if opOf[i] == -999 : continue
            builderDep[ ipTo[i][k] ].add( opOf[i] )
    #print "BUILDERDEPS", builderDep, '\n\n\n'
    topoOrder, topoParent, topoParentIndex = findTopologicalOrder(builderDep)
    if not writeDigraph : return topoOrder, topoParent, topoParentIndex
    digfile = open(writeDigraph, 'w')
    # write digraph and tree in graphviz format
    names = [] ;
    for b in blist : names.append( str(b.name()) )
    print >> digfile, writeDotGraph(builderDep, names)
    #print topoOrder
    names = [] ;
    for bi in topoOrder : names.append( str(blist[bi].name()) )
    #print topoParent
    #print topoParentIndex
    uniqueParent = []
    for t in topoParentIndex :
        if t == -1 : uniqueParent.append( [] )
        else : uniqueParent.append( [t] )
    print >> digfile, writeDotGraph(uniqueParent, names, order=1)
    digfile.close()
    return topoOrder, topoParent, topoParentIndex


def findFallback(bi, builderFallIndex, numtrials, numtrialsUsed) :
    fi = builderFallIndex[bi]
    candidate = -1
    if fi == -1 : # have no-one to fall upon, find earlier bi with fi = -1
        for bk in range(bi-1, -1, -1) :
            if builderFallIndex[bk] == -1 :
                candidate = bk
                break
    else : candidate = fi
    
    if candidate == -1 : return candidate
    if numtrialsUsed[candidate] >= numtrials[candidate] : 
        return findFallback(candidate, builderFallIndex, numtrials, numtrialsUsed)
    return candidate
 
def Build(bi, grid, blist, builderFallIndex, numtrials, numtrialsUsed, rlist, clashExclInds) :
    while bi < len(blist) :
        done = 0
        if numtrialsUsed[bi] == 0 : blist[bi].endSession() ; blist[bi].startSession()
        for ti in range(numtrialsUsed[bi] , numtrials[bi]) :
            numtrialsUsed[bi] = numtrialsUsed[bi] + 1
            print "----------BUILDER index %d %d/%d ----------" % (bi, ti, numtrials[bi]),
            blist[bi].describe()
            if not blist[bi].build() : # call builder
                print "Builder could not build" ; continue
            restraintSatisfied = 1 # check restraints
            for ri in range(len( rlist[bi] )) :
                restraintSatisfied = rlist[bi][ri].satisfied()
                if restraintSatisfied == 0 : break
            if restraintSatisfied == 0 :
                print "restraints unsatisfied" ; continue
            if grid.add(blist[bi].getOP(), clashExclInds[bi]) == 0 : # clash-check
                print "clashcheck failed" ; continue
            done = 1
            break
        if done == 0 :
            print "fallback of", bi, "=", builderFallIndex[bi]
            newbi = findFallback(bi, builderFallIndex, numtrials, numtrialsUsed)
            for bk in range(bi, newbi, -1) :
                cbop = blist[bk].getOP()
                for i in range(cbop.size()) : grid.remove(cbop[i])
                numtrialsUsed[bk] = 0
            if newbi < 0 : return 0
            bi = newbi
        else : bi = bi + 1
    return 1


# check which points are positionally-inited
# determine topo-sorted order of builders
# determine which restraint can be checked after a builder
def findBuilderRestraintOrder(blist, irlist) :
    newrestraints = []
    for i in range(len(blist)) : newrestraints.append( set([]) )
    # find source builder for each point
    ptinds = set() ## find number of points
    for b in blist :
        for i in range(b.getIP().size()) : ptinds.add( b.getIP()[i] )
        for i in range(b.getOP().size()) : ptinds.add( b.getOP()[i] )
    srcBldr = {}
    for k in ptinds : srcBldr[k] = -999
    for bi in range(len(blist)) :
        b = blist[bi]
        bop = b.getOP()
        for i in range(bop.size()) :
            assert srcBldr[ bop[i] ] == -999
            srcBldr[ bop[i] ] = bi
    
    for ri in range(len(irlist)) :
        r = irlist[ri]
        rinds, highestBI = r.getInds(), -999
        for i in range(rinds.size()) :
            bi = srcBldr[ rinds[i] ]
            if highestBI < bi : highestBI = bi
        assert highestBI != -999
        newrestraints[highestBI].add(ri)
    restraints = []
    for riset in newrestraints :
        rs = []
        for ri in riset : rs.append(irlist[ri])
        restraints.append(rs)
    return restraints

