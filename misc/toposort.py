import random, math

# given a digraph in form of adj list of incoming and outgoing edges, topologically sort
# mychildren and myparents are lists of sets
# a -> b means b depends on a. we want a to be earlier to b in sorted order
# while 1:
#  find all parentless unvisited node. if cant, break
#  mark all visited with a visit-order
#  find its children with all parents visited and mark them as endpoints
def topoSortBFS(myparents, mychildren) :
    assert len(myparents) == len(mychildren)
    visited = [-999] * len(myparents)
    visitOrder = 0
    layers = []
    print "myparents", myparents
    print "mychildren", mychildren
    while 1 :
	candidates = []
	for i in range(len(myparents)) :
	    if visited[i] != -999 : continue
	    allParentsVisited = 1
	    for pi in myparents[i] :
		if visited[pi] == -999 : #unvisited parent
		    allParentsVisited = 0
		    break
	    if allParentsVisited == 0 : continue
	    candidates.append(i)
	print "candidates" , candidates
	if len(candidates) == 0 : break # done!
	# sort candidates on descending number of parents
	for ci0 in range(len(candidates)) :
	    for ci1 in range(ci0+1,len(candidates)) :
		if myparents[ candidates[ci0] ] < myparents[ candidates[ci0] ] :
		    temp = candidates[ci0]
		    candidates[ci0] = candidates[ci1]
		    candidates[ci1] = temp
	print "sorted candidates" , candidates
	layers.append(candidates)
	for ci in candidates :
	    visited[ci] = visitOrder
	    visitOrder = visitOrder + 1
    return layers

# given a digraph in form of adj list of incoming and outgoing edges, topologically sort
# mychildren and myparents are lists of sets
# a -> b means b depends on a. we want a to be earlier to b in sorted order
# while 1:
#  find a parentless unvisited node. if cant, break
#  mark it visited with a visit-order
#  find its children with all parents visited and mark them as endpoints
def topoSortDFS(myparents, mychildren, randomized=None) :
    assert len(myparents) == len(mychildren)
    visited = [-999] * len(myparents)
    visitOrder = 0
    parentless = [0] * len(myparents)
    for i in range(len(myparents)) :
	if len(myparents[i]) == 0 : parentless[i] = 1
    topoOrder = []
    while 1 :
	candidates = []
	for i in range(len(myparents)) :
	    if parentless[i] == 1 and visited[i] == -999 :
		candidates.append(i)
	if len(candidates) == 0 : break # done!
	if randomized != None :
	    rind = (int) (math.floor( random.random() * len(candidates) ))
	    curi = candidates[rind]
	else : curi = candidates[0]
	topoOrder.append(curi)
	visited[curi] = visitOrder
	visitOrder = visitOrder + 1
	for ci in mychildren[curi] :
	    allParentsVisited = 1
	    for k in myparents[ci] :
		if visited[k] == -999 : allParentsVisited = 0
	    if allParentsVisited == 1 : parentless[ci] = 1
    print topoOrder
    return topoOrder

def findTreeSize(ni, treesizes, newchildren) :
    if treesizes[ni] >= 0 : return
    treesizes[ni] = 1
    for ci in newchildren[ni] :
        findTreeSize(ci, treesizes, newchildren)
        treesizes[ni] = treesizes[ni] + treesizes[ci]

def findTopologicalOrder(myparents, bwt) :
    # make mychildren from myparents
    mychildren = []
    for i in range(len(myparents)) : mychildren.append( set() )
    for i in range(len(myparents)) :
        for k in myparents[i] : mychildren[k].add(i)
    #print "myparents", myparents
    #print "mychildren", mychildren
    # use dfs and get assign 0/1 parent to all nodes
    visited = [0] * len(mychildren)
    newchildren = []
    for i in mychildren : newchildren.append( set() )
    stack = [] #push parentless nodes into stack & call dfs
    for i in range(len(myparents)) :
        if len(myparents[i]) == 0 : stack.append(i)
    newDFSchildren(myparents, mychildren, newchildren, visited, stack)
    #print "newchildren", newchildren
    # newchildren --> newparents
    newparents, newp = [], []
    for i in range(len(newchildren)) : newparents.append( set() )
    for i in range(len(newchildren)) :
        for k in newchildren[i] : newparents[k].add(i)
    for i in range(len(newchildren)) :
        if len(newparents[i]) == 0 : newp.append(-1)
        elif len(newparents[i]) == 1 : newp.append( list(newparents[i])[0] )
        else : assert 1==0 # each node has to have 0/1 parent
    newparents = newp
    #print "newparents", newparents
    # find size of tree rooted at each node
    treesizes = list(bwt)
    for ni in range(len(newchildren)) : findTreeSize(ni, treesizes, newchildren)
    # find sorted order, smallest tree first
    stack, topoOrder = [], []
    for ni in range(len(newchildren)) :
        if len(myparents[ni]) == 0 : stack.append(ni)
    for i in range(len(stack)) :
        for j in range(i+1,len(stack)) :
            if treesizes[ stack[i] ] < treesizes[ stack[j] ] :
                temp = stack[i]
                stack[i] = stack[j]
                stack[j] = temp
    while 1 :
        if len(stack) == 0 : break
        ni = stack.pop()
        topoOrder.append(ni)
        children = list(newchildren[ni])
        for i in range(len(children)) :
            for j in range(i+1,len(children)) :
                if treesizes[ children[i] ] < treesizes[ children[j] ] :
                    temp = children[i]
                    children[i] = children[j]
                    children[j] = temp
        stack = stack + children
    topoParent = []
    for t in topoOrder : topoParent.append( newparents[t] )
    topoParentIndex = [-1] * len(topoParent) # topoParentIndex
    for i in range( len(topoParent)-1, -1, -1 ) :
        par = topoParent[i]
        parind = -1
        try : parind = topoOrder.index(par)
        except ValueError : pass
        if parind == -1 :
            for k in range( i-1, -1, -1 ) :
                if topoParent[k] == -1 :
                    parind = k; break
        topoParentIndex[i] = parind
    return topoOrder, topoParent, topoParentIndex

def newDFSchildren(myparents, mychildren, newchildren, visited, stack) :
    # pop a node from stack & mark it visited until u find node with >0 unvisited children
    goodni = -1
    while 1 :
        #print "stack", stack, len(stack)
        if len(stack) == 0 : break
        ni = stack.pop()
        visited[ni] = 1
        #print "popping", ni, " has children" , mychildren[ni]
        goodni = -1
        for ci in mychildren[ni] :
            if visited[ci] == 0 : goodni = ni ; break
        if goodni >= 0 : break
    if goodni >= 0 : ni = goodni
    else : return
    # for an unvisited child in its children,
    #   if its parents are all visited, push it on stack, mark it in newchildren & recurse
    for ci in mychildren[ni] :
        if visited[ci] == 1 : continue
        allParentsVisited = 1
        assert len(myparents[ci]) > 0
        for pi in myparents[ci] :
            if visited[pi] == 0 : allParentsVisited = 0 ; break
        if allParentsVisited == 0 : continue
        stack.append(ci)
        newchildren[ni].add(ci)
        newDFSchildren(myparents, mychildren, newchildren, visited, stack)
    #find a node with no parent or all visited parents
    #if none found, return
    #mark it visited
    #for each of its child,
    #    mark it visited
    

## keep nodepairs together if possible in topological sorted order
def topoSort(parents, nodepairs, randomized=None) :
    children = []
    for i in range(len(parents)) : children.append( set() )
    for i in range(len(parents)) :
	for k in parents[i] :
	    children[k].add(i)
    return newTopoSort(parents, children)
    #print "CHILDREN", children
    layers = topoSortBFS(parents, children)
    topoOrder = []
    for li in range(len(layers)-1) :
	thislayer, nextlayer = layers[li], layers[li+1]
	for np in nodepairs :
	    if not np[0] in thislayer : continue
	    if not np[1] in nextlayer : continue
	    thislayer.remove(np[0])
	    thislayer.append(np[0])
	    nextlayer.remove(np[1])
	    nextlayer.insert(0,np[1])
	    nodepairs.remove(np)
	    break
	for ti in thislayer : topoOrder.append(ti)
    for ti in layers[ len(layers)-1 ] : topoOrder.append(ti)
    return topoOrder
    #return topoSortDFS(parents, children, randomized)

if __name__ == '__main__' :
    parents, children = [], []

    parents.append( set([]) )
    parents.append( set([0]) )
    parents.append( set([1]) )
    parents.append( set([1]) )
    parents.append( set([0]) )
    parents.append( set([4]) )
    parents.append( set([]) )
    parents.append( set([6]) )
    parents.append( set([6]) )

    topoOrder, topoParent, topoParentIndex = findTopologicalOrder( parents )
    assert len(topoOrder) == len(topoParent)
    assert len(topoOrder) == len(topoParentIndex)
    for i in range(len(topoOrder)) :
	print topoOrder[i], topoParent[i], topoParentIndex[i]
