
## without any sidechain atoms
res, resids, resnums, resns, chids, inscodes, pts
aiSC, ai2pos, radii

for ri in resids.keys() :
    for altloc in ["A", "B", "C"] :
        for an in resAtoms [resns[ri]] :
            if an in [N,CA,C,O,CB] : continue
            res[ ri ][ an+altloc ] = len(pts)
            aiSC[ len(pts) ] = 
            ai2pos[ len(pts) ] =
            radii.append( an )
            pts.append( [-999., -999., -999.] )

for ri in resids.keys() :
    ## code for finding multimers
    #############
    for arot in multimers :
        
        

