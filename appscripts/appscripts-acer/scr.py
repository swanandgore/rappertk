
        from prepareChainV5 import buildAcb
        import checkProtChainsV4



        unbuiltpts = [] ; guidedSamplingRadius = None; nbuilt = 0 ; all_built = {}; all_unbuilt  = {} ; 
        for wer in range(0,s.natt):
            if wer != 0 :  
                unbuiltpts = [] ; nbuilt = 0 ; all_built = {};        all_unbuilt  = {} ; 
                res, resids, resnums, resns, chids, inscodes, pts =  s.removeDummies(res, resids, resnums, resns, chids, inscodes, pts, dummies)

                
            blist, numtrials, bfoll, rlist, dummies = [], [], {}, [], {} ; bandInfrastructure = []
            badids = []

            for bd in s.badresids :
                print "\n          ========>",bd

            for k ,v in resids.items():
                if v in s.badresids :
                    badids.append(k)
    
    
            if s.mconly :
                res, pts = removeSC(res, pts,res.keys(),badids)
    
    
            print "\nGrouping residues into bands...."
            loops, looptypes = locateRegionsRandomize3(resids, chids, s.badresids,s.modelN2C,s.poor) # locate loops and order them randomly

    
            print "\n\nChecking for missing atoms and chain breaks....."
            print "\n\n=====================  MISSING ATOMS SUMMARY ====================================================="

            mcmiss, scmiss, chainBreaks = checkProtChainsV4.check(res, resids, resnums, resns, chids, inscodes, pts, s.cacaCutoff,s.null,1,s.mconly)
    
            if s.caRad != None :
                bands = checkProtChainsV4.getbands(res, resids, resnums, resns, chids, inscodes, pts, s.cacaCutoff,s.null,1,s.mconly)
    
            print "\n\n================== END MISSING ATOMS SUMMARY ====================================================="
            print
            print
    
            


            chainBreaks = list(chainBreaks) ;         chainBreaks.sort() ;
            wholenewloops = [];        wholenewlooptypes = []
            modloops , modlooptypes ,  notbuilding = [], [], []
            bands.keys().sort()
            for k in range(len(loops)):
                loopstart = loops[k][0] ;             loopend = loops[k][1]
                if s.caRad != None :
                    for dol in range(loops[k][0],loops[k][1]+1):
                        if ' CA ' not in res[dol].keys():
                            print "If using Ca-restraints, all residues in the region that needs to be rebuilt have to contain a Calpha atom"
                            import sys ; sys.exit()
    
                    for bb in bands.values() :
                        bstart = bb[0] ; bstop = bb[1]
                        
                        if s.caRad != None :
                            if loopstart == bstart and (bstop - bstart) < 2 and loopend >= bstop :
                                for xx in range(loops[k][0] ,bstop + 1) : 
                                    notbuilding.append(resids[xx])
                            elif loopstart == bstart and (bstop - bstart) >= 2 and loopend >= bstop :
                                wholenewlooptypes.append("NtoC")
                                wholenewloops.append([bstart, bstop])
    
                            elif loopstart == bstart and  loopend < bstop and  (bstop - bstart > 0):
                                wholenewlooptypes.append("Nter")
                                wholenewloops.append([bstart, loopend])                            

                            elif loopstart == bstart and  loopend <= bstop and  (bstop - bstart == 0):
                                for xx in range(loopstart , loopend+1) : 
                                    notbuilding.append(resids[xx])
    
                            elif loopstart > bstart  and loopend >= bstop and loopstart <= bstop and (bstop - bstart > 0):
                                wholenewlooptypes.append("Cter")
                                wholenewloops.append([loopstart,bstop])                            
    
                            elif loopstart > bstart  and loopend >= bstop and loopstart <= bstop and (bstop - bstart == 0):
                                for xx in range(loopstart , bstop+1) : 
                                    notbuilding.append(resids[xx])
                            
                            elif loopstart > bstart and loopend < bstop:
                                wholenewlooptypes.append("loop")
                                wholenewloops.append([loopstart,loopend])
    
    
                            elif loopstart < bstart and loopend == bstart and (bstop - bstart > 0) :
                                wholenewlooptypes.append("Nter")
                                wholenewloops.append([bstart,bstart])

                                
                            elif loopstart < bstart and loopend == bstart and (bstop - bstart == 0) :
                                for xx in range(bstart , bstart+1) : 
                                    notbuilding.append(resids[xx])
    
                            elif loopstart < bstart and loopend < bstop and loopend in range(bstart,bstop) and (bstop - bstart > 0):
                                wholenewlooptypes.append("Nter")
                                wholenewloops.append([bstart,loopend])
                                
                            elif loopstart < bstart and loopend < bstop and loopend in range(bstart,bstop) and (bstop - bstart == 0):
                                for xx in range(bstart , loopend+1) : 
                                    notbuilding.append(resids[xx])
                                    
                            elif loopstart < bstart and loopend >= bstop and bstop-bstart >= 2:
                                wholenewlooptypes.append("NtoC")
                                wholenewloops.append([bstart,bstop])

                            elif loopstart < bstart and loopend >= bstop and bstop-bstart < 2:
                                for xx in range(bstart ,bstop + 1) : 
                                    notbuilding.append(resids[xx]);

                                    
                            else :
                                print "Skipping ... "
                                
                            
                        print "--------------------------------------------------------------------------------------"
                        print
                        print
    
    
            if s.caRad != None : 
                loops = wholenewloops
                looptypes = wholenewlooptypes
    
    
            min_band_size = 5 ; max_band_size = 50
    
            for kx in range(len(loops)):
                if loops[kx][1] - loops[kx][0] > 500  and (looptypes[kx] == "NtoC" or looptypes[kx]== 'Cter' or looptypes[kx] == 'Nter') and s.caRad != None :
                    
                    modloops,modlooptypes = break_chain_in_bands( loops[kx][0], loops[kx][1], pts, min_band_size ,max_band_size ,res,resns,resids)
    
                    loops.remove(loops[kx])
                    looptypes.remove(looptypes[kx])
                    for ml in modloops:
                        loops.append(ml)
                    for mlt in modlooptypes :
                        looptypes.append(mlt)
                else :
                    if loops[kx][1] - loops[kx][0] > 500  and s.caRad == None and looptypes[kx] != "loop":
                        print "Rappertk cannot build a segment greater than 500 residues"
                        print "Either (1) set --use-ca-restraints True "
                        print "                  or                    "
                        print "       (2) Restrict residue range to < 500 residues"
                        
                        import sys ; sys.exit()
                        
            
    
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
            
            for xd in s.missingpts.keys():
                missid = rev_resids[xd]
                for aid in s.missingpts[xd]:
                    ptid  = res[missid][aid]
                    missingpointids.append(ptid)
                    
                
            for li in range(len(loops)) :
                print "BANDS ARE ",resids[loops[li][0]],resids[loops[li][1]],looptypes[li]
                startindex, endindex = loops[li]
                if s.scRad == None:
                    for x in range(startindex,endindex+1):
                        scMissInds.append(x)


            

                        
            for li in range(len(loops)) :
                startindex, endindex = loops[li]
                chid = chids[startindex]
                assert chid == chids[endindex]
                bl,bf,nt,rl,dum, reorderBuilders = None,None,None,None,None, 1
    
                
                if s.caRad == None : 
                    for x in range(startindex,endindex+1):
                        caMissInds.append(x)
    
                        
                if s.scRad == None : 
                    s.scRad  = 1. 
    
    
                        
    
                if looptypes[li] == "NtoC" :
                    if resns[ loops[li][0] ] in [ "  A", "  T", "  C", "  G", "  U", ] :
                        bl,nt,bf,rl,dum = prepareRNAchain(chid, res, resids, resnums, resns, chids, inscodes, pts, s.caRad, s.scRad)
                    else :
                
                        rev_resnum = {} ; chainkey = []
                        chainid = chids[startindex]

                        for k,v in chids.items():
                            if chids[k] == chainid: chainkey.append(k)

                        for k in chainkey : rev_resnum[int(resnums[k])] = int(k)
    

                        startResidue = int(resnums[startindex]) ; endResidue = int(resnums[endindex])


                        
                        if (startindex not in chainkey) :
                            print "(1) N-ter residue [%s] for building chain needs to be present in input coordinate file"%resids[startindex],s.null
                            import sys  ; sys.exit()
    
                        if ' CA ' not in res[startindex].keys():
                            print "C-alpha atom of anchor N residue [%s] for chain to be modelled needs to be present in the pdbfile"%resids[startindex]
                            import sys  ; sys.exit()
    
    
    
                        bl,nt,bf,rl,dum = s.prepC.preparePeptideChain(startindex,endindex,chid, res, resids, resnums, resns, chids, inscodes, pts, s.mconly, 1, guidedSamplingRadius, s.caRad, s.scRad, scMissInds,caMissInds)
    
                elif looptypes[li] == "Nter" :
    
                    chainkey = [] ;  chainid = chids[startindex]
                    for k,v in chids.items():
                        if chids[k] == chainid:
                            chainkey.append(k)



                    if endindex+1 not in chainkey : 
                        print " Anchor residue for %s  needs to be present in the coordinate file"%resids[endindex]
                        import sys ; sys.exit()
                        
                        
    
                    if " CA " not in res[endindex+1].keys() :
                        print " C-alpha atom of anchor residue %s  needs to be present in the coordinate file"%resids[endindex+1]
                        
                        import sys ;  sys.exit()
                        
                    dum, firstindex = addNdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts)
    
    



                    bl,bf,nt,rl = s.prepC.prepareChainTerminal("Nterm", startindex, endindex, firstindex, res, resids, resnums, resns, chids, inscodes, pts, s.mconly, s.caRad, s.scRad, guidedSamplingRadius, scMissInds,caMissInds)
    



    
                elif looptypes[li] == "Cter" :
                    chainkey = []
                    chainid = chids[startindex]
                    for k,v in chids.items():
                        if chids[k] == chainid:
                            chainkey.append(k)

    


                    if startindex-1 not in chainkey :
                        print "Anchor residue for  %s  C-terminal to be modelled needs to be present in the pdbfile"%resids[startindex]
                        import sys ;  sys.exit()
    
    
                    if " CA " not in res[startindex-1].keys():
                        print "CA atom of Anchor  C residue  for C-terminal to be modelled needs to be present in the pdbfile" ; import sys ;  sys.exit()
    
                        
                    dum, lastindex = addCdummyGly(chid, res, resids, resnums, resns, chids, inscodes, pts)
                    bl,bf,nt,rl = s.prepC.prepareChainTerminal("Cterm", startindex, endindex, lastindex, res, resids, resnums, resns, chids, inscodes, pts, s.mconly, s.caRad, s.scRad, guidedSamplingRadius, scMissInds,caMissInds)
    
    
    
                elif looptypes[li] == "loop" :
                    chainkey = []
                    chainid = chids[startindex]
                    for k,v in chids.items():
                        if chids[k] == chainid:
                            chainkey.append(k)

                        
                    if (startindex-1 not in chainkey):
                        print "Anchor N residue for loop to be modelled needs to be present in the pdbfile"
                        import sys ;  sys.exit()
    
                    if ((startindex-2) not in chainkey ):
                        print "Anchor N-terminal residues,start-1 and start-2,  for loop to be modelled needs to be present in the pdbfile"
                        import sys ;  sys.exit()                    
            
                    if ((endindex + 1 ) not in chainkey ):
                        print "Anchor C-terminal residue for loop to be modelled needs to be present in the pdbfile"
                        import sys ; sys.exit()
    
                    if ((endindex + 2 ) not in chainkey ):
                        print "Anchor C-terminal residue for loop to be modelled needs to be present in the pdbfile"
                        import sys ; sys.exit()                    
                        
    
                    if " CA " not in res[startindex-1].keys() or  " CA " not in res[endindex+1].keys():
                        print "CA atom of Anchor N + C residue  for loop to be modelled needs to be present in the pdbfile"
                        import sys ; sys.exit()
    
                    

                    if endindex-startindex <= 5  :   #or random.random() > 1.5 :
    
                        bl,bf,nt,rl = s.prepC.preparePeptideLoop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, s.mconly, s.caRad, s.scRad, scMissInds, guidedSamplingRadius,None,caMissInds,s.loopclosure)
    
    
                    else :
                        
                        bl,bf,nt,rl = s.prepC.preparePeptideMidloop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, s.mconly, s.caRad, s.scRad, scMissInds, guidedSamplingRadius,caMissInds,s.loopclosure)
                        
                        reorderBuilders = None
                
    
                else :
                    print "Dont understand looptype", looptypes[li] ; assert None
    
    

                mergeBRlists(blist, bfoll, numtrials, rlist, bl, bf, nt, rl, dummies, dum)
                bandInfrastructure.append( (bl,nt,bf,[], reorderBuilders) )

    


        
    
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
                if s.edopt == 1 : 
                    for ori in range(len(xrlist)) :
                        extraRestrOpt.append(rsize + ori)
    



            for xr in rlist :
                extraRestraints.append(xr)

            if s.allopt   == 1 :
                rsize = len(extraRestraints)
                for ori in range(len(rlist)) :
                    extraRestrOpt.append(rsize + ori)
    

            modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, dummies.keys(), s.outpdb)

    

            while (nbuilt < s.nmodels) :

                all_built[nbuilt] = [] ;              all_unbuilt[nbuilt] = []
                built, unbuilt,unbuiltpts1 = bandBuild(loops, looptypes, bandInfrastructure, res, resids, resnums, resns, chids, inscodes, pts, extraRestraints, extraRestrOpt,    s.backtrack, s.popsize, radii, gridHelper, s.ranker,missingpointids,s.natt)
                for r in unbuiltpts1 :
                   unbuiltpts.append(r)
                
                if len(unbuiltpts1) == 0 :
                    all_built[nbuilt].append(built);        all_unbuilt[nbuilt].append(unbuilt)
                    modelRenderer.render(pts,None,missingpointids,unbuiltpts) ;
                    nbuilt += 1
                else:
                    all_built[nbuilt].append(built);        all_unbuilt[nbuilt].append(unbuilt)
                    modelRenderer.render(pts,None,missingpointids,unbuiltpts) ;
                    nbuilt += 1

                
            if nbuilt == s.nmodels :
                return nbuilt,all_built,all_unbuilt

        return nbuilt,all_built,all_unbuilt





###        else :
            freeRlabel = None

2..... 237
2..... 238
2..... 239
2..... 240
2..... 233
