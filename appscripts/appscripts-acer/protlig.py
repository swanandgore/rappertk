from peptidebuild import ModelRenderer, BufferModelRenderer
from pdbr import protein
import prot2res
import re, string, sys, os, shutil
from geometry import calcDist, CAtraceGH, MapIntIntFloat
from builders import VecInt, VecFloat, VecBuilder, VecVecFloat, Rotator, TransRotator, BuilderGroup
from restraints import DistanceRestraint
from data import vdwr, consts, PROBE_DISULFIDE_OVERLAP_MARGIN
from prepareChain import  incompleteSCcorrection
import data
from data import vdwr, consts
from data import sgtable , long2shortHM
import random
from checkfornone import fixdirnames
import prepareChainV5
from PopulationStrategy import PopulationStrategy
from xray import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS, fft, omitmap, mapman, mapdump , mtzdump , fftnew
from prepareChainV5 import removeSC, mergeBRlists, addNdummyGly, addCdummyGly, printResAtoms, PrepareChain, prepareRNAchain, incompleteSCcorrection , incompleteMCSCcorrection  , incompleteMCcorrection , incompleteSCcorrection2
from loopbuild12 import Multiloop, incompleteSCcorrection2 ,  locateRegionsRandomize3

def findLigDists(ligname, res, pts, resns,hetind) :
    ligdists = {}
    ligind = None
    print "ligname",ligname,len(ligname)
#    import sys ; sys.exit()
    print res.keys()
    for ind in res.keys() :
        if resns[ind] == ligname : ligind = ind ; break
        print resns[ind]
    if ligind == None:
        print "%s not found in ip PDB file"%ligname
        sys.exit()
    for ind in res.keys() :
        if resns[ind] == ligname or 'HOH' in resns[ind] or hetind[ind] == 1:
            print "Continuing for" , ligname
            continue
        for name,ai in res[ind].items() :
            mindist = 1e10
            for lan,lai in res[ligind].items() :
                ldist = calcDist( VecFloat(pts[lai]), VecFloat(pts[ai]) )
                #print "ldist", lan, ldist, pts[lai], pts[ai]
                if ldist < mindist : mindist = ldist
            if mindist < 1e9 :

                ligdists[ind] = mindist
    return ligdists

def parseAtomnames(l) :
    '''from line containing tokens of type [aaaa], extract 4-character atomnames and return the list'''
    i = 0 ; anames = []
    while i <= len(l)-4 :
        if l[i] == '[' :
            assert l[i+5] == ']'
            anames.append( l[i+1:i+5] )
            i = i + 6
        else : i = i + 1
    return anames

def makeLigBuilder(res, resns, pts, ligfile) :
    '''make a series of builders and group into a GroupBuilder.
        Currently recongnize only the lines starting with mindist, rotbond and init'''
    ligind, restraints, builders = None, [], []
    tol1, tol2, init1, init2 = None, None, None, None
    allAnames, ligname = set(), None
    for l in open(ligfile, 'r').readlines() :
        l = re.sub("\n", "", l)
        l = re.sub(" *#.*", "", l)
        if len(l) == 0 : continue
        flds = l.split()
        anames = parseAtomnames(l)
        for an in anames : allAnames.add(an)
        if 'mindist' in l:
            for fi in range(1,len(flds)) :
                if flds[fi][0] == '[' : break
                tol = string.atof(flds[fi])
                restraints.append( DistanceRestraint(VecInt([res[ligind][anames[0]], res[ligind][anames[1]]]), "DistanceRestraint within ligand", tol, 1e10) )
        elif  'ligname' in l :
            for k,v in resns.items() :
                if v == flds[1] : ligname = flds[1] ; ligind = k ; break
            if ligind == None:
                print "ligname entry not found in ",ligfile
                import sys; sys.exit()
            assert ligind != None
        elif  'rotbond'  in l :
            min, max, step = string.atof(flds[1]), string.atof(flds[2]), string.atof(flds[3]),
            IPs, OPs = [ res[ligind][anames[0]],res[ligind][anames[1]] ], []
            for an in anames[2:] : OPs.append(res[ligind][an])
            builders.append( Rotator(VecInt(IPs), VecInt(OPs), None, "Rotator", min, max, step) )
        elif 'init' in l :
            tol1, tol2 = string.atof(flds[1]), string.atof(flds[2])
            init1, init2 = anames[0], anames[1]
        else :
            print "Unrecognized file format for %s, exiting....." % ligfile ;
            import sys; sys.exit()
            assert None
    assert tol1
    IPs, OPs = [], [res[ligind][init1], res[ligind][init2]]
    allAnames.remove(init1) ; allAnames.remove(init2)
    for an in allAnames : OPs.append(res[ligind][an])
    builders = [ TransRotator(VecInt(IPs),VecInt(OPs),None,"TransRotator", VecFloat(pts[res[ligind][init1]]),tol1, VecFloat(pts[res[ligind][init2]]),tol2) ] + builders
    return builders, restraints, ligname



def main(pdbfile,pdbout,dir_xyzout,mapfn,mtzfn,caRes,scRes,caRad,scRad,scReduction,popsize,verbose,backtrack,rotLib,nmodels,mconly,sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freeRlabel,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt,mapformat,ranker, xrayRestGen, ligfile , closeCutoff):

    prot = protein(pdbfile, read_hydrogens=0, read_waters= 0, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts, hetind = prot2res.readProtRes3(prot)
    null = []
    
    ## find ligand name
    ligname = ""
    for l in open(ligfile,'r').readlines() :
        if l[0:7] == "ligname" :
            for xl in l[8:]:
                if xl in [" ","\t","\n"] :
                    continue
                else:

                    ligname = ligname + xl
            
    if ligname == "" or ligname == None :
        print "Cannot find ligand in coordinate file"
        import sys ; sys.exit()




    leastDists = findLigDists(ligname, res, pts, resns,hetind)
    closeRes = [None] * (len(leastDists))

    for ind,d in leastDists.items() :
        if d < closeCutoff :
            closeRes[ind] = 1

    segstart, segs = None, []
    for ind in range(len(closeRes)) :
        if closeRes[ind] and not segstart : segstart = ind
        elif closeRes[ind] == None and segstart :
            segs.append([segstart,ind-1])
            segstart = None
    if segstart : segs.append([segstart,len(closeRes)-1])

    ## correct segments that are separated by 1 residue

    segs2remove = []
    for si in range(1,len(segs)) :
        if segs[si][0] - segs[si-1][1] <= 2 :
            if segs[si][1]-segs[si][0] == 0 and segs[si-1][1]-segs[si-1][0] == 0 : # both lengths 1, merge segs
                segs[si][0] = segs[si-1][0]

                segs2remove.append(si-1)
                si += 1
            elif segs[si][1]-segs[si][0] == 0 : # shrink si-1
                segs[si-1][1] -= 1
            elif segs[si-1][1]-segs[si-1][0] == 0 : # shrink si
                segs[si][0] += 1
            else : segs[si][0] += 1
    newsegs = []
    for si in range(len(segs)) :
        if si in segs2remove : continue
        newsegs.append(segs[si])

    segs = newsegs
    badresids = []
    badids = [] 
    for seg in segs:
        for i in range(seg[0],seg[1]+1) :
            badresids.append(resids[i])
            badids.append(i)
    ### added
    from prepareChainV5 import buildAcb
    import checkProtChainsV4

    unbuiltpts = [] ; guidedSamplingRadius = None; nbuilt = 0 ;  
    blist, numtrials, bfoll, rlist, dummies = [], [], {}, [], {}
    knownPositions = range(len(pts))

    bldrs, irlist, ligname = makeLigBuilder(res, resns, pts, ligfile)
    blist = [ BuilderGroup(VecBuilder(bldrs), "LigandBuilder") ]
    for ind in res.keys() :
        if resns[ind] != ligname : continue
        for an,pi in res[ind].items() : knownPositions.remove(pi)
    numtrials.append(1000)

    
    
    if mconly :
        res, pts = removeSC(res, pts,res.keys(),badids)

    print "\nGrouping residues into bands...."
    loops, looptypes = locateRegionsRandomize3(resids, chids, badresids,"False","False") # locate loops and order them randomly
    print "\n\nChecking for missing atoms and chain breaks....."
    print "\n\n=====================  MISSING ATOMS SUMMARY ====================================================="
    mcmiss, scmiss, chainBreaks = checkProtChainsV4.check(res, resids, resnums, resns, chids, inscodes, pts, cacaCutoff,null,1,mconly)
    if caRad != None :
        bands = checkProtChainsV4.getbands(res, resids, resnums, resns, chids, inscodes, pts, cacaCutoff,null,1,mconly)
    print "\n\n================== END MISSING ATOMS SUMMARY ====================================================="
    print
    print

    prepC = prepareChainV5.PrepareChain(rotLib)
    chainBreaks = list(chainBreaks) ;         chainBreaks.sort() ;
    wholenewloops = [];        wholenewlooptypes = []
    modloops , modlooptypes ,  notbuilding = [], [], []
    bands.keys().sort()
    for k in range(len(loops)):
        loopstart = loops[k][0] ;             loopend = loops[k][1]
        if caRad != None :
            for dol in range(loops[k][0],loops[k][1]+1):
                if ' CA ' not in res[dol].keys():
                    print "If using Ca-restraints, all residues in the region that needs to be rebuilt have to contain a Calpha atom"
                    import sys ; sys.exit()
            for bb in bands.values() :
                bstart = bb[0] ; bstop = bb[1]
                
                if caRad != None :
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
    if caRad != None : 
        loops = wholenewloops
        looptypes = wholenewlooptypes
    min_band_size = 5 ; max_band_size = 50
    for kx in range(len(loops)):
        if loops[kx][1] - loops[kx][0] > 500  and (looptypes[kx] == "NtoC" or looptypes[kx]== 'Cter' or looptypes[kx] == 'Nter') and caRad != None :
            
            modloops,modlooptypes = break_chain_in_bands( loops[kx][0], loops[kx][1], pts, min_band_size ,max_band_size ,res,resns,resids)
            loops.remove(loops[kx])
            looptypes.remove(looptypes[kx])
            for ml in modloops:
                loops.append(ml)
            for mlt in modlooptypes :
                looptypes.append(mlt)
        else :
            if loops[kx][1] - loops[kx][0] > 500  and caRad == None and looptypes[kx] != "loop":
                print "Rappertk cannot build a segment greater than 500 residues"
                print "Either (1) set --use-ca-restraints True "
                print "                  or                    "
                print "       (2) Restrict residue range to < 500 residues"
                
                import sys ; sys.exit()
                
    
    mcMissInds = []; scMissInds = [] ; caMissInds = []; badresids = []  ; badids = [] ; missingpts = {}
    mcmiss, scmiss, chainBreaks1 = checkProtChainsV4.check(res, resids, resnums, resns, chids, inscodes, pts, cacaCutoff,[],0)
    
    for startindex,endindex in loops :
        for ri in range(startindex, endindex+1) :
            badresids.append( resids[ri] )
            badids.append(ri)
            
    for bd in badids :
        if resids[bd] in notbuilding:
            print "Error in banding", notbuilding,resids[bd]
            import sys ; sys.exit()
    if mconly != 1:
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
                if km in missingpts.keys() :
                    for vvm in vm : 
                        missingpts[km].append(vvm)
                else :
                    missingpts[km] = []
                    for vvm in vm : 
                        missingpts[km].append(vvm)
                    
    if mconly != 1:
        for ri in scMissInds :
            misspts = incompleteSCcorrection2(res, resns, pts ,ri,resids)
            for km,vm in misspts.items() :
                if km in missingpts.keys() :
                    for vvm in vm : 
                        missingpts[km].append(vvm)
                else :
                    missingpts[km] = []
                    for vvm in vm : 
                        missingpts[km].append(vvm)
                        
            if ri not in mcMissInds : 
                print "Building temporary sidechain for",resids[sc]
                buildAcb(res, resids, resnums, resns, chids, inscodes, pts, ri)
        vvpts = VecVecFloat(pts)
        for ri in scMissInds :
            prepC.makeChiBuilder(ri,ri,ri, res,resns,resids).buildSample(vvpts, 0)
        for pi in range(len(pts)) : pts[pi] = [ vvpts[pi][0], vvpts[pi][1], vvpts[pi][2], ]
        for r in range(len(loops)):
            if looptypes[r] == "Nter":
                for sc in scmiss :
                    if sc == int(loops[r][1]) + 1  :
                        misspts = incompleteSCcorrection2(res, resns, pts ,sc,resids)
                        scMissInds.append(sc)
                        for km,vm in misspts.items() :
                            if km in missingpts.keys() :
                                for vvm in vm : 
                                    missingpts[km].append(vvm)
                            else :
                                missingpts[km] = []
                                for vvm in vm : 
                                    missingpts[km].append(vvm)
                                    
                            
                        if sc not in mcmiss : 
                            print "Building temporary sidechain for",resids[sc] #,sc,mcMissInds
                            buildAcb(res, resids, resnums, resns, chids, inscodes, pts, sc)
                        vvpts = VecVecFloat(pts)
                        prepC.makeChiBuilder(sc,sc,sc, res,resns,resids).buildSample(vvpts, 0)
                        for pi in range(len(pts)) :pts[pi] = [ vvpts[pi][0], vvpts[pi][1], vvpts[pi][2], ]
                        print "Addining missing sidechain atoms:",sc
                        
            if looptypes[r] == "Cter":
                for sc in scmiss :
                    if sc == int(loops[r][0]) - 1  :
                        misspts = incompleteSCcorrection2(res, resns, pts ,sc,resids)
                        scMissInds.append(sc)
                        for km,vm in misspts.items() :
                            if km in missingpts.keys() :
                                for vvm in vm : 
                                    missingpts[km].append(vvm)
                            else :
                                missingpts[km] = []
                                for vvm in vm : 
                                    missingpts[km].append(vvm)
                                    
                        if sc not in mcmiss: 
                            print "Building temporary sidechain for",resids[sc]
                            buildAcb(res, resids, resnums, resns, chids, inscodes, pts, sc)
                            #                            vvpts = VecVecFloat(pts)                                
                            #                          prepC.makeChiBuilder(sc,sc,sc, res,resns,resids).buildSample(vvpts, 0)
                            #                           for pi in range(len(pts)) : pts[pi] = [ vvpts[pi][0], vvpts[pi][1], vvpts[pi][2], ]
            if looptypes[r] == "loop":
                for sc in scmiss :
                    if sc == int(loops[r][0]) - 1  :
                        misspts = incompleteSCcorrection2(res, resns, pts ,sc,resids)
                        scMissInds.append(sc)
                        for km,vm in misspts.items() :
                            if km in missingpts.keys() :
                                for vvm in vm : 
                                    missingpts[km].append(vvm)
                            else :
                                missingpts[km] = []
                                for vvm in vm : 
                                    missingpts[km].append(vvm)
                                    
                        if sc not in mcmiss : 
                            print "Building temporary sidechain for",resids[sc]
                            buildAcb(res, resids, resnums, resns, chids, inscodes, pts, sc)
                        vvpts = VecVecFloat(pts)
                        prepC.makeChiBuilder(sc,sc,sc, res,resns,resids).buildSample(vvpts, 0)
                        for pi in range(len(pts)) :    pts[pi] = [ vvpts[pi][0], vvpts[pi][1], vvpts[pi][2], ]
                    if sc == int(loops[r][1]) + 1  :
                        misspts = incompleteSCcorrection2(res, resns, pts ,sc,resids) 
                        scMissInds.append(sc)
                        for km,vm in misspts.items() :
                            if km in missingpts.keys() :
                                for vvm in vm : 
                                    missingpts[km].append(vvm)
                            else :
                                missingpts[km] = []
                                for vvm in vm : 
                                    missingpts[km].append(vvm)
                        if sc not in mcmiss : 
                            print "Building temporary sidechain for",resids[sc]
                            buildAcb(res, resids, resnums, resns, chids, inscodes, pts, sc)
                        ###vvpts = VecVecFloat(pts)
                        #s.prepC.makeChiBuilder(sc,sc,sc, res,resns,resids).buildSample(vvpts, 0)
                        #for pi in range(len(pts)) : pts[pi] = [ vvpts[pi][0], vvpts[pi][1], vvpts[pi][2], ]                            
                        
                        
                        ## print to check: What abt restraints
            
    rev_resids = {} ; missingpointids = []
    for king, queen in resids.items():
        rev_resids[queen] = king
    keyt = missingpts.keys()
    
    for xd in missingpts.keys():
        missid = rev_resids[xd]
        for aid in missingpts[xd]:
            ptid  = res[missid][aid]
            missingpointids.append(ptid)
            
        
    for li in range(len(loops)) :
        print "BANDS ARE ",resids[loops[li][0]],resids[loops[li][1]],looptypes[li]
        startindex, endindex = loops[li]
        if scRad == None:
            for x in range(startindex,endindex+1):
                scMissInds.append(x)
    
                
    for li in range(len(loops)) :
        startindex, endindex = loops[li]
        chid = chids[startindex]
        assert chid == chids[endindex]
        bl,bf,nt,rl,dum, reorderBuilders = None,None,None,None,None, 1
        
        if caRad == None : 
            for x in range(startindex,endindex+1):
                caMissInds.append(x)
                
        if scRad == None : 
            scRad  = 1. 
                
        if looptypes[li] == "NtoC" :
            if resns[ loops[li][0] ] in [ "  A", "  T", "  C", "  G", "  U", ] :
                bl,nt,bf,rl,dum = prepareRNAchain(chid, res, resids, resnums, resns, chids, inscodes, pts, caRad, scRad)
            else :
        
                rev_resnum = {} ; chainkey = []
                chainid = chids[startindex]
                for k,v in chids.items():
                    if chids[k] == chainid: chainkey.append(k)
                for k in chainkey : rev_resnum[int(resnums[k])] = int(k)
                startResidue = int(resnums[startindex]) ; endResidue = int(resnums[endindex])
                
                if (startindex not in chainkey) :
                    print "(1) N-ter residue [%s] for building chain needs to be present in input coordinate file"%resids[startindex],null
                    import sys  ; sys.exit()
                if ' CA ' not in res[startindex].keys():
                    print "C-alpha atom of anchor N residue [%s] for chain to be modelled needs to be present in the pdbfile"%resids[startindex]
                    import sys  ; sys.exit()
                bl,nt,bf,rl,dum = prepC.preparePeptideChain(startindex,endindex,chid, res, resids, resnums, resns, chids, inscodes, pts, mconly, 1, guidedSamplingRadius, caRad, scRad, scMissInds,caMissInds)
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
            bl,bf,nt,rl = prepC.prepareChainTerminal("Nterm", startindex, endindex, firstindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, guidedSamplingRadius, scMissInds,caMissInds)
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
            bl,bf,nt,rl = prepC.prepareChainTerminal("Cterm", startindex, endindex, lastindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, guidedSamplingRadius, scMissInds,caMissInds)
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
            
            #if endindex-startindex <= 5  :   #or random.random() > 1.5 :
            bl,bf,nt,rl = prepC.preparePeptideLoop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, scMissInds, guidedSamplingRadius,None,caMissInds,1)
            #else :
                
            #    bl,bf,nt,rl = prepC.preparePeptideMidloop(startindex, endindex, res, resids, resnums, resns, chids, inscodes, pts, mconly, caRad, scRad, scMissInds, guidedSamplingRadius,caMissInds,1)
                
            #    reorderBuilders = None
        
        else :
            print "Dont understand looptype", looptypes[li] ; assert None
        mergeBRlists(blist, bfoll, numtrials, rlist, bl, bf, nt, rl, dummies, dum)



        modelRenderer = ModelRenderer(res, resns, chids, resnums, inscodes, dummies.keys(), pdbout)
            
    


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
        for an in dummies[k] : aiSC[res[k][an]] = -1

    print "---------------------------------------\n\n\n"        





    radii = [0] * len(pts) ## make radii for all atoms including dummies
    for index, val in res.items() :
        for aname, pi in val.items() :
            if resns[index] in vdwr.keys() : radii[pi] = vdwr[resns[index]][aname]
            else :
                try : radii[pi] = vdwr['XXX'][aname[1]]
                except KeyError : radii[pi] = vdwr['XXX']['C']
    assert len(radii) == len(pts)
    pts = VecVecFloat(pts)
    print "---------------------------------------\n\n\n"

    ssReduction = data.PROBE_DISULFIDE_OVERLAP_MARGIN
    gridHelper = CAtraceGH(VecInt(aiSC), VecInt(ai2pos), scReduction, ssReduction)

    if a != None and b != None and c != None and alpha != None and beta != None and gamma != None : gridHelper.setCellsym(a, b, c, alpha, beta , gamma , sg)
    
    extraRestraints, extraRestrOpt = [], []
    if not type(xrayRestGen) == type([1,2]) : xrayRestGen = [xrayRestGen]


    ## density restraints         
    for rg in xrayRestGen :
        if rg == None : continue
        xrlist, optional = rg.generate(blist, aiSC, res, resids, pts)
        rsize = len(extraRestraints)
        for xr in xrlist : extraRestraints.append(xr)
        if edOpt == 1 : 
            for ori in range(len(xrlist)) :
                extraRestrOpt.append(rsize + ori)
    


    ### Positional and closure restraints
    for xr in rlist :
        extraRestraints.append(xr)

    if allOpt   == 1 :
        rsize = len(extraRestraints)
        for ori in range(len(rlist)) :
            extraRestrOpt.append(rsize + ori)

    ### ligand restraints 
    for xr in irlist :
        extraRestraints.append(xr)

    nbuilt = 0 
    while (nbuilt < nmodels) :
        knownPositions = range(len(pts))
        for b in blist :
            bop = b.getOP()
            for i in range(bop.size()) :
                if bop[i] in knownPositions :
                    knownPositions.remove(bop[i])
    
    
        unbuiltpts= [] ;     built = [] ; unbuilt = [] ;     
    
    
        
        buffRenderer = BufferModelRenderer(pts)

        strategy = PopulationStrategy(backtrack, popsize, pts, radii, knownPositions, gridHelper, blist, numtrials, bfoll, 1 , extraRestraints, extraRestrOpt, buffRenderer, natt, 1, None)
        strategy.ranker = ranker
    
        if strategy.execute() == 0 : ## have to restart, restore knownPositions OR just let the band be same as given
            for b in bl : ## modify knownPositions and pts
                for i in range(b.getOP().size()) :
                    if  b.getOP()[i]  not in missingpts:
                        unbuiltpts.append( b.getOP()[i] )
                
                
        if len(unbuiltpts) == 0 :
            modelRenderer.render(pts,None,missingpointids,unbuiltpts) ;
            nbuilt += 1
        else:
            modelRenderer.render(pts,None,missingpointids,unbuiltpts) ;
            nbuilt += 1

    from copyheader import copyheader
    copyheader(pdbout,pdbfile)

def randomize(doRandomize) :
    import misc, random
    if doRandomize :
        random.seed(doRandomize)
        misc.RanGen.instance().seedme( doRandomize )
    else : 
        random.seed(1973) # for Paul!
        misc.RanGen.instance().seedme(1942) # for Tom!


def callmain() :
    from commonOptions import makeParser , parseOptions
    import optparse ; parser = optparse.OptionParser()
    import sys

    ################# I/O Files and directories ##########################################
    
    parser.add_option("--xyzin", action='store', type='string', dest='pdbfile', help='Name of input coordinates file')
    parser.add_option("--xyzout", action='store', type='string', dest='pdbout', help='Name of output coordinates file',default="modelout.pdb")
    parser.add_option("--dir-xyzout", action='store', type='string', dest='dir_xyzout', help='Name of directory in which the output files will be created. Complete PATH needed')
    parser.add_option("--mapin", action='store', type='string', dest='mapfn', help='Name of input MAP file (complete PATH needed)',default=None)
    parser.add_option("--hklin", action='store', type='string', dest='mtzfn', help='Name of input (phased) MTZ file (complete PATH needed)',default=None)


    ### Assignment of restraints and radii #############################################
    
    parser.add_option("--use-ca-restraints", action='store', dest='caRes', help='[True/False], Apply positional restraints on the C-alpha atoms',default="True")
    parser.add_option("--use-sc-restraints", action='store', dest='scRes',type= 'string', help='[True/False],  Apply positional restraints on the centroid of the sidechain atoms',default="True",)
    parser.add_option("--ca-restraint-radius", action='store', type='float', dest='caRad', help='Radius of spherical restraint ( Angstrom ) on the C-alpha atom position', default=1.0)
    parser.add_option("--sc-centroid-restraint-radius", action='store', type='float', dest='scRad', help='Radius of spherical restraint ( Angstrom ) on the centroid of the sidechain atoms', default=2.0)

    ############ General parameters #####################################################
    
    parser.add_option("--sidechain-vdw-reduction", action='store', type='float', dest='scReduction', help='Factor by which to reduce effective Van der Waals distance for sidechain atoms', default=0.75)
    parser.add_option("--population-size", action='store', type='int', dest='popsize', help='Population size', default=100)
    parser.add_option("--verbose", action='store', type='int', dest='verbose', help='Level of printed output [1-10], 1:Concise output log, 10: Detailed output log', default=6)
    parser.add_option("--backtrack", action='store', type='string', dest='backtrack', help='[numsteps]X[stepsize] e.g. when set to 4X5, In case of failure to build at a residue the program will  backtrack 4 times by 5 residues. For detailed help see Rappertk wiki/webpage', default="4X5")
    parser.add_option("--rotamerlib", action='store', type='string', dest='rotLib', help='[PRL/SCL1.0/SCL0.5/SCL0.2] Name of rotamer library to use when building side chains ', default='PRL')        
    parser.add_option("--num-models", action='store', type='int', dest='nmodels', help='Number of models wanted ', default=1)


    #################### Build parameters ################################################
    
    parser.add_option("--mconly", action='store', type='string', dest='mconly', help='[True/False] Build mainchain only', default="False")
    parser.add_option("--sconly", action='store', type='string', dest='sconly', help='[True/False] Build side chains only, can only be used when MAP/MTZ file is given. See web page for further details', default="False")

    parser.add_option("--opsax", action='store', type='string', dest='opsax', help='[True/False] Reassign side chains with OPSAX, will only be used when MTZ or MAP file is given', default="True")

    parser.add_option("--attempts", action='store', type='int', dest='natt', help='Number of attempts made to build section', default=5)

    parser.add_option("--cacaCutoff", action='store', type='float', dest='cacaCutoff', help='Minimum distance ( angstrom ) between adjacent Calpha atoms in order to detect a chain-break', default=5.)


    #################    Electron density parameters ####################################
    
    parser.add_option("--a", action='store', type='float', dest='a', help='Cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='Cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='Cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='Cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='Cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='Cell spacegroup, in CCP4 notation')
    parser.add_option("--resolution", action='store', type='float', dest='resolution', help='Resolution of the data')
    
    parser.add_option("--FP", action='store', type='string', dest='f1label', help='Column label for FP in MTZ file', default="FP")
    parser.add_option("--SIGFP", action='store', type='string', dest='sigf1label', help='Column label for sigFP in MTZ file', default="SIGFP")
    parser.add_option("--FC", action='store', type='string', dest='f2label', help='Column label for FC in MTZ file', default=None)
    parser.add_option("--PHI", action='store', type='string', dest='phiclabel', help='Column label for PHI in MTZ file', default=None)
    parser.add_option("--use-FREER", action='store', type='string', dest='usefreer', help='[True/False] Use FreeR set ? ', default="False")
    parser.add_option("--FREER", action='store', type='string', dest='freerlabel', help='Column label for FreeR in MTZ file', default=None)




    


    ######### Ouptut parameters #############################################

    parser.add_option("--models-get-native-bfactors", action='store', type='string', dest='nativeBfac', help='[True/False] Assign B-factors of remodelled atoms to original values', default="False")
    parser.add_option("--default-mainchain-b-factor", action='store', type='float', dest='mcBfac', help='The value of B-factor assigned to the newly built main chain atoms', default=20.)
    parser.add_option("--default-sidechain-b-factor", action='store', type='float', dest='scBfac', help='The value of B-factor assigned to the newly built side chain atoms', default=30.)



    ### Electron density parametets #########################################

    parser.add_option("--minimum-sig", action='store', type='float', dest='minXSig', help='Minimum sigma ', default=0.25)
    parser.add_option("--maximum-sig", action='store', type='float', dest='maxXSig', help='Maximum sigma ', default=2.0)

    parser.add_option("--use-ranking", action='store', type='string', dest='userank', help='Use ranking ', default="False")


    ########## Optional restraints ##########################################
    
    parser.add_option("--make-ed-optional", action='store', type='string', dest='edOpt', help='[True/False]  If False, then the mainchain will be unconditionally forced to lie in positive density. If True then positive density restraint on the mainchain will be made optional.This is useful when tracing through a structure with regions in very poor (non-existent) density', default= "False")


    parser.add_option("--make-all-restraints-optional", action='store', type='string', dest='allOpt', help='[True / False ]  If True, then all  restraints will be made optional', default="False")    


    parser.add_option("--maptype", action='store', type='string', dest='mapformat', help='Options are : [2fofc/omit] ', default="2fofc")


    
    
    parser.add_option("--ligfile", action='store', type='string', dest='ligfile', help='ligand description, see msq.ligdesc for description of MSQ in 1di9.pdb')
    parser.add_option("--around-ligand", action='store', type='float', dest='closeCutoff', help='min-dist between residue/ligand to be considered close to ligand', default=10)

    options = parseOptions(parser)
    options.modelN2C  = 'False'
    (options, args) = parser.parse_args()

    from checkfornone import checkfornone_ppi 


    pdbfile,pdbout,dir_xyzout,mapfn,mtzfn,caRes,scRes,caRad,scRad,scReduction,popsize,verbose,backtrack,rotLib,nmodels,mconly,sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freeRlabel,poordum,ccdum,loopdum,startdum,startcdum,stopdum,stopcdum,chaindum,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt,closuredum,addscdum,userotdum =    checkfornone_ppi(options.pdbfile,options.pdbout,options.dir_xyzout,options.mapfn,options.mtzfn,options.caRes,options.scRes,options.caRad,options.scRad,options.scReduction,options.popsize,options.verbose,options.backtrack,options.rotLib,options.nmodels,options.mconly,options.sconly,options.opsax,options.natt,options.cacaCutoff,options.a,options.b,options.c,options.alpha,options.beta,options.gamma,options.sg,options.resolution,options.f1label,options.sigf1label,options.f2label,options.phiclabel,options.usefreer,options.freerlabel,"False",0.9,None,None,None,None,None,None,options.nativeBfac,options.mcBfac,options.scBfac,options.minXSig,options.maxXSig,options.edOpt,options.allOpt,"True","False","False")
    


    print "	--around-ligand ".ljust(30), 
    print "%10f "%options.closeCutoff

    print "	--ligfile ".ljust(30), 
    print "%10s "%options.ligfile



    import misc;    misc.setVerbosity(options.verbose)


    if dir_xyzout == None : print "Directory for output files --dir-xyzout needs to be set. ",dir_xyzout; import sys ; sys.exit()
    if  os.path.isfile(dir_xyzout) : print dir_xyzout,"Is a file, Please rename output directory" ;  import sys ; sys.exit()
    dir_xyzout = fixdirnames(dir_xyzout)


    head,tail  = os.path.split(dir_xyzout)
    if head and os.path.isdir(head):
        print "Checking for directory ............",head
    else:
        print "%s does not exists"%head
        import sys ; sys.exit()
    if not os.path.isdir(dir_xyzout) :
        os.mkdir(dir_xyzout)
    os.chdir(dir_xyzout)
    
    if pdbout == None : print "Output filename needs to be set."; import sys ; sys.exit()    
    if (os.path.isdir(pdbout) == True )  : print "Output Coordinate file %s is a directory"%pdbout; import sys ; sys.exit()
    pdbout = fixdirnames(pdbout)
    head,tail  = os.path.split(pdbout)
    if head != "" : pdbout = tail

    if (os.path.isfile("%s/%s" % (dir_xyzout,"modelinit.pdb"))==True) :
        os.remove("%s/%s" % (dir_xyzout,"modelinit.pdb"))
    if  pdbfile == None : print "Please specify input coordinates file %s " % pdbfile ; import sys ;        sys.exit()
    if (os.path.isdir(pdbfile) == True )  : print "Input Coordinate file %s is a directory"%pdbfile; import sys ; sys.exit()
    pdbfile = fixdirnames(pdbfile)
    if (os.path.isfile(pdbfile)==False)  :        print "Cannot find file %s " %pdbfile ; import sys ;        sys.exit()
    shutil.copyfile(pdbfile, "%s/%s" % (dir_xyzout,"modelinit.pdb"))
    modelIn =  "modelinit.pdb"



    if  options.ligfile == None : print "Please specify ligand description file %s " % options.ligfile ; import sys ;        sys.exit()
    if (os.path.isdir(options.ligfile) == True )  : print "Ligand description file %s is a directory"%options.ligfile; import sys ; sys.exit()
    options.ligfile = fixdirnames(options.ligfile)
    if (os.path.isfile(options.ligfile)==False)  :        print "Cannot find file %s " %options.ligfile ; import sys ;        sys.exit()
    shutil.copyfile(options.ligfile, "%s/%s" % (dir_xyzout,"lig.desc"))
    ligfilepath = "%s/lig.desc" % (options.dir_xyzout)






    opsaxoff = None ; guidedsampling = None;    scvdwr = scReduction ;   
    xrayRestGen = [] ;    ranker = None ; badresids = [] ;         nullres = [] ; missingpts = {}
    misc.setVerbosity(verbose) 
    randomize(1)
    str2bool = {"False": 0 , "True":1}

    mconly =  str2bool[mconly]
    sconly =  str2bool[sconly]
    opsax =  str2bool[opsax]
    usefreer =  str2bool[usefreer]
    nativeBfac =  str2bool[nativeBfac]

    caRes =  str2bool[caRes]
    scRes =  str2bool[scRes]
    edOpt =  str2bool[edOpt]
    allOpt =  str2bool[allOpt]

    
    if mconly == 0 :        mconly = None
    if allOpt == 1 :        edOpt = 1
    if caRes == 0 :         scRes = 0
    

    if caRes == 0 and sconly == 0 : print "--use-ca-restraints to be set to True" ; import sys ; sys.exit()



    if (mapfn !=  None or mtzfn !=None ):
        from stump import getCRYST , getRESO

        if (a == None or b == None or c == None or alpha == None or beta == None or gamma == None) :
            print "Getting cell paramters from coordinate file....."
            a,b,c,alpha , beta , gamma,d1  = getCRYST(pdbfile)
            if (a == None or b == None or c == None or alpha == None or beta == None or gamma == None ):
                print "CRYST card cannot be read from coordinate file. Please input cell paramater a, b , c , alpha, beta , gamma = ",a , b , c , alpha , beta  , gamma 
                import sys ; sys.exit()
            print "Found a b c alpha beta gamma  ", a , b , c , alpha , beta  , gamma 

        if sg == None : 
            print "Getting space group from coordinate file....."
            d1,d2,d3,d4 , d5 , d6, sg  = getCRYST(pdbfile)
            if sg == None : 
                print "Please input space group " , sg ; import sys ; sys.exit()
        ss = ""
        for sg1 in sg:
            if sg1 in ["\n","\t","\s"]:
                continue
            else :
                ss = ss+sg1
        sg = ss
        if sg  in long2shortHM.keys(): shortsg = long2shortHM[sg];  sg = shortsg
        if sg not in sgtable.keys(): print "Check --sg , Not recognised [%s]"%sg ;            import sys ; sys.exit()
        print "Setting Space Group to",sg
        
        if resolution == None : 
            print "Getting resolution limit from coordinate file........"
            resolution = getRESO(pdbfile)
            if (resolution == None):
                print "Please input resolution using --resolution value, currently" , resolution ; import sys ; sys.exit()
            print "Resolution = [ " , resolution, " ] "


    if caRes == 0 and sconly != 1 : print "\n********** WARNING!! No C-alpha positional restraints will be used";caRad = None; 
    if scRes == 0 : print "********** WARNING!! No sidechain centroid positional restraints will be used"; scRad = None

    if  mapfn != None and mtzfn != None :
        print "********** WARNING!! both mtz and map file given, please only give either MTZ file or MAP file  *************"
        import sys ; sys.exit()

    if  mapfn != None:
        n = 2
        if (os.path.isfile("%s/%s" % (dir_xyzout,"init.map"))==True) :
            os.remove("%s/%s" % (dir_xyzout,"init.map"))        
        useomitmap,usecnsmap,useccp4map  = None,None,None
        esmin, esmax, esmean, rcmult = minXSig, maxXSig, .0, 5

        mapfilepath= mapfn
        if (os.path.isfile(mapfilepath)==False) :
            print "Cannot find file %s " %mapfn ; import sys ;            sys.exit()
        if  "map" not in mapfn : 
            print "Cannot understand i/p map format,the file should be  (*.map) format " ; import sys ; sys.exit()
        from xray import mapdump
        if  mapdump(mapfn)  == None :
            print "Check format of map file, needs to be CCP4 format";
            import sys ;              sys.exit()            


        shutil.copyfile(mapfilepath, "%s/init.map" % (dir_xyzout))
        ccp4map = mapfn

    elif mtzfn != None :
        n = 2
        if (os.path.isfile("%s/%s" % (dir_xyzout,"init.mtz"))==True) :
            os.remove("%s/%s" % (dir_xyzout,"init.mtz"))            
        mtzfilepath= "%s"%(mtzfn)
        esmin, esmax, esmean, rcmult  = minXSig, maxXSig, .0, 5
        
        if (os.path.isfile(mtzfilepath)==False) :
            print "Cannot find file %s " %mtzfn ; import sys ;sys.exit()
        if  "mtz" not in mtzfn : 
            print "Cannot understand i/p HKLIN format,the file should be in mtz (*.mtz) format " ; import sys ; sys.exit()


        shutil.copyfile(mtzfilepath, "%s/init.mtz" % (dir_xyzout))
        mtzfilepath = "%s/init.mtz" % (dir_xyzout)


        ################   COPY MtZ COLOUMN LABELS  #########################################
        from xray import mtzdump

        if (f1label == None):
            print "Please specify FP label  " , f1label ; import sys ; sys.exit()
        if (sigf1label == None):
            print "Please specify SIGFP label  " , sigf1label ; import sys ; sys.exit()
            
        if usefreer == 1 : 
            if freeRlabel == None :
                print "If you would like to use the freeR set, please specify column name for freeR"
                import sys ; sys.exit()
            
        else :
            freeRlabel = None


        if  mtzdump(mtzfn,f1label,sigf1label,freeRlabel,f2label,phiclabel)  == None :
            print "Check structure factor file, format needs to be MTZ"
            import sys ;              sys.exit()            
            

        if (f2label == None or phiclabel == None ):
            print "Column labels for  FC and PHI  not specified, will use input coordinate structure to obtain FC and PHI"

            
            if sfall(pdbfile, mtzfn, "phased.mtz",resolution,usefreer,f1label,sigf1label,freeRlabel) == None :
                
                print "Structure factor file cannot be phased , please enter column labels for  FC and PHI"
                import sys ; sys.exit()
            
            mtzfn  = "phased.mtz"
            f2label  = "FC" ; phiclabel = "PHIC"



        ### MAKE UNWEIGHTED DIFF MAP : todo sigmaA weighted maps #########################

        from xray import fftnew
        useomitmap,usecnsmap,useccp4map  = None,None,None
        mapformat  = "2fofc"
        if mapformat == "cns": usecnsmap = 1
        elif mapformat == "omit": useomitmap = 1
        elif mapformat == "2fofc": useccp4map = 1
        else : print "Unrecognised mapformat ,exiting" ; import sys; sys.exit()
            
        if useccp4map == 1 or useomitmap   == 1 :
            fofcmap = "%s%dFo-%dFC.map"%(mtzfn,2,(1))
            fcmap = "%sFC.map"%mtzfn
            mapfn = fofcmap

            if usefreer == 1 : 
                fftnew(mtzfn, fofcmap, pdbfile,f1label,f2label,phiclabel,sigf1label,n,(n)-1,freeRlabel)
                fftnew(mtzfn, fcmap, pdbfile ,f1label,f2label,phiclabel,sigf1label,0,1,freeRlabel,1)
            else :
                fftnew(mtzfn, fofcmap, pdbfile , f1label,f2label,phiclabel,sigf1label,n,(n)-1,None)
                fftnew(mtzfn, fcmap, pdbfile ,f1label,f2label,phiclabel,sigf1label,0,1,freeRlabel,1)



        if useomitmap == 1 :
            omit_mapfile = "%s%dFo-%dFC.omit.map"%(mtzfn,n,((n)-1))
            mapfn = omit_mapfile 
            if usefreer == 1 : 
                omitsucc = omitmap(mtzfn, omit_mapfile, f1label,f2label,phiclabel,sigf1label,n,(n)-1,freeRlabel)
            else :
                omitsucc = omitmap(mtzfn, omit_mapfile, f1label,f2label,phiclabel,sigf1label,n,(n)-1)
            if omitsucc == None :
                useomitmap = 0
                useccp4map = 1
                mapfn = fofcmap
        if usecnsmap == 1 :
            fofcmap = mtzfn+"2fofc.map"
            fcmap = mtzfn+"fc.map"
            mapfn = fofcmap
            
        
        

    ############### SET UP   ######################
    if (mapfn !=  None or mtzfn !=None ):
        esmin, esmax, esmean, rcmult = minXSig, maxXSig, .0, 5
        mapcoeff = "%dF1-%dF2"%(n, n-1)
        xrayRestGen = []
        xrayRestGen.append( prepareChainV5.XrayRestraintsGenerator(mapfn, "map", f2label, phiclabel, mapcoeff, esmin, esmax, esmean, [], edOpt ) )

        if options.userank not in ["True", "False"] :
            print "Incorrect value entered for --use-rank. Options are [True/False]"
        if options.userank == "True" :
            ranker = XrayRanker(mapfn, "map", f2label, f1label, mapcoeff, esmin, esmax)
            ranker.rankChildren = rcmult ; ranker.rankRecurse = 1
            ranker.rankLeaderBuilderOnly = None ; ranker.rankGivenEnsemble = None
        

    #os.chdir(options.dir_xyzout)

    if (options.mconly == 1 and  options.sconly == 1) :
        print "Mainchain only and sidechain only modes cannot be used together"
        import sys ; sys.exit()
        
        

    main(pdbfile,pdbout,dir_xyzout,mapfn,mtzfn,caRes,scRes,caRad,scRad,scReduction,popsize,verbose,backtrack,rotLib,nmodels,mconly,sconly,opsax,natt,cacaCutoff,a,b,c,alpha,beta,gamma,sg,resolution,f1label,sigf1label,f2label,phiclabel,usefreer,freeRlabel,nativeBfac,mcBfac,scBfac,minXSig,maxXSig,edOpt,allOpt,mapformat,ranker, xrayRestGen,ligfilepath, options.closeCutoff)


#    if options.mtzfn != None : 
#
#        xrayRestGen = XrayRestraintsGenerator(options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype, options.sigXmin, options.sigXmax, options.sigXmean)
##        esmin, esmax, esmean, rcmult = .000, 5., .0, 5 

#        ranker = XrayRanker("phased.mtzFPFCPHIC_A.map", "map", "FC", "PHIC", "2F1-F2", esmin, esmax)
#        ranker.rankChildren = rcmult ; ranker.rankRecurse = 1
#        ranker.rankLeaderBuilderOnly = None ; ranker.rankGivenEnsemble = None
        


if __name__ == "__main__" : callmain()
