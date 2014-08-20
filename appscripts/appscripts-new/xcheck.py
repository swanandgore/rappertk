import math, os, sys, string, re
from restraints import EDrestraint
from builders import VecInt, VecVecFloat
from pdbr import protein, isAAres
import prot2res
import data, misc
from data import vdwr
from geometry import CAtraceGH
from stdev import findMeanStdev

def main(pdbfile, mtzfn, f1label, f2label, philabel, maptype, aasc=None) :
    mt = {"2F1-F2":0, "F1":1}
    prot = protein(pdbfile, read_hydrogens=0, read_waters=1, read_hets=1)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    pts = VecVecFloat(pts)
    Xcorrs = {}
#    print len(res.keys()), len(set(res.keys()))
    allkeys = list(res.keys()) ; allkeys.sort()
    for ri in allkeys :
        key = resids[ri]
        if aasc and (not isAAres(resns[ri]) or resns[ri] in ["GLY","ALA"]) : continue
        if aasc :
            ptinds = []
            for k,v in res[ri].items() :
                if not k in [" N  "," CA "," C  "," O  "," CB "] : ptinds.append(v)
            ptinds = VecInt(ptinds)
        else :
            ptinds = VecInt( res[ri].values() )

        xcorr = EDrestraint.findBadfit(mtzfn, f1label, f2label, philabel, pts, ptinds)
        Xcorrs[ key ] = xcorr
#        print ri, "RESID", key        print "XcorrZ [%s]"%key, xcorr
    #mean, stdev = findMeanStdev(Xcorrs.values())
    #for key,val in Xcorrs.items() : print "XcorrZ", key, val, (val-mean) / stdev
    return Xcorrs

class XrayRanker :
    def __init__(s, phasedmtz, fplabel, fclabel, philabel, maptype, esmin, esmax) :
        s.mt = {"2F1-F2":0, "F1":1}
        ## Change to make universal maptype
        s.phasedmtz, s.fplabel, s.fclabel, s.philabel, s.maptype, s.esmin, s.esmax = phasedmtz, fplabel, fclabel, philabel, maptype, esmin, esmax
        s.rankChildren, s.rankRecurse, s.rankGivenEnsemble = 1, None, None
        
    def score(s, crdlist, pis=None) :
        crdlist1 = VecVecFloat(crdlist)
        #print "-----------------------------------"
        #for ci in range(crdlist1.size()) :
        #    print crdlist1[ci][0] , crdlist1[ci][0] , crdlist1[ci][0]
        if not pis : pis = VecInt( range(crdlist1.size()) )
        else : pis = VecInt( range(pis) )
        if re.compile("\.map$").search(s.phasedmtz) :
            er = EDrestraint.makeEDrestraintFromMap(pis, "Xranker", s.phasedmtz, s.esmin, s.esmax, s.esmax)
        else :
            er = EDrestraint.makeEDrestraintFromMTZ(pis, "Xranker", s.phasedmtz, s.fplabel, s.fclabel, s.philabel, s.maptype, s.esmin, s.esmax, s.esmax)
        return -1*er.scoreAround(crdlist1)


class XrayScorer :
    def __init__(s, plotfn, cutoff) :
        s.mt = {"2F1-F2":0, "F1":1}
        s.plotfn = plotfn
        s.cutoff = cutoff
        s.scores, s.keyorder = {}, []
    def tplot(s) :
        for ki in range(len(s.keyorder)) :
            key, val = s.keyorder[ki], s.scores[s.keyorder[ki]]
            print "[%s]" % key,
            newval = list(val) ; newval.reverse()
            for nv in newval :
                val = int(round(nv * 10))
                for i in range(val) : print "_",
                for i in range(10-val) : print " ",
            print "|"
        return
    def gplot(s) :
        if not s.plotfn : return
        ofp = open(s.plotfn+".dat", 'w')
        for key in s.keyorder :
            for val in s.scores[key] : print >> ofp, val,
            print >> ofp, ""
        ofp.close()
        cmd, inputs = "gnuplot", ['set term postscript color enhanced  "Helvetica" 4']
        inputs.append('set out "' + s.plotfn + '"')
        xt = "set xtics ("
        for ki in range(len(s.keyorder)) : xt += ' "%s" %d,' % (s.keyorder[ki], ki)
        inputs.append( xt[0:len(xt)-1] + ")")
        inputs.append( "set xtics rotate by 90" )
        inputs.append("set size ratio 0.1") ; inputs.append("set grid")
        pl = "plot"
        for i in range(len(s.scores.values()[0])) : pl += " '" + s.plotfn + ".dat' u %d w l," % (i+1)
        inputs.append( pl[0:len(pl)-1] )
        for ip in inputs : print ip
        from procrun import proc_run_exitOnError as execCmd
#        ]        execCmd(cmd, inputs)

    def score(s, pdbfile, mtzfn, f1label="FP", f2label="FC", philabel="PHIC", maptype="2F1-F2", aasc="mcsc",verbose=5) :
        assert aasc in ["pept","mc","sc","mcsc"]
        prot = protein(pdbfile, read_hydrogens=0, read_waters=1, read_hets=1)
        res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
        pts = VecVecFloat(pts)
        retkeys = []
        badlines = []
        retids = []
        allkeys = list(res.keys()) ; allkeys.sort()

        for ai in range(len(allkeys)) :
            ri = allkeys[ai]
            key = resids[ri]
            if not isAAres(resns[ri]) : continue
            if aasc=="sc" and resns[ri] in ["GLY","ALA"] : continue

            if aasc=="sc" and isAAres(resns[ri]) :
                ptinds = []
                for k,v in res[ri].items() :
                    if not k in [" N  "," CA "," C  "," O  "," CB "] : ptinds.append(v)
                ptinds = VecInt(ptinds)

            elif aasc=="pept" and isAAres(resns[ri]) and ai+1 < len(allkeys) and isAAres(resns[allkeys[ai+1]]) :
                ptinds = []
                for k,v in res[ri].items() :
                    if k in [" C  "," O  "] : ptinds.append(v)
                for k,v in res[ allkeys[ai+1] ].items() :
                    if k in [" N  ",] : ptinds.append(v)
                ptinds = VecInt(ptinds)
            elif aasc=="mc" and isAAres(resns[ri]) :
                ptinds = []
                for k,v in res[ri].items() :
                    if k in [" N  "," CA "," C  "," O  "," CB "] : ptinds.append(v)
                ptinds = VecInt(ptinds)

            elif aasc=="mcsc" and isAAres(resns[ri]) :
                ptinds = VecInt( res[ri].values() )

            else:
                ptinds = VecInt( res[ri].values() )
                
            if re.compile("\.map$").search(mtzfn) and re.compile("\.map$").search(f1label) :
                xcorr = EDrestraint.findBadfit1(mtzfn, f1label, pts, ptinds)
            else :
                xcorr = EDrestraint.findBadfit(mtzfn, f1label, f2label, philabel, pts, ptinds)
#            print "XcorrZ %s [%s] [%f]"%(aasc,key, xcorr)
            if not s.scores.has_key(key) : s.scores[key] = [] ; s.keyorder.append(key)
            s.scores[key].append(xcorr)
            if verbose > 5 :
                print "[%s], [%s] correlation coefficient = %f"%(key,aasc,xcorr)
            if xcorr < s.cutoff :
                retkeys.append(key) ;
                bl = "[%s]" % key
                #for xc in s.scores[key] :
                #    if xc < s.cutoff : bl += " ? "
                #    else : bl += " ! "
                badlines.append(bl)
                retids.append(ri)
#        print len(badlines), "BADRESIDS", aasc, len(s.scores.values()[0]), s.cutoff, "------------------------------------" ;
#        for bl in badlines : print bl
        s.gplot() #s.tplot()
        return retkeys,retids

def callmain() :
    from commonOptions import makeParser, parseOptions, addXrayOptions
    parser = makeParser()
    parser = addXrayOptions(parser)
    options = parseOptions(parser)

    main(options.pdbfile, options.mtzfn, options.f1label, options.f2label, options.philabel, options.maptype)

if __name__ == "__main__" : callmain()



        #if f1label == "map" : xed = EDrestraint.makeEDrestraintFromMap(ptinds, "xcheckEDrestraint", mtzfn, 0.5, 5., 2.)
        #else : xed = EDrestraint.makeEDrestraintFromMTZ(ptinds, "xcheckEDrestraint", mtzfn, f1label, f2label, philabel, mt[maptype], 0.5, 5, 2)
        #xcorr = xed.scoreSigma(pts) / ptinds.size()
