import string

class ProtData : pass

def readprot(lines) :
    from pdbr import protein
    import prot2res
    prot = protein(lines, read_hydrogens=None, read_waters=None, read_hets=None)
    res, resids, resnums, resns, chids, inscodes, pts = prot2res.readProtRes(prot)
    resid2ri = {}
    for ri in resids.keys() : resid2ri[resids[ri]] = ri
    pd = ProtData()
    pd.res, pd.resids, pd.resnums, pd.resns, pd.chids, pd.inscodes, pd.pts, pd.resid2ri \
            = res, resids, resnums, resns, chids, inscodes, pts, resid2ri
    return pd
    
def readfile(pdbfile) :
    lines = open(pdbfile,'r').readlines()
    newlines, prots = [], []
    for l in lines :
        if l[0:6] == "ENDMDL" :
            prots.append( readprot(newlines) )
            newlines = []
        elif l[0:3] == "END" : continue
        elif l[0:5] == "MODEL" : continue
        else : newlines.append(l)
    if len(newlines) > 0 : prots.append( readprot(newlines) )
    return prots

def findMCrmsd(pdi, pdj, loopres=None) :
    rmsd, cnt = 0., 0
    for resid in pdi.resids.values() :
        if loopres != None and not resid in loopres : continue
        rii, rij = pdi.resid2ri[resid], pdj.resid2ri[resid]
        #for an in [" N  "," CA "," C  "," O  ",] :
        for an in [" CA ",] :
            p1 = pdi.pts[ pdi.res[rii][an] ]
            p2 = pdj.pts[ pdj.res[rij][an] ]
            rmsd += (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2])
            cnt += 1
    import math
    return math.sqrt( rmsd/cnt )

def findSCrmsd(pdi, pdj, rii, rij) :
    assert set(pdi.res[rii].keys()) == set(pdj.res[rij].keys())
    rmsd, cnt = 0., 0
    for an in pdi.res[rii].keys() :
        if an in [" N  "," CA "," C  "," O  ",] : continue
        p1 = pdi.pts[ pdi.res[rii][an] ]
        p2 = pdj.pts[ pdj.res[rij][an] ]
        rmsd += (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2])
        cnt += 1
    import math
    return math.sqrt( rmsd/cnt )

if __name__=="__main__":
    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--e1", action='store', type='string', dest='e1', help='')
    parser.add_option("--e2", action='store', type='string', dest='e2', help='')
    parser.add_option("--loopres", action='store', type='string', dest='loopres', help='', default=None)

    (options, args) = parser.parse_args()

    options.e1 = string.split( options.e1, ':' )
    options.e2 = string.split( options.e2, ':' )

    print "set 1 consists of", options.e1
    print "set 2 consists of", options.e2

    if options.loopres != None :
        loopres = []
        for l in open(options.loopres, 'r').readlines() : loopres.append( l[1:len(l)-2] )
        options.loopres = loopres

    prots1, prots2 = [], []
    for fn in options.e1 : prots1 += readfile(fn) ## the underlying heterogeneity
    for fn in options.e2 : prots2 += readfile(fn) ## the model of heterogeneity

    for resid in prots1[0].resids.values() :
        for pd in prots1[1:] : assert resid in pd.resids.values()
    for resid in prots2[0].resids.values() :
        for pd in prots2[1:] : assert resid in pd.resids.values()
    for resid in prots2[0].resids.values() :
        assert resid in prots1[0].resids.values()


    for resid in prots1[0].resids.values() :
        if options.loopres != None and not resid in options.loopres : continue
        hetmod, onemod = 0., 999.
        for pdi in prots1 :
            hm = 999.
            for pdj in prots2 :
                rii, rij = pdi.resid2ri[resid], pdj.resid2ri[resid]
                if pdi.resns[rii] in ["GLY","ALA"] : rmsd = 0.
                else : rmsd = findSCrmsd(pdi, pdj, rii, rij)
                if onemod > rmsd : onemod = rmsd
                if hm > rmsd : hm = rmsd
            hetmod += hm
        print resid, "%6.3f %6.3f" % (hetmod, onemod)

    hetmodMC = 0.
    for pdi in prots1 :
        hmMC = 999.
        for pdj in prots2 :
            rmsd = findMCrmsd(pdi,pdj, options.loopres)
            print rmsd
            if hmMC > rmsd : hmMC = rmsd
        hetmodMC += hmMC
    print "HETMOD-MC", hetmodMC
