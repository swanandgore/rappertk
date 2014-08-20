import string, os, re
from xray import cif2mtz, uniqueify, sfall, mtz2hkl, cns_generate, cns_anneal, sgCCP4toCNS, mapman


def fixCNSop(pdbfile) :
    from pdbr import isPdbAAline, isPdbAtomLine, line2atomname, line2resn, changeResn, changeAtomname
    lines = []
    for l in open(pdbfile, 'r').readlines() :
        if not isPdbAtomLine(l) : lines.append(l)
        elif line2resn(l) == "CU " : lines.append( changeResn(l, " CU") )
        elif line2resn(l) == "MSE" and line2atomname(l) == " SE " : lines.append( changeAtomname(l, "SE  ") )
        elif line2resn(l) == "ILE" and line2atomname(l) == " CD " : lines.append(re.sub(" CD ", " CD1", l))
        elif line2atomname(l) == " OT1" : lines.append(re.sub(" OT1", " O  ", l))
        elif line2atomname(l) == " OXT" : continue
        else : lines.append(l)
    op = open(pdbfile, 'w')
    for l in lines : op.write(l)
    op.close()
def cnsRefinement(mtzin,pdbin, mtzout,pdbout, a,b,c,alpha,beta,gamma,sg,reso, cnsArgs,cycle, extraTOPfile=None, extraPARfile=None) :
    mtz2hkl(mtzin, "cns.hkl")
    cns_generate(pdbin, "generate.mtf", "generate.pdb", extraTOPfile, extraPARfile, "generate.log")
    wa = -1 ;
    cns_anneal(a, b, c, alpha, beta, gamma, sgCCP4toCNS[sg], reso,
        "cns.hkl", "generate.mtf", "generate.pdb", extraPARfile, "anneal%d.log"%cycle, wa, cnsArgs[cycle]["num_cycles"], cnsArgs[cycle]["temperature"])
    fixCNSop("anneal.pdb")
    os.rename("anneal.pdb", pdbout)
    sfall(pdbout, "rfree.mtz", mtzout, reso)
    mapman("anneal_2fofc.map", mtzout+"2fofc.map")
    mapman("anneal_fc.map", mtzout+"fc.map")
    #moleman(pdbout)

def findRfreeCNSfile(cnsfn) :
    for l in open(cnsfn,'r').readlines() :
        if re.compile("final.*free").search(l) :
            return string.atof( re.sub(".*free_r=", "", l) )
    assert None

if __name__=="__main__" :
    import optparse ; parser = optparse.OptionParser()
    parser.add_option("--collection", action='store', type='string', dest='collection', help='string of single-conformer files separated by :')
    parser.add_option("--sf", action='store', type='string', dest='sf', help='deposited str factors')
    parser.add_option("--scratchdir", action='store', type='string', dest='scratchdir', help='', default=None)
    parser.add_option("--a", action='store', type='float', dest='a', help='cell dimension a')
    parser.add_option("--b", action='store', type='float', dest='b', help='cell dimension b')
    parser.add_option("--c", action='store', type='float', dest='c', help='cell dimension c')
    parser.add_option("--alpha", action='store', type='float', dest='alpha', help='cell angle alpha')
    parser.add_option("--beta", action='store', type='float', dest='beta', help='cell angle beta')
    parser.add_option("--gamma", action='store', type='float', dest='gamma', help='cell angle gamma')
    parser.add_option("--sg", action='store', type='string', dest='sg', help='cell spacegroup, in CCP4 notation')
    parser.add_option("--resolution", action='store', type='float', dest='resolution', help='resolution of the data')

    (options, args) = parser.parse_args()

    collection = string.split( options.collection, ':' )

    if re.compile("\.pdb$").search(collection[0]) : pass ## meaning that all pdbs are given , else assuming that directories are given
    else :
        newcoll = []
        for adir in collection :
            minrfree = 999. ; minpdb = None
            for fn in os.listdir(adir) :
                if not re.compile("^cns.*pdb$").search(fn) : continue
                if re.compile("^cns0.pdb$").search(fn) : continue
                if re.compile("^cns1.pdb$").search(fn) : continue
                rfree = findRfreeCNSfile(adir+"/"+fn)
                if minrfree > rfree : minrfree = rfree ; minpdb = fn
            assert minpdb != None
            newcoll.append(adir + "/" + minpdb)
        collection = newcoll
    print collection

    rfrees = []
    for coll in collection : rfrees.append( findRfreeCNSfile(coll) )
    for ri in range(len(rfrees)) : ## order collection based on Rfree
        for rj in range(ri+1,len(rfrees)) :
            if rfrees[ri] > rfrees[rj] :
                temp = rfrees[ri] ; tempfn = collection[ri]
                rfrees[ri] = rfrees[rj] ; collection[ri] = collection[rj]
                rfrees[rj] = temp ; tempfn = collection[rj] = tempfn
    print collection


    if not os.path.isdir(options.scratchdir) : os.mkdir(options.scratchdir)
    os.chdir(options.scratchdir)


    ## find intra-collection ca-rmsd, aa-rmsd, chi-1, chi-12 variations
    carmsd, mcrmsd, aarmsd, chi1cons, chi12cons, count = 0.,0.,0.,0.,0., 0
    from pdbr import protein
    from evalCAtrace import comparePhiPsiOmegaChi
    for ci in range(len(collection)) :
        for cj in range(ci+1, len(collection)) :
            prot1 = protein( collection[ci], read_hydrogens=0, read_waters=0, read_hets=0 )
            prot2 = protein( collection[cj], read_hydrogens=0, read_waters=0, read_hets=0 )
            results = comparePhiPsiOmegaChi( prot1, prot2 )
            chi1cons += results[0]
            chi12cons += results[2]
            carmsd += results[4]
            mcrmsd += results[5]
            aarmsd += results[6]
            count += 1
    print "INTRA-COMPARISONS %6.3f %6.3f %6.3f %6.3f %6.3f" % ( carmsd/count, mcrmsd/count, aarmsd/count, chi1cons/count, chi12cons/count )

    cif2mtz(options.sf, "base.mtz", options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg)
    uniqueify("base.mtz", "rfree.mtz")

    for colsize in range(2, len(collection)+1) :
        from multProtref import joinPDBs
        joinPDBs("coll.pdb", collection[0:colsize])

        sfall("coll.pdb", "rfree.mtz", "phased.mtz")

        cnsArgs = {}
        for cycle in range(20) : cnsArgs[cycle] = {} ; cnsArgs[cycle]["num_cycles"] = 1 ; cnsArgs[cycle]["temperature"] = 50

        cnsRefinement("phased.mtz", "coll.pdb", "phasedout.mtz", "collout.pdb",
            options.a, options.b, options.c, options.alpha, options.beta, options.gamma, options.sg, options.resolution,
            cnsArgs, 0)

        os.rename("collout.pdb", "coll%d.pdb" % colsize)
        os.rename("phasedout.mtz2fofc.map", "phasedout%d.map" % colsize)

    for colsize in range(2, len(collection)+1) :
        print "MINRFREE", colsize, findRfreeCNSfile("coll%d.pdb" % colsize)
